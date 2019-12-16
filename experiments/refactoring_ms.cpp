/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2019  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <mockturtle/algorithms/extract_subnetwork.hpp>
#include <mockturtle/algorithms/node_resynthesis/composed.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
#include <mockturtle/algorithms/node_resynthesis/dsd.hpp>
#include <mockturtle/algorithms/refactoring_inplace.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/views/fanout_view.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <lorina/aiger.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, double, float, bool> exp( "refactoring", "benchmark", "size_before", "size_after", "diff", "diff[%]", "runtime", "equivalent" );

  /* refactoring parameters */
  refactoring_inplace_params ps;
  ps.progress = true;
  ps.max_pis = 10;

  dsd_resynthesis_params dsd_ps;
  dsd_ps.dsd_ps.with_xor = false;
  dsd_ps.prime_input_limit = 6u;
  auto cexact_resyn = cached_exact_xag_resynthesis<aig_network>( "/tmp/cache_exact.json", 10e5, *dsd_ps.prime_input_limit );
  dsd_resynthesis<aig_network, decltype( cexact_resyn )> dsd_resyn( cexact_resyn, dsd_ps );
  cached_resynthesis<aig_network, decltype( dsd_resyn )> cdsd_resyn( dsd_resyn, ps.max_pis, "/tmp/cache_dsd.json" );

  for ( auto const& benchmark : epfl_benchmarks( ~experiments::hyp ) )
  {
    using aig_view_t = fanout_view2<depth_view<aig_network>>;

    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      fmt::print( "[i] could not read benchmark {}\n", benchmark );
      continue;
    }

    uint32_t const size_before = aig.num_gates();

    refactoring_inplace_stats st;
    depth_view<aig_network> depth_aig{aig};
    aig_view_t aig_view{depth_aig};

    /* cut computing function */
    xcut<aig_view_t> cut_comp( aig_view, ps.max_pis );
    refactoring_inplace( aig_view, cut_comp, cdsd_resyn, ps, &st );
    aig = cleanup_dangling( aig ); // invalidates aig_view

    uint32_t const size_after = aig.num_gates();

    auto const cec = ( benchmark == "hyp" ) ? true : abc_cec( aig, benchmark );
    exp( benchmark,
         size_before,
         aig.num_gates(),
         size_before - size_after,
         100.0*(size_before - size_after) / size_before,
         to_seconds( st.time_total ),
         cec );

    std::cout << "cec = " << cec << std::endl;

    cexact_resyn.save();
    cdsd_resyn.save();
  }

  exp.save();
  exp.table();

  cexact_resyn.report();
  cdsd_resyn.report();

  return 0;
}
