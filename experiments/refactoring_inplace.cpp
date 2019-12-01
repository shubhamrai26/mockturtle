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

#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/refactoring_inplace.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <lorina/aiger.hpp>

#include <experiments.hpp>

template<typename Ntk>
class cut_compute
{
public:
  using node = typename Ntk::node;

public:
  explicit cut_compute()
  {
  }

  std::vector<std::vector<node>> operator()( node const& n )
  {
    (void)n;
    return {{}};
  }
}; /* cut_compute */

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, float, bool> exp( "cut_rewriting", "benchmark", "size_before", "size_after", "runtime", "equivalent" );

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
      continue;

    xag_npn_resynthesis<aig_network> resyn;
    cut_compute<aig_network> cut_comp;

    refactoring_inplace_params ps;
    ps.progress = true;

    uint32_t size_before = aig.num_gates();

    refactoring_inplace_stats st;
    refactoring_inplace( aig, resyn, cut_comp, ps, &st );
    aig = cleanup_dangling( aig );

    auto cec = abc_cec( aig, benchmark );

    exp( benchmark, size_before, aig.num_gates(), to_seconds( st.time_total ), cec );
  }

  exp.save();
  exp.compare( {}, {}, {"size_after"});

  return 0;
}
