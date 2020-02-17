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

#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/verilog.hpp>
#include <mockturtle/algorithms/xmg_resub.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/resubstitution.hpp>
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/algorithms/node_resynthesis/xmg3_npn.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/xmg.hpp>
#include <mockturtle/properties/xmgcost.hpp>
#include <mockturtle/io/write_verilog.hpp>


#include <experiments.hpp>

int main()
{
    using namespace experiments;
    using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, float, float, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool > exp( "xmg_resubstitution", "benchmark", "size_before_cr", "size_before_resub", "size_after", "runtime_resub", "runtime_rewrite", "total_xor3", "actual_xor3", "actual_xor2", "total_maj", "actual_maj", "remaining_maj","iteration #", "improv_rw", "improv_resub", "eq_rw", "eq_resub" );

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    if (benchmark != "adder" && benchmark != "div") 
        continue;
    fmt::print( "[i] processing {}\n", benchmark );
    xmg_network xmg;
    lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xmg ) );

    xmg_cost_params xmg_ps2;
    uint32_t size_before_resub;
    uint32_t size_before_cr;
    uint32_t num_iters = 0;

    // XMG resubstitution parameter set
    resubstitution_params resub_ps;
    resubstitution_stats resub_st;
    resub_ps.max_pis = 8u;
    resub_ps.max_inserts = 1u;  // Discuss with Heinz once.

    // XMG rewriting parameter set
    cut_rewriting_params cr_ps;
    cut_rewriting_stats cr_st;
    cr_ps.cut_enumeration_ps.cut_size = 4;
    cr_ps.progress = true;


    do 
    {
        num_iters++;
        size_before_cr = xmg.num_gates();
	xmg_ps2.reset();

        xmg3_npn_resynthesis<xmg_network> resyn;
        cut_rewriting( xmg, resyn, cr_ps, &cr_st );
        xmg = cleanup_dangling( xmg );
        //num_gate_profile(xmg,xmg_ps2);
        
        const auto cec2 = benchmark == "hyp" ? true : abc_cec( xmg, benchmark );
        auto improvement_after_rewriting = size_before_cr - xmg.num_gates();
        size_before_resub = xmg.num_gates();

        xmg_resubstitution(xmg, resub_ps, &resub_st);
        xmg = cleanup_dangling( xmg );
        num_gate_profile(xmg,xmg_ps2);
        xmg_ps2.report();
    
        const auto cec = benchmark == "hyp" ? true : abc_cec( xmg, benchmark );
        auto improvement_after_resub = size_before_resub - xmg.num_gates();
        std::cout <<  "For benchmark "<< benchmark << " improvement after rewrite" << improvement_after_rewriting << " and improvement after resub " << improvement_after_resub << " at iteration # " << num_iters << std::endl;

        exp( benchmark, size_before_cr, size_before_resub, xmg.num_gates(), to_seconds( resub_st.time_total ), to_seconds( cr_st.time_total), xmg_ps2.total_xor3, xmg_ps2.actual_xor3, xmg_ps2.actual_xor2, xmg_ps2.total_maj, xmg_ps2.actual_maj, xmg_ps2.remaining_maj, num_iters, improvement_after_rewriting, improvement_after_resub, cec2, cec );

    } while ( (size_before_cr - xmg.num_gates()) > 0 );

    mockturtle::write_verilog(xmg,benchmark+".v");
  }
  
  exp.save();
  exp.table();
  return 0;
}
