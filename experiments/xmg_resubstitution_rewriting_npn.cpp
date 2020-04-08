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
#include <mockturtle/algorithms/xmg_optimization.hpp>
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

  experiment<std::string, uint32_t, uint32_t, uint32_t, float, float, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool > exp3( "xmg_resubstitution", "benchmark", "size_before_cr", "size_before_resub", "size_after", "runtime_resub", "runtime_rewrite", "total_xor3", "actual_xor3", "actual_xor2", "total_maj", "actual_maj", "remaining_maj","iteration #", "improv_rw", "improv_resub", "eq_rw", "eq_resub" );
  
  experiment<std::string, uint32_t, float, float, float, float, uint32_t, uint32_t, uint32_t, uint32_t, float, float, bool> exp( "xmg_resubstituion", "benchmark", "iter.", "imp_rw", "imp_rs", "time_rw", "time_rs", "xor3", "xor3'", "maj", "maj'", "xor3_imp", "maj_imp", "equivalent" );

  experiment<std::string, float> exp2 ( "xmg_resubstituion", "benchmark", "area_imp" );
  for ( auto const& benchmark : epfl_benchmarks() )
  {
    //if (benchmark != "adder" && benchmark != "div") 
    //    continue;
    fmt::print( "[i] processing {}\n", benchmark );
    xmg_network xmg;
    lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xmg ) );
    std::string const genlib_path = "/afs/pd.inf.tu-dresden.de/users/shubham/abc-vlsi-cad-flow/std_libs//date_lib_count_tt_4.genlib";

    float area_before = abc_map( xmg, genlib_path );

    xmg_cost_params xmg_ps, xmg_ps2;
    uint32_t size_before_resub;
    uint32_t size_before_cr;
    uint32_t num_iters = 0;
    float improv_after_rewrite = 0;
    float total_imp = 0;
    float rel_xor3 = 0;
    float rel_maj = 0;

    // XMG resubstitution parameter set
    resubstitution_params resub_ps;
    resubstitution_stats resub_st;
    resub_ps.max_pis = 8u;
    //resub_ps.progress = true;
    resub_ps.max_inserts = 1u;  // Discuss with Heinz once.
    resub_ps.use_dont_cares = true;  // Discuss with Heinz once.
    resub_ps.window_size = 16u;  // Discuss with Heinz once.

    // XMG rewriting parameter set
    cut_rewriting_params cr_ps;
    cut_rewriting_stats cr_st;
    cr_ps.cut_enumeration_ps.cut_size = 4;
    //cr_ps.progress = true;


    do 
    {
        num_iters++;
        size_before_cr = xmg.num_gates();
        xmg_ps2.reset();
        xmg_ps.reset();
        num_gate_profile(xmg,xmg_ps);

        xmg3_npn_resynthesis<xmg_network> resyn;
        cut_rewriting( xmg, resyn, cr_ps, &cr_st );
        xmg = cleanup_dangling( xmg );
        
        const auto cec2 = benchmark == "hyp" ? true : abc_cec( xmg, benchmark );
        if (size_before_cr == 0u)
          improv_after_rewrite = 0;
        else
        {
          int diff = size_before_cr - xmg.num_gates();
          improv_after_rewrite = 100 * (double(std::abs(diff))/size_before_cr);
        }

        size_before_resub = xmg.num_gates();

	//xmg_dont_cares_optimization(xmg);
        xmg_resubstitution(xmg, resub_ps, &resub_st);
        xmg = cleanup_dangling( xmg );
        num_gate_profile(xmg,xmg_ps2);
        xmg_ps2.report();
    
        const auto cec = benchmark == "hyp" ? true : abc_cec( xmg, benchmark );

        if (size_before_cr == 0u)
          total_imp = 0;
        else
        {
          int diff = size_before_cr - xmg.num_gates();
          total_imp = 100 * (double(std::abs(diff))/size_before_cr);
        }

        if(xmg_ps.actual_xor3 == 0u )
          rel_xor3 = 0;
        else 
        {
          int diff = xmg_ps.actual_xor3 - xmg_ps2.actual_xor3;
          rel_xor3 = 100 * (double(std::abs(diff))/xmg_ps.actual_xor3); 
          std::cout << "rel_xor " << rel_xor3 <<  std::endl;
        }
        if (xmg_ps.actual_maj == 0u )
          rel_maj = 0;
        else
        {
          int diff = xmg_ps.actual_maj - xmg_ps2.actual_maj; 
          rel_maj  = 100 * (double(std::abs(diff)) / xmg_ps.actual_maj) ; 
          std::cout << "rel_maj " << rel_maj <<  std::endl;
        }

        std::cout <<  "For benchmark "<< benchmark << " improvement after rewrite " << improv_after_rewrite << " and improvement after resub " << total_imp << " at iteration # " << num_iters << std::endl;

        exp3( benchmark, size_before_cr, size_before_resub, xmg.num_gates(), to_seconds( resub_st.time_total ), to_seconds( cr_st.time_total), xmg_ps2.total_xor3, xmg_ps2.actual_xor3, xmg_ps2.actual_xor2, xmg_ps2.total_maj, xmg_ps2.actual_maj, xmg_ps2.remaining_maj, num_iters, improv_after_rewrite, total_imp, cec2, cec );
        //
      exp( benchmark, num_iters, improv_after_rewrite, total_imp, to_seconds( cr_st.time_total ), to_seconds( resub_st.time_total ), xmg_ps.actual_xor3, xmg_ps2.actual_xor3, xmg_ps.actual_maj, xmg_ps2.actual_maj, rel_xor3, rel_maj, cec );
      float area_after = abc_map( xmg, genlib_path );
      exp2 ( benchmark, ( ( ( area_before - area_after ) / area_before ) * 100 ) ) ;

    //} while ( (size_before_cr - xmg.num_gates()) > 0 );
    } while ( total_imp > 0.5 );

    //float area_after = abc_map( xmg, genlib_path );
    //exp2 ( benchmark, ( ( ( area_before - area_after ) / area_before ) * 100 ) ) ;
  }
  
  exp.save();
  exp.table();
  exp2.save();
  exp2.table();
  exp3.save();
  exp3.table();
  return 0;
}
