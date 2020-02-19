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
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/verilog.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/resubstitution.hpp>
#include <mockturtle/algorithms/xmg_resub.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/xmg.hpp>
#include <mockturtle/properties/xmgcost.hpp>

#include <experiments.hpp>

int main()
{
    using namespace experiments;
    using namespace mockturtle;

  experiment<std::string, uint32_t, float, float,  uint32_t, uint32_t, uint32_t, uint32_t, float, float, bool> exp( "xmg_resub", 
            "benchmark", "iter.", "rel_imp", "runtime", "xor3", "xor3'", "maj", "maj'", "xor3_imp", "maj_imp", "equivalent" );


  for ( auto const& benchmark : epfl_benchmarks() )
  {
    if (benchmark != "adder") 
        continue;
    fmt::print( "[i] processing {}\n", benchmark );
    xmg_network xmg;
    lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xmg ) );

    int area_before = abc_map(xmg);
    // Preparing for the xmg_cost calculation 
    xmg_cost_params xmg_ps,xmg_ps2;

    // For XMG Resubstitution
    resubstitution_params ps;
    resubstitution_stats st;
    ps.max_pis = 8u;
    ps.max_inserts = 1u;  // Discuss with Heinz once.
    ps.progress = true;

    uint32_t size_before; 
    uint32_t num_iters = 0u;
    float improvements = 0;
    float rel_xor3 = 0;
    float rel_maj = 0;

    do 
    {
	    xmg_ps.reset();
      xmg_ps2.reset();
      num_iters++;
      size_before = xmg.num_gates();
      
      num_gate_profile(xmg,xmg_ps);
      xmg_resubstitution(xmg, ps, &st);

      xmg = cleanup_dangling( xmg );

      num_gate_profile(xmg,xmg_ps2);

      const auto cec = benchmark == "hyp" ? true : abc_cec( xmg, benchmark );
      std::cout << "size_before " <<  size_before << std::endl;
      std::cout << "xmg num_gates " <<  xmg.num_gates() << std::endl;
      if (size_before == 0u)
        improvements = 0;
      else 
      {
        int diff = size_before - xmg.num_gates();
        improvements = 100 * (double(std::abs(diff))/size_before);
        std::cout << " improvements " << improvements <<  std::endl;
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

      std::cout <<  "For benchmark "<< benchmark << "improvement " << improvements << "at iteration # " << num_iters << std::endl;
        
      exp( benchmark, num_iters, improvements, to_seconds( st.time_total ), xmg_ps.actual_xor3, xmg_ps2.actual_xor3,xmg_ps.actual_maj, xmg_ps2.actual_maj, rel_xor3, rel_maj, cec );

      std::cout << "Trying out Mapping" << std::endl;

    } while ((size_before - xmg.num_gates()) > 0);
    
    int area_after = abc_map(xmg);
    std:: cout << "improvement in area after mapping "  << (area_after - area_before) << std::endl;

    // Figure out how to integrate the xmg_cost.hpp as well  
    
  }
  
  exp.save();
  exp.table();
  return 0;
}
