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

#include "experiments.hpp"
#include <mockturtle/algorithms/resubstitution.hpp>
#include <mockturtle/algorithms/xmg_resub.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/collapse_mapped.hpp>
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/algorithms/lut_mapping.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/xmg_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/xmg4_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/xmg3_npn.hpp>
#include <mockturtle/properties/xmgcost.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/index_list.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/views/mapping_view.hpp>
#include <lorina/lorina.hpp>
#include <fmt/format.h>

#include <string>
#include <vector>

std::string const genlib_path = "/home/shubham/My_work/abc-vlsi-cad-flow/std_libs/date_lib_count_tt_4.genlib";
template<typename Ntk>
mockturtle::klut_network lut_map( Ntk const& ntk, uint32_t k = 4 )
{
  mockturtle::write_verilog( ntk, "/tmp/ex_network.v" );
  system( fmt::format( "abc -q \"/tmp/ex_network.v; &get; &if -a -K {}; &put; write_blif /tmp/ex_output.blif\"", k ).c_str() );
  
  mockturtle::klut_network klut;
  if ( lorina::read_blif( "/tmp/ex_output.blif", mockturtle::blif_reader( klut ) ) != lorina::return_code::success )
  {
    std::cout << "ERROR 1" << std::endl;
    std::abort();
    return klut;
  }
  return klut;
}

int main()
{
  using namespace experiments;

  /* Winston & Mathias's results from ASP-DAC'17 */
  std::map<std::string, std::tuple<uint32_t, uint32_t>> aspdac17_xmg;
  aspdac17_xmg.emplace( "adder", std::make_tuple( 639, 251 ) );
  aspdac17_xmg.emplace( "bar", std::make_tuple( 3281, 888 ) );
  aspdac17_xmg.emplace( "div", std::make_tuple( 29607, 12094 ) );
  aspdac17_xmg.emplace( "hyp", std::make_tuple( 155349, 50835 ) );
  aspdac17_xmg.emplace( "log2", std::make_tuple( 27936, 8438 ) );
  aspdac17_xmg.emplace( "max", std::make_tuple( 2296, 745 ) );
  aspdac17_xmg.emplace( "multiplier", std::make_tuple( 17508, 5700 ));
  aspdac17_xmg.emplace( "sin", std::make_tuple( 5100, 1655 ) );
  aspdac17_xmg.emplace( "sqrt", std::make_tuple( 20130, 6595 ) );
  aspdac17_xmg.emplace( "square", std::make_tuple( 15070, 3969 ) );
  
  /* load database from file */
  mockturtle::xmg_network db;
  if ( read_verilog( "xmg3_without_sd.v", mockturtle::verilog_reader( db ) ) != lorina::return_code::success )
  {
    std::cout << "ERROR" << std::endl;
    std::abort();
    return -1;
  }
  else
  {
    std::cout << "[i] DB loaded" << std::endl;
  }
  
  /* generate resynthesis strategy */

  /* option 1: X3MG strategy using databse from file */
  mockturtle::xmg4_npn_resynthesis<mockturtle::xmg_network> npn_resyn( mockturtle::detail::to_index_list( db ) );

  /* option 2: X2MG strategy */
  //mockturtle::xmg_npn_resynthesis npn_resyn;
  experiments::experiment<std::string, double, double, double, double, bool >
    exp( "RFET_area", "benchmark", "c2rs_area", "init_area", "final_area", "klut6_size", "equiv" );
  
  //experiments::experiment<std::string, uint32_t, uint32_t, std::string, uint32_t, uint32_t, std::string>
  //  exp( "cut_rewriting", "benchmark", "size_aspdac", "size_ours", "xmg_improv", "klut6_aspdac", "klut6_ours", "klut_improv" );
  //for ( auto const& benchmark : epfl_benchmarks( experiments::arithmetic ) )
  
  for ( auto const& benchmark : epfl_benchmarks( ) )
  {
      if (benchmark != "div" && benchmark != "voter" )
          continue;
    fmt::print( "[i] processing {}\n", benchmark );

    /* read aig */
    mockturtle::aig_network aig;
    if ( lorina::read_aiger( experiments::benchmark_path( benchmark ), mockturtle::aiger_reader( aig ) ) != lorina::return_code::success )
    {
      std::cout << "ERROR 2" << std::endl;
      std::abort();
      return -1;
    }

    /* LUT map AIG into k-LUT network */
    auto klut = lut_map( aig );
    double c2rs_area = abc_map_compress2rs( klut, genlib_path ); 
    double init_area = abc_map( klut, genlib_path);

    mockturtle::xmg_cost_params ps1;
    ps1.reset();

    mockturtle::xmg_network xmg;
    while ( true )
    {
      auto const klut_size_before = klut.size();
      xmg = mockturtle::node_resynthesis<mockturtle::xmg_network>( klut, npn_resyn );

      /* resubstitution */ 
      mockturtle::resubstitution_params resub_ps;
      mockturtle::resubstitution_stats resub_st;
      resub_ps.max_pis = 8u;
      //resub_ps.progress = true;
      resub_ps.max_inserts = 1u;  
      resub_ps.use_dont_cares = true; 
      resub_ps.window_size = 12u;  
      mockturtle::xmg_resubstitution(xmg, resub_ps, &resub_st);
      xmg = mockturtle::cleanup_dangling( xmg );

      /* option 3: ABC if-mapping flow */
      //auto const new_klut = lut_map( xmg );

      /* option 4: mockturtle mf-mapping flow */
       mockturtle::lut_mapping_params ps;
       ps.cut_enumeration_ps.cut_size = 4;
       ps.cut_enumeration_ps.cut_limit = 16;
       
       mockturtle::mapping_view<mockturtle::xmg_network, true> mapped_xmg{xmg};
       mockturtle::lut_mapping<decltype( mapped_xmg ), true>( mapped_xmg, ps );
       const auto new_klut = *mockturtle::collapse_mapped_network<mockturtle::klut_network>( mapped_xmg );
      
      if ( new_klut.size() >= klut_size_before )
        break;

      klut = new_klut;
    }

    num_gate_profile( xmg, ps1 );
    ps1.report( );

    //auto const klut6 = lut_map( xmg, 6u );
    mockturtle::lut_mapping_params ps;
    ps.cut_enumeration_ps.cut_size = 6;
    ps.cut_enumeration_ps.cut_limit = 16;

    mockturtle::mapping_view<mockturtle::xmg_network, true> mapped_xmg{xmg};
    mockturtle::lut_mapping<decltype( mapped_xmg ), true>( mapped_xmg, ps );
    const auto new_klut6 = *mockturtle::collapse_mapped_network<mockturtle::klut_network>( mapped_xmg );
      

      
    std::cout << "final XMG size = " << xmg.size() << std::endl;
    std::cout << "final KLUT-6 size = " << new_klut6.size() << std::endl;

    const auto cec42 = benchmark == "hyp" ? true : abc_cec( xmg, benchmark );
    std::cout <<"Equivalence " << (cec42 ? "true" : "false") << std::endl;
    
    double final_area = abc_map( xmg, genlib_path );

    exp( benchmark, c2rs_area, init_area, final_area, new_klut6.size(), cec42 );

    //exp( benchmark,
    //     std::get<0>( aspdac17_xmg[benchmark] ), xmg.size(),
    //     fmt::format( "{:3.2f}", (double(std::get<0>( aspdac17_xmg[benchmark] ) ) - xmg.size())/std::get<0>( aspdac17_xmg[benchmark] ) ),
    //     std::get<1>( aspdac17_xmg[benchmark] ), klut6.size(),
    //     fmt::format( "{:3.2f}", (double(std::get<1>( aspdac17_xmg[benchmark] ) ) - klut6.size())/std::get<1>( aspdac17_xmg[benchmark] ) ) );
  exp.save();
  exp.table();
  }
  
  exp.save();
  exp.table();

  return 0;
}
