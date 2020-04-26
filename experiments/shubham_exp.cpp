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

#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/collapse_mapped.hpp>
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/algorithms/detail/database_generator.hpp>
#include <mockturtle/algorithms/lut_mapping.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
#include <mockturtle/algorithms/node_resynthesis/cached.hpp>
#include <mockturtle/algorithms/node_resynthesis/xmg4_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/xmg_npn.hpp>
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

template<typename Ntk>
mockturtle::klut_network lut_map( Ntk const& ntk, uint32_t k = 4 )
{
  mockturtle::write_verilog( ntk, "/tmp/network.v" );
  system( fmt::format( "abc -q \"/tmp/network.v; &get; &if -a -K {}; &put; write_blif /tmp/output.blif\"", k ).c_str() );

  mockturtle::klut_network klut;
  if ( lorina::read_blif( "/tmp/output.blif", mockturtle::blif_reader( klut ) ) != lorina::return_code::success )
  {
    std::cout << "ERROR 1" << std::endl;
    std::abort();
    return klut;
  }
  return klut;
}

void example1()
{
  /* enumerate NPN representatives */
  std::unordered_set<kitty::dynamic_truth_table, kitty::hash<kitty::dynamic_truth_table>> classes;
  kitty::dynamic_truth_table tt( 4u );
  do
  {
    const auto res = kitty::exact_npn_canonization( tt );
    classes.insert( std::get<0>( res ) );
    kitty::next_inplace( tt );
  } while ( !kitty::is_const0( tt ) );
  std::cout << "[i] enumerated "
            << ( 1 << ( 1 << tt.num_vars() ) ) << " functions into "
            << classes.size() << " classes." << std::endl;

  /* generate database with exact XMG synthesis */
  mockturtle::xmg_network xmg;
  mockturtle::exact_xmg_resynthesis<mockturtle::xmg_network> exact( {.use_only_self_dual_gates = true} );
  mockturtle::detail::database_generator dbgen( xmg, exact, {} );
  for ( const auto& f : classes )
  {
    dbgen.add_function( f );

    std::cout << ".";
    std::cout.flush();
  }
  mockturtle::write_verilog( xmg, "db.v" );
}

void example2()
{
  using namespace experiments;

  /* load database from file */
  mockturtle::xmg_network xmg2_db, xmg3_db, xmgs_db, techlib_db;
  read_verilog( "xmg2_db.v", mockturtle::verilog_reader( xmg2_db ) );
  read_verilog( "xmg3_db.v", mockturtle::verilog_reader( xmg3_db ) );
  read_verilog( "xmgs_db.v", mockturtle::verilog_reader( xmgs_db ) );
  read_verilog( "techlib.v", mockturtle::verilog_reader( techlib_db ) );

  /* option 1: XMG strategy using databse from file */
  mockturtle::xmg4_npn_resynthesis<mockturtle::xmg_network> xmg2_resyn( mockturtle::detail::to_index_list( xmg2_db ) );
  mockturtle::xmg4_npn_resynthesis<mockturtle::xmg_network> xmg3_resyn( mockturtle::detail::to_index_list( xmg3_db ) );
  mockturtle::xmg4_npn_resynthesis<mockturtle::xmg_network> xmgs_resyn( mockturtle::detail::to_index_list( xmgs_db ) );
  mockturtle::xmg4_npn_resynthesis<mockturtle::xmg_network> techlib_resyn( mockturtle::detail::to_index_list( techlib_db ) );

  experiments::experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, bool>
    exp( "shubham", "benchmark", "LUTs", "XMG2", "XMG3", "XMGs", "TechLib", "CEC" );

  for ( auto const& benchmark : epfl_benchmarks( experiments::arithmetic ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    /* read aig */
    mockturtle::aig_network aig;
    if ( lorina::read_aiger( experiments::benchmark_path( benchmark ), mockturtle::aiger_reader( aig ) ) != lorina::return_code::success )
    {
      std::cout << "ERROR 2" << std::endl;
      std::abort();
      return;
    }

    /* LUT map AIG into k-LUT network */
    auto klut = lut_map( aig, 4u );

    /* resynthesize klut with resynthesis strategies */
    mockturtle::xmg_network xmg2, xmg3, xmgs, techlib;
    xmg2 = mockturtle::node_resynthesis<mockturtle::xmg_network>( klut, xmg2_resyn );
    xmg3 = mockturtle::node_resynthesis<mockturtle::xmg_network>( klut, xmg3_resyn );
    xmgs = mockturtle::node_resynthesis<mockturtle::xmg_network>( klut, xmgs_resyn );
    techlib = mockturtle::node_resynthesis<mockturtle::xmg_network>( klut, techlib_resyn );

    bool cec = true;
    if ( benchmark != "hyp" )
    {
      cec &= abc_cec( xmg2, benchmark );
      cec &= abc_cec( xmg3, benchmark );
      cec &= abc_cec( xmgs, benchmark );
      cec &= abc_cec( techlib, benchmark );
    }

    exp( benchmark, klut.size(), xmg2.size(), xmg3.size(), xmgs.size(), techlib.size(), cec );

    exp.save();
    exp.table();
  }

  exp.save();
  exp.table();
}

void example3()
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
  if ( read_verilog( "xmg_npn4_db.v", mockturtle::verilog_reader( db ) ) != lorina::return_code::success )
  {
    std::cout << "ERROR" << std::endl;
    std::abort();
    return;
  }
  else
  {
    std::cout << "[i] DB loaded" << std::endl;
  }

  /* generate resynthesis strategy */

  /* option 1: X3MG strategy using databse from file */
  mockturtle::xmg4_npn_resynthesis<mockturtle::xmg_network> npn_resyn( mockturtle::detail::to_index_list( db ) );

  /* option 2: X2MG strategy */
  // mockturtle::xmg_npn_resynthesis npn_resyn;

  experiments::experiment<std::string, uint32_t, uint32_t, std::string, uint32_t, uint32_t, std::string>
    exp( "cut_rewriting", "benchmark", "size aspdac", "size ours", "xmg improv", "klut6 aspdac", "klut6 ours", "klut improv" );
  for ( auto const& benchmark : epfl_benchmarks( experiments::arithmetic ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    /* read aig */
    mockturtle::aig_network aig;
    if ( lorina::read_aiger( experiments::benchmark_path( benchmark ), mockturtle::aiger_reader( aig ) ) != lorina::return_code::success )
    {
      std::cout << "ERROR 2" << std::endl;
      std::abort();
      return;
    }

    /* LUT map AIG into k-LUT network */
    auto klut = lut_map( aig );

    mockturtle::xmg_network xmg;
    while ( true )
    {
      auto const klut_size_before = klut.size();
      xmg = mockturtle::node_resynthesis<mockturtle::xmg_network>( klut, npn_resyn );

      /* option 3: ABC if-mapping flow */
      auto const new_klut = lut_map( xmg );

      /* option 4: mockturtle mf-mapping flow */
      // mockturtle::lut_mapping_params ps;
      // ps.cut_enumeration_ps.cut_size = 4;
      // ps.cut_enumeration_ps.cut_limit = 16;
      //
      // mockturtle::mapping_view<mockturtle::xmg_network, true> mapped_xmg{xmg};
      // mockturtle::lut_mapping<decltype( mapped_xmg ), true>( mapped_xmg, ps );
      // const auto new_klut = *mockturtle::collapse_mapped_network<mockturtle::klut_network>( mapped_xmg );

      if ( new_klut.size() >= klut_size_before )
        break;

      klut = new_klut;
    }

    auto const klut6 = lut_map( xmg, 6u );

    std::cout << "final XMG size = " << xmg.size() << std::endl;
    std::cout << "final KLUT-6 size = " << klut6.size() << std::endl;

    exp( benchmark,
         std::get<0>( aspdac17_xmg[benchmark] ), xmg.size(),
         fmt::format( "{:3.2f}", (double(std::get<0>( aspdac17_xmg[benchmark] ) ) - xmg.size())/std::get<0>( aspdac17_xmg[benchmark] ) ),
         std::get<1>( aspdac17_xmg[benchmark] ), klut6.size(),
         fmt::format( "{:3.2f}", (double(std::get<1>( aspdac17_xmg[benchmark] ) ) - klut6.size())/std::get<1>( aspdac17_xmg[benchmark] ) ) );
  }

  exp.save();
  exp.table();
}

namespace mockturtle
{

struct exact_techmap_params
{
  uint32_t conflict_limit = 1000;
};

template<class Ntk = xmg_network>
class exact_techmap_resynthesis
{
public:
  using signal = mockturtle::signal<Ntk>;

public:
  explicit exact_techmap_resynthesis( exact_techmap_params const& ps = {} )
    : ps( ps )
  {
  }

  template<typename LeavesIterator, typename TT, typename Fn>
  bool operator()( Ntk& ntk, TT const& function, LeavesIterator begin, LeavesIterator end, Fn&& fn ) const
  {
    std::cout << function.num_vars() << " "; kitty::print_hex( function ); std::cout << std::endl;

    static_assert( kitty::is_complete_truth_table<TT>::value, "Truth table must be complete" );
    auto const tt = function.num_vars() < 3u ? kitty::extend_to( function, 3u ) : function;
    bool const normal = kitty::is_normal( tt );

    percy::chain chain;
    percy::spec spec;
    spec.conflict_limit = int32_t( ps.conflict_limit );
    spec.verbosity = 0;
    spec.fanin = 3;

    /* specify local normalized gate primitives */
    kitty::dynamic_truth_table const0{3};
    kitty::dynamic_truth_table a{3};
    kitty::dynamic_truth_table b{3};
    kitty::dynamic_truth_table c{3};
    kitty::create_nth_var( a, 0 );
    kitty::create_nth_var( b, 1 );
    kitty::create_nth_var( c, 2 );

    spec.add_primitive( const0 ); // 00
    spec.add_primitive( a ); // aa
    spec.add_primitive( b ); // cc
    spec.add_primitive( c ); // f0

    spec.add_primitive( a & b ); // 88
    spec.add_primitive( ~a & b ); // 44
    spec.add_primitive( a & ~b ); // 22
    spec.add_primitive( a & c ); // a0
    spec.add_primitive( ~a & c ); // 50
    spec.add_primitive( a & ~c ); // 0a
    spec.add_primitive( b & c ); // c0
    spec.add_primitive( ~b & c ); // 30
    spec.add_primitive( b & ~c ); // 0c

    spec.add_primitive( kitty::ternary_majority( a, b, c ) ); // e8
    spec.add_primitive( kitty::ternary_majority( ~a, b, c ) ); // d4
    spec.add_primitive( kitty::ternary_majority( a, ~b, c ) ); // b2
    spec.add_primitive( kitty::ternary_majority( a, b, ~c ) ); // 8e

    spec.add_primitive( a ^ b ); // 66
    spec.add_primitive( a ^ c ); // 5a
    spec.add_primitive( b ^ c ); // 3c
    spec.add_primitive( a ^ b ^ c ); // 96

    percy::bsat_wrapper solver;
    percy::ssv_encoder encoder(solver);

    spec[0] = normal ? tt : ~tt;

    for ( auto i = 0u; i < 1u; ++i )
    {
      auto const result = percy::next_struct_solution( spec, chain, solver, encoder );
      if ( result != percy::success )
        break;

      assert( result == percy::success );

      auto const sim = chain.simulate();
      assert( chain.simulate()[0] == spec[0] );

      std::vector<signal> signals( tt.num_vars(), ntk.get_constant( false ) );
      std::copy( begin, end, signals.begin() );

      for ( auto i = 0; i < chain.get_nr_steps(); ++i )
      {
        auto const c1 = signals[chain.get_step( i )[0]];
        auto const c2 = signals[chain.get_step( i )[1]];
        auto const c3 = signals[chain.get_step( i )[2]];

        switch( chain.get_operator( i )._bits[0] )
        {
        case 0x00:
          signals.emplace_back( ntk.get_constant( false ) );
          break;
        case 0xe8:
          signals.emplace_back( ntk.create_maj( c1,  c2,  c3 ) );
          break;
        case 0xd4:
          signals.emplace_back( ntk.create_maj( !c1,  c2,  c3 ) );
          break;
        case 0xb2:
          signals.emplace_back( ntk.create_maj( c1,  !c2,  c3 ) );
          break;
        case 0x8e:
          signals.emplace_back( ntk.create_maj( c1,  c2,  !c3 ) );
          break;
        case 0x66:
          signals.emplace_back( ntk.create_xor( c1,  c2 ) );
          break;
        case 0x5a:
          signals.emplace_back( ntk.create_xor( c1,  c3 ) );
          break;
        case 0x3c:
          signals.emplace_back( ntk.create_xor( c2,  c3 ) );
          break;
        case 0x96:
          signals.emplace_back( ntk.create_xor3( c1,  c2,  c3 ) );
          break;
        case 0x88:
          signals.emplace_back( ntk.create_and( c1, c2 ) );
          break;
        case 0x44:
          signals.emplace_back( ntk.create_and( !c1, c2 ) );
          break;
        case 0x22:
          signals.emplace_back( ntk.create_and( c1, !c2 ) );
          break;
        case 0xa0:
          signals.emplace_back( ntk.create_and( c1, c3 ) );
          break;
        case 0x50:
          signals.emplace_back( ntk.create_and( !c1, c3 ) );
          break;
        case 0x0a:
          signals.emplace_back( ntk.create_and( c1, !c3 ) );
          break;
        case 0xc0:
          signals.emplace_back( ntk.create_and( c2, c3 ) );
          break;
        case 0x30:
          signals.emplace_back( ntk.create_and( !c2, c3 ) );
          break;
        case 0x0c:
          signals.emplace_back( ntk.create_and( c2, !c3 ) );
          break;
        default:
          std::cerr << "[e] unsupported operation " << kitty::to_hex( chain.get_operator( i ) ) << "\n";
          assert( false );
          break;
        }
      }

      assert( chain.get_outputs().size() > 0u );
      uint32_t const output_index = ( chain.get_outputs()[0u] >> 1u );
      auto const output_signal = output_index == 0u ? ntk.get_constant( false ) : signals[output_index - 1];
      if ( !fn( chain.is_output_inverted( 0 ) ^ normal ? output_signal : !output_signal ) )
      {
        return false; /* quit */
      }
    }
    return false;
  }

private:
  exact_techmap_params const& ps;
}; /* exact_xmg_resynthesis */

} /* mockturtle */

void example4()
{
  /* enumerate NPN representatives */
  std::unordered_set<kitty::dynamic_truth_table, kitty::hash<kitty::dynamic_truth_table>> classes;
  kitty::dynamic_truth_table tt( 4u );
  do
  {
    const auto res = kitty::exact_npn_canonization( tt );
    classes.insert( std::get<0>( res ) );
    kitty::next_inplace( tt );
  } while ( !kitty::is_const0( tt ) );
  std::cout << "[i] enumerated "
            << ( 1 << ( 1 << tt.num_vars() ) ) << " functions into "
            << classes.size() << " classes." << std::endl;

  /* generate database with exact XMG synthesis */
  mockturtle::xmg_network xmg;
  mockturtle::exact_techmap_resynthesis<mockturtle::xmg_network> exact;
  mockturtle::detail::database_generator dbgen( xmg, exact, {} );
  for ( const auto& f : classes )
  {
    dbgen.add_function( f );

    std::cout << ".";
    std::cout.flush();
  }
  mockturtle::write_verilog( xmg, "db.v" );
}

void example5()
{
  using namespace experiments;

  uint32_t const size = 6u;
  
  /* prepare XMG DB */
  mockturtle::exact_xmg_resynthesis<mockturtle::xmg_network> xmg2_exact( {.use_xor3 = false} );
  mockturtle::cached_resynthesis<mockturtle::xmg_network, decltype( xmg2_exact )> cached_xmg2_exact( xmg2_exact, size, "exact_xmg2_cache6.v" ); 
  
  mockturtle::exact_xmg_resynthesis<mockturtle::xmg_network> xmg3_exact( {.use_xor3 = true} );
  mockturtle::cached_resynthesis<mockturtle::xmg_network, decltype( xmg3_exact )> cached_xmg3_exact( xmg3_exact, size, "exact_xmg3_cache6.v" ); 
  
  /* prepare exact techlib resyn */
  mockturtle::exact_techmap_resynthesis<mockturtle::xmg_network> techlib_exact;
  mockturtle::cached_resynthesis<mockturtle::xmg_network, decltype( techlib_exact )> cached_techlib_exact( techlib_exact, size, "exact_techlib_cache6.v" );

  experiments::experiment<std::string, uint32_t, uint32_t, bool>
    exp( "shubham", "benchmark", "LUTs", "TechLib", "CEC" );

  for ( auto const& benchmark : epfl_benchmarks( experiments::arithmetic & ~experiments::hyp ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    /* read aig */
    mockturtle::aig_network aig;
    if ( lorina::read_aiger( experiments::benchmark_path( benchmark ), mockturtle::aiger_reader( aig ) ) != lorina::return_code::success )
    {
      std::cout << "ERROR 2" << std::endl;
      std::abort();
      return;
    }

    /* LUT map AIG into k-LUT network */
    auto klut = lut_map( aig, size );

    /* translate to XMG */
    mockturtle::xmg_network xmg;
    xmg = mockturtle::node_resynthesis<mockturtle::xmg_network>( klut, cached_xmg3_exact );

    /* rewrite KLUT into tech-mapped KLUT */
    mockturtle::cut_rewriting_params ps;
    ps.cut_enumeration_ps.cut_size = size;

    auto new_xmg = cut_rewriting( xmg, cached_techlib_exact, ps );
    auto cec = abc_cec( xmg, benchmark );

    exp( benchmark, klut.size(), new_xmg.size(), cec );

    exp.save();
    exp.table();
  }

  exp.save();
  exp.table();
}

int main()
{
  example5();
  return 0;
}
