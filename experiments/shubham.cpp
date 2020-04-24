#include <lorina/lorina.hpp>
#include <kitty/kitty.hpp>
#include <percy/percy.hpp>

#include <mockturtle/networks/xmg.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/io/index_list.hpp>

template<class Ntk = mockturtle::xmg_network>
class exact_xmg_resynthesis
{
public:
  explicit exact_xmg_resynthesis() = default;

  template<typename LeavesIterator, typename Fn>
  void operator()( Ntk& ntk, kitty::dynamic_truth_table const& function, LeavesIterator begin, LeavesIterator end, Fn&& fn ) const
  {
    using node = mockturtle::node<Ntk>;
    using signal = mockturtle::signal<Ntk>;
    
    percy::chain chain;
    percy::spec spec;
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

    spec.add_primitive( const0 );
    spec.add_primitive( a );
    spec.add_primitive( b );
    spec.add_primitive( c );
    spec.add_primitive( kitty::ternary_majority(  a,  b,  c ) );
    spec.add_primitive( kitty::ternary_majority( ~a,  b,  c ) );
    spec.add_primitive( kitty::ternary_majority(  a, ~b,  c ) );
    spec.add_primitive( kitty::ternary_majority(  a,  b, ~c ) );
    spec.add_primitive( a ^ b ^ c );

    percy::bsat_wrapper solver;
    percy::ssv_encoder encoder(solver);
    
    spec[0] = kitty::is_normal( function ) ? function : ~function;

    for ( auto i = 0u; i < 10u; ++i )
    {
      // auto const result = percy::synthesize( spec, chain );      
      auto const result = percy::next_struct_solution( spec, chain, solver, encoder );
      if ( result != percy::success )
        break;
      
      assert( result == percy::success );
      
      auto const sim = chain.simulate();
      assert( chain.simulate()[0] == spec[0] );
      
      std::vector<signal> signals( begin, end );
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
      
        case 0x96:
          signals.emplace_back( ntk.create_xor3( c1,  c2,  c3 ) );
          continue;
          
        default:
          std::cerr << "[e] unsupported operation " << kitty::to_hex( chain.get_operator( i ) ) << "\n";
          assert( false );
          break;
        }
      }

      fn( chain.is_output_inverted( 0 ) ? !signals.back() : signals.back() );
    }
  }
}; /* exact_xmg */

template<typename Ntk, typename ResynFn>
class exact_database_generator
{
public:
  using node = mockturtle::node<Ntk>;
  using signal = mockturtle::signal<Ntk>;
    
  explicit exact_database_generator( Ntk& ntk, ResynFn const& resyn, uint32_t num_vars )
    : ntk( ntk )
    , resyn( resyn )
    , num_vars( num_vars )
  {
    for ( auto i = 0u; i < num_vars; ++i )
    {
      pis.emplace_back( ntk.create_pi() );
    }
  }

  void add_function( kitty::dynamic_truth_table tt )
  {
    /* normalize the function first if necessary */
    if ( !kitty::is_normal( tt ) )
    {
      tt = ~tt;
    }

    /* resynthesize the function and add it to the database */
    resyn( ntk, tt, std::begin( pis ), std::end( pis ),
           [&]( const signal& s )
           {
             std::cout << "[i] function: ";
             kitty::print_binary( tt );
             std::cout << " stored at PO #" << ntk.num_pos() << std::endl;
             ntk.create_po( s );
           } );
  }

  Ntk& ntk;
  ResynFn const& resyn;
  
  uint32_t num_vars;
  std::vector<signal> pis;
}; /* exact_database_generator */

int main()
{
  /* compute NPN classe */
  std::unordered_set<kitty::dynamic_truth_table, kitty::hash<kitty::dynamic_truth_table>> classes;
  kitty::dynamic_truth_table tt( 4u );
  do
  {
    /* apply NPN canonization and add resulting representative to set */
    const auto res = kitty::exact_npn_canonization( tt );
    classes.insert( std::get<0>( res ) );

    /* increment truth table */
    kitty::next_inplace( tt );
  } while ( !kitty::is_const0( tt ) );

  std::cout << "[i] enumerated "
            << ( 1 << ( 1 << tt.num_vars() ) ) << " functions into "
            << classes.size() << " classes." << std::endl;

  /* synthesize database */
  mockturtle::xmg_network xmg;
  exact_xmg_resynthesis<mockturtle::xmg_network> resyn;
  exact_database_generator<mockturtle::xmg_network, decltype( resyn )> gen( xmg, resyn, 4u );
  for ( const auto& cls : classes )
  {
    gen.add_function( cls );
  }

  std::cout << xmg.size() << std::endl;

  /* write into a file */
  mockturtle::write_verilog( xmg, "shared_xmg.v" );

  std::cout << mockturtle::detail::to_index_list_string( xmg ) << std::endl;
  
  return 0;
}
