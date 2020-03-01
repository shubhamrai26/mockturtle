#pragma once

#include <kitty/dynamic_truth_table.hpp>

#include <fmt/format.h>
#include <fstream>
#include <iomanip>
#include <map>
#include <mockturtle/networks/xmg.hpp>
#include <string>
#include <unordered_map>
#include <utility>

namespace percy
{

template<typename ResultType>
class npn_db_populate
{
public:
  explicit npn_db_populate( chain const& c );

  ResultType operator()() const;
};

template<>
class npn_db_populate<std::string>
{
public:
  using signal = mockturtle::xmg_network::signal;

public:
  explicit npn_db_populate( chain const& c, mockturtle::xmg_network& xmg, bool has_constants = false )
      : c( c ), xmg( xmg ), has_constants( has_constants )
  {
  }

  std::string operator()() const
  {
    assert( xmg.num_pis() == c.get_nr_inputs() );
    assert( c.get_nr_outputs() == 1u );

    auto const output_literal = c.get_outputs()[0u];
    auto const output_variable = output_literal >> 1u;
    if ( output_variable == 0u )
    {
      /* special case of constant 0 */
        std::cout << "Constant Condition " << std::endl;
        index_to_signal.emplace_back( xmg.get_constant( false ) );
        return "";
    }

    for ( auto i = 0u; i < c.get_nr_inputs(); ++i )
      index_to_signal.emplace_back( xmg.make_signal( xmg.pi_at( i ) ) );

    std::cout << "Number of steps " << c.get_nr_steps() << std::endl;
    for ( auto i = c.get_nr_inputs(); i < c.get_nr_inputs() + c.get_nr_steps(); ++i )
    {
      auto const nr_in = c.get_nr_inputs();

      auto const& tt = c.get_operator( i - nr_in );

      build_up_xmg( i, tt, true );
    }

    std::cout << "size of index_to_signal " << index_to_signal.size() << std::endl;
    return "";
  }

  void build_up_xmg( uint32_t index, kitty::dynamic_truth_table const& tt, uint32_t const& compl_flag ) const
  {

    kitty::dynamic_truth_table in1( 3 ), in2( 3 ), in3( 3 );
    kitty::dynamic_truth_table cnst0( 3 );

    kitty::create_nth_var( in1, 0 );
    kitty::create_nth_var( in2, 1 );
    kitty::create_nth_var( in3, 2 );
    kitty::create_from_hex_string( cnst0, "00" );

    auto const fanins = c.get_step( index - c.get_nr_inputs() );
    assert( fanins.size() == 3u );

    std::array<signal, 3u> fanin_signals;

    for ( auto i = 0u; i < 3u; ++i )
    {
      assert( fanins.at( i ) < index_to_signal.size() );
      fanin_signals[i] = index_to_signal[fanins.at( i )];
    }

    auto const a_i = fanin_signals[0u];
    auto const b_i = fanin_signals[1u];
    auto const c_i = fanin_signals[2u];

    /* For Majority */
    const auto tt1 = ( in1 & in2 ) | ( in2 & in3 ) | ( in1 & in3 );
    const auto tt2 = ( ~in1 & in2 ) | ( in2 & in3 ) | ( ~in1 & in3 );
    const auto tt3 = ( in1 & ~in2 ) | ( ~in2 & in3 ) | ( in1 & in3 );
    const auto tt4 = ( in1 & in2 ) | ( in2 & ~in3 ) | ( in1 & ~in3 );

    const auto tt5 = ( in1 & in2 ) | ( in2 & cnst0 ) | ( in1 & cnst0 );
    const auto tt6 = ( ~in1 & in2 ) | ( in2 & cnst0 ) | ( ~in1 & cnst0 );
    const auto tt7 = ( in1 & ~in2 ) | ( ~in2 & cnst0 ) | ( in1 & cnst0 );
    const auto tt8 = ( ~in1 & ~in2 ) | ( ~in2 & cnst0 ) | ( ~in1 & cnst0 );

    /* For XOR3 */
    const auto tt9 = in1 ^ in2 ^ in3;
    const auto tt10 = ~in1 ^ in2 ^ in3;

    const auto tt11 = ( ~in1 & ~cnst0 ) | ( ~cnst0 & cnst0 ) | ( ~in1 & cnst0 );
    const auto tt12 = ( in1 & in3 ) | ( in3 & ~cnst0 ) | ( in1 & ~cnst0 );
    const auto tt13 = ( ~in1 & in3 ) | ( in3 & ~cnst0 ) | ( ~in1 & ~cnst0 );
    const auto tt14 = ( in1 & ~in3 ) | ( ~in3 & ~cnst0 ) | ( in1 & ~cnst0 );
    const auto tt15 = ( ~in1 & ~in3 ) | ( ~in3 & ~cnst0 ) | ( ~in1 & ~cnst0 );

    const auto tt16 = ( in2 & in3 ) | ( in3 & ~cnst0 ) | ( in2 & ~cnst0 );
    const auto tt17 = ( ~in2 & in3 ) | ( in3 & ~cnst0 ) | ( ~in2 & ~cnst0 );
    const auto tt18 = ( in2 & ~in3 ) | ( ~in3 & ~cnst0 ) | ( in2 & ~cnst0 );
    const auto tt19 = ( ~in2 & ~in3 ) | ( ~in3 & ~cnst0 ) | ( ~in2 & ~cnst0 );

    const auto tt20 = in1 ^ in2 ^ cnst0;
    const auto tt21 = in2 ^ in3 ^ cnst0;
    const auto tt22 = in1 ^ in3 ^ cnst0;

    signal f;

    if ( tt == tt1 || tt == ~tt1 )
      f = (tt == tt1 ) ? xmg.create_maj( a_i, b_i, c_i ) : ! xmg.create_maj( a_i, b_i, c_i ) ;
    else if ( tt == tt2 || tt == ~tt2 )
      f = (tt == tt2) ? xmg.create_maj( !a_i, b_i, c_i ) : ! xmg.create_maj( !a_i, b_i, c_i ) ;
    else if ( tt == tt3 || tt == ~tt3 )
      f = (tt == tt3) ? xmg.create_maj( a_i, !b_i, c_i ) : ! xmg.create_maj( a_i, !b_i, c_i ) ;
    else if ( tt == tt4 || tt == ~tt4 )
      f = (tt == tt4) ? xmg.create_maj( a_i, b_i, !c_i ) : ! xmg.create_maj( a_i, b_i, !c_i );
    else if ( tt == tt5 || tt == ~tt5 )
      f = (tt == tt5) ? xmg.create_maj( a_i, b_i, xmg.get_constant( false )) : ! xmg.create_maj( a_i, b_i, xmg.get_constant( false )  );
    else if ( tt == tt6 || tt == ~tt6 )
      f = (tt == tt6 ) ? xmg.create_maj( !a_i, b_i, xmg.get_constant( false ) ) : ! xmg.create_maj( !a_i, b_i, xmg.get_constant( false ) );
    else if ( tt == tt7 || tt == ~tt7 )
      f = (tt == tt7) ? xmg.create_maj( a_i, !b_i, xmg.get_constant( false ) ) : !  xmg.create_maj( a_i, !b_i, xmg.get_constant( false ) );
    else if ( tt == tt8 || tt == ~tt8 )
      f = (tt == tt8) ? xmg.create_maj( !a_i, !b_i, xmg.get_constant( false ) ) : ! xmg.create_maj( !a_i, !b_i, xmg.get_constant( false ) );
    else if ( tt == tt9 || tt == ~tt9 )
      f = (tt == tt9) ? xmg.create_xor3( a_i, b_i, c_i ) : ! xmg.create_xor3( a_i, b_i, c_i );
    else if ( tt == tt10 || tt == ~tt10 )
      f = (tt == tt10) ? xmg.create_xor3( !a_i, b_i, c_i ) : ! xmg.create_xor3( !a_i, b_i, c_i );
    else if ( tt == tt11 || tt == ~tt11 )
      f = (tt == tt11) ? xmg.create_maj( !a_i, xmg.get_constant( true ), xmg.get_constant( false ) ) : ! xmg.create_maj( !a_i, xmg.get_constant( true ), xmg.get_constant( false ) );
    else if ( tt == tt12 || tt == ~tt12 )
    {
        std::cout << "Check  " << ( (tt == tt12) ? "normal " : "compl" ) << std::endl;
      f = (tt == tt12) ? xmg.create_maj( a_i, xmg.get_constant( true ), c_i ) : ! xmg.create_maj( a_i, xmg.get_constant( true ), c_i );
    }
    else if ( tt == tt13 || tt == ~tt13 )
      f = (tt == tt13) ? xmg.create_maj( !a_i, xmg.get_constant( true ), c_i ) :  ! xmg.create_maj( !a_i, xmg.get_constant( true ), c_i ) ;
    else if ( tt == tt14 || tt == ~tt14 )
      f = (tt == tt14) ? xmg.create_maj( a_i, xmg.get_constant( true ), !c_i ) : ! xmg.create_maj( a_i, xmg.get_constant( true ), !c_i );
    else if ( tt == tt15 || tt == ~tt15 )
      f = (tt == tt15) ? xmg.create_maj( !a_i, xmg.get_constant( true ), !c_i ) : ! xmg.create_maj( !a_i, xmg.get_constant( true ), !c_i );
    else if ( tt == tt16 || tt == ~tt16 )
      f = (tt == tt16) ? xmg.create_maj( xmg.get_constant( true ), b_i, c_i ) : ! xmg.create_maj( xmg.get_constant( true ), b_i, c_i );
    else if ( tt == tt17 || tt == ~tt17 )
      f = (tt == tt17) ? xmg.create_maj( xmg.get_constant( true ), !b_i, c_i ) : ! xmg.create_maj( xmg.get_constant( true ), !b_i, c_i );
    else if ( tt == tt18 || tt == ~tt18 )
      f =  (tt == tt18) ?  xmg.create_maj( xmg.get_constant( true ), b_i, !c_i ) : !xmg.create_maj( xmg.get_constant( true ), b_i, !c_i );
    else if ( tt == tt19 || tt == ~tt19 )
      f = (tt == tt19) ?  xmg.create_maj( xmg.get_constant( true ), !b_i, !c_i ) : !  xmg.create_maj( xmg.get_constant( true ), !b_i, !c_i );
    else if ( tt == tt20 || tt == ~tt20 )
      f = (tt == tt20) ? xmg.create_xor3( a_i, b_i, xmg.get_constant( false ) ) : ! xmg.create_xor3( a_i, b_i, xmg.get_constant( false ) );
    else if ( tt == tt21 || tt == ~tt21 )
      f = (tt == tt21) ? xmg.create_xor3( xmg.get_constant( false ), b_i, c_i ) : ! xmg.create_xor3( xmg.get_constant( false ), b_i, c_i );
    else if ( tt == tt22 || tt == ~tt22 )
      f = (tt == tt22) ?  xmg.create_xor3( a_i, xmg.get_constant( false ), c_i ) : ! xmg.create_xor3( a_i, xmg.get_constant( false ), c_i );
    else
    {
      kitty::print_binary( tt );
      std::cout << " No match is found " << std::endl;
    }

    uint32_t next_index = index_to_signal.size();
    assert( next_index == index );
    index_to_signal.emplace_back( f );
  }

protected:
  chain const& c;
  mockturtle::xmg_network& xmg;
  bool has_constants = false;
  mutable std::vector<signal> index_to_signal;

  /* maps the function to their string representation, e.g.,  2 --> "(%s%s)" */
  std::map<uint32_t, std::string> masks;
}; /* npn_db_populate */

} // namespace percy
