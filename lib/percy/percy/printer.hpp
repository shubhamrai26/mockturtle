#pragma once

#include <kitty/dynamic_truth_table.hpp>

#include <fmt/format.h>

#include <iomanip>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>

namespace percy
{

template<typename ResultType>
class printer
{
public:
  explicit printer( chain const& c );

  ResultType operator()() const;
};

template<>
class printer<std::string>
{
public:
  explicit printer( chain const& c, bool has_constants = false )
      : c( c ), has_constants( has_constants )
  {
  }

  std::string operator()() const
  {
    assert( c.get_nr_outputs() == 1u );

    std::string result;
    auto const output_literal = c.get_outputs()[0u];

    if ( output_literal & 1 )
      result += "!";

    auto const output_variable = output_literal >> 1u;
    if ( output_variable == 0u )
    {
      /* special case of constant 0 */
      result += "0";
    }
    else
    {
      result += step_to_expression( output_variable - 1u );
    }

    return result;
  }

  std::string step_to_expression( uint32_t index ) const
  {
    auto const nr_in = c.get_nr_inputs();
    if ( has_constants && index == 0 )
    {
      return "0";
    }
    else if ( index < nr_in )
    {

      return std::string( 1u, 'a' + index );
    }

    const auto& step = c.get_step( index - nr_in );
    auto const word = word_from_tt( c.get_operator( index - nr_in ) );
    auto const it = masks.find( word );
    if ( it == masks.end() )
    {
      assert( false );
      return "error";
    }

    auto const mask = it->second;

    std::string result;
    switch ( c.get_fanin() )
    {
    case 2u:
      return fmt::format( mask,
                          step_to_expression( step.at( 0 ) ),
                          step_to_expression( step.at( 1 ) ) );
    case 3u:
      return fmt::format( mask,
                          step_to_expression( step.at( 0 ) ),
                          step_to_expression( step.at( 1 ) ),
                          step_to_expression( step.at( 2 ) ) );
    case 4u:
      return fmt::format( mask,
                          step_to_expression( step.at( 0 ) ),
                          step_to_expression( step.at( 1 ) ),
                          step_to_expression( step.at( 2 ) ),
                          step_to_expression( step.at( 3 ) ) );
    case 5u:
      return fmt::format( mask,
                          step_to_expression( step.at( 0 ) ),
                          step_to_expression( step.at( 1 ) ),
                          step_to_expression( step.at( 2 ) ),
                          step_to_expression( step.at( 3 ) ),
                          step_to_expression( step.at( 4 ) ) );
    default:
      return mask;
    }
  }

  void add_function( kitty::dynamic_truth_table const& tt, std::string const& mask )
  {
    masks[word_from_tt( tt )] = mask;
  }

protected:
  static inline uint32_t word_from_tt( kitty::dynamic_truth_table const& tt )
  {
    auto word = 0;
    for ( int i = 0; i < ( 1 << tt.num_vars() ); ++i )
    {
      if ( kitty::get_bit( tt, i ) )
        word |= ( 1 << i );
    }
    return word;
  }

  bool is_xor3( kitty::dynamic_truth_table const& tt ) const
  {
    kitty::dynamic_truth_table xor3_tt( 3 );
    kitty::create_parity( xor3_tt );
    if ( tt != xor3_tt )
    {
      return false;
    }
    return true;
  }

protected:
  chain const& c;

  bool has_constants = false;

  /* maps the function to their string representation, e.g.,  2 --> "(%s%s)" */
  std::map<uint32_t, std::string> masks;
}; /* printer */

} // namespace percy
