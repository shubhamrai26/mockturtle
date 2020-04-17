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

/*!
  \file xmg4_npn.hpp
  \brief Replace with size-optimum XAGs and AIGs from NPN (from ABC rewrite)

  \author Mathias Soeken
*/

#pragma once

#include "../../algorithms/simulation.hpp"
#include "../../io/index_list.hpp"
#include "../../networks/xmg.hpp"
#include "../../utils/node_map.hpp"
#include "../../utils/stopwatch.hpp"

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/print.hpp>
#include <kitty/static_truth_table.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cstdint>
#include <vector>

namespace mockturtle
{

struct xmg4_npn_resynthesis_params
{
  /*! \brief Be verbose. */
  bool verbose{false};
};

struct xmg4_npn_resynthesis_stats
{
  stopwatch<>::duration time_classes{0};
  stopwatch<>::duration time_db{0};

  uint32_t db_size;
  uint32_t covered_classes;

  void report() const
  {
    std::cout << fmt::format( "[i] build classes time = {:>5.2f} secs\n", to_seconds( time_classes ) );
    std::cout << fmt::format( "[i] build db time      = {:>5.2f} secs\n", to_seconds( time_db ) );
  }
};

/*! \brief Resynthesis function based on pre-computed AIGs.
 *
 * This resynthesis function can be passed to ``cut_rewriting``.  It will
 * produce a network based on pre-computed XAGs with up to at most 4 variables.
 * Consequently, the nodes' fan-in sizes in the input network must not exceed
 * 4.
 *
   \verbatim embed:rst

   Example

   .. code-block:: c++

      const aig_network aig = ...;
      xmg4_npn_resynthesis<aig_network> resyn;
      aig = cut_rewriting( aig, resyn );

   .. note::

      The implementation of this algorithm was heavily inspired by the rewrite
      command in AIG.  It uses the same underlying database of subcircuits.
   \endverbatim
 */
template<class Ntk, class DatabaseNtk = xmg_network>
class xmg4_npn_resynthesis
{
public:
  xmg4_npn_resynthesis( std::vector<uint32_t> const& subgraphs, xmg4_npn_resynthesis_params const& ps = {}, xmg4_npn_resynthesis_stats* pst = nullptr )
      : ps( ps )
      , pst( pst )
      , subgraphs( subgraphs )
      , _repr( 1u << 16u )
  {
    static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
    static_assert( has_get_constant_v<Ntk>, "Ntk does not implement the get_constant method" );
    static_assert( has_create_and_v<Ntk>, "Ntk does not implement the create_and method" );
    static_assert( has_create_xor_v<Ntk>, "Ntk does not implement the create_xor method" );
    static_assert( has_create_not_v<Ntk>, "Ntk does not implement the create_not method" );

    static_assert( is_network_type_v<DatabaseNtk>, "DatabaseNtk is not a network type" );
    static_assert( has_get_node_v<DatabaseNtk>, "DatabaseNtk does not implement the get_node method" );
    static_assert( has_is_complemented_v<DatabaseNtk>, "DatabaseNtk does not implement the is_complemented method" );
    static_assert( has_is_xor_v<DatabaseNtk>, "DatabaseNtk does not implement the is_xor method" );
    static_assert( has_size_v<DatabaseNtk>, "DatabaseNtk does not implement the size method" );
    static_assert( has_create_pi_v<DatabaseNtk>, "DatabaseNtk does not implement the create_pi method" );
    static_assert( has_create_and_v<DatabaseNtk>, "DatabaseNtk does not implement the create_and method" );
    static_assert( has_create_xor_v<DatabaseNtk>, "DatabaseNtk does not implement the create_xor method" );
    static_assert( has_foreach_fanin_v<DatabaseNtk>, "DatabaseNtk does not implement the foreach_fanin method" );
    static_assert( has_foreach_node_v<DatabaseNtk>, "DatabaseNtk does not implement the foreach_node method" );
    static_assert( has_make_signal_v<DatabaseNtk>, "DatabaseNtk does not implement the make_signal method" );

    build_classes();
    build_db();
  }

  virtual ~xmg4_npn_resynthesis()
  {
    if ( ps.verbose )
    {
      st.report();
    }

    if ( pst )
    {
      *pst = st;
    }
  }

  template<typename LeavesIterator, typename Fn>
  void operator()( Ntk& ntk, kitty::dynamic_truth_table const& function, LeavesIterator begin, LeavesIterator end, Fn&& fn ) const
  {
    kitty::static_truth_table<4u> tt = kitty::extend_to<4u>( function );

    /* get representative of function */
    const auto [repr, phase, perm] = _repr[*tt.cbegin()];

    /* check if representative has circuits */
    const auto it = _repr_to_signal.find( repr );
    if ( it == _repr_to_signal.end() )
    {
      return;
    }

    std::vector<signal<Ntk>> pis( 4, ntk.get_constant( false ) );
    std::copy( begin, end, pis.begin() );

    std::unordered_map<node<DatabaseNtk>, signal<Ntk>> db_to_ntk;
    db_to_ntk.insert( {0, ntk.get_constant( false )} );
    for ( auto i = 0; i < 4; ++i )
    {
      db_to_ntk.insert( {i + 1, ( phase >> perm[i] & 1 ) ? ntk.create_not( pis[perm[i]] ) : pis[perm[i]]} );
    }

    for ( auto const& cand : it->second )
    {
      const auto f = copy_db_entry( ntk, _db.get_node( cand ), db_to_ntk );
      if ( !fn( _db.is_complemented( cand ) != ( phase >> 4 & 1 ) ? ntk.create_not( f ) : f ) )
      {
        return;
      }
    }
  }

private:
  signal<Ntk>
  copy_db_entry( Ntk& ntk, node<DatabaseNtk> const& n, std::unordered_map<node<DatabaseNtk>, signal<Ntk>>& db_to_ntk ) const
  {
    if ( const auto it = db_to_ntk.find( n ); it != db_to_ntk.end() )
    {
      return it->second;
    }

    std::array<signal<Ntk>, 3> fanin{};
    _db.foreach_fanin( n, [&]( auto const& f, auto i ) {
      const auto ntk_f = copy_db_entry( ntk, _db.get_node( f ), db_to_ntk );
      fanin[i] = _db.is_complemented( f ) ? ntk.create_not( ntk_f ) : ntk_f;
    } );

    const auto f = _db.is_xor3( n ) ? ntk.create_xor3( fanin[0], fanin[1], fanin[2] ) : ntk.create_maj( fanin[0], fanin[1], fanin[2] );
    db_to_ntk.insert( {n, f} );
    return f;
  }

  void build_classes()
  {
    stopwatch t( st.time_classes );

    kitty::static_truth_table<4u> tt;
    do {
      _repr[*tt.cbegin()] = kitty::exact_npn_canonization( tt );
      kitty::next_inplace( tt );
    } while ( !kitty::is_const0( tt ) );
  }

  void build_db()
  {
    stopwatch t( st.time_db );

    _db = create_from_ternary_index_list<DatabaseNtk>( subgraphs.begin() );
    const auto sim_res = simulate_nodes<kitty::static_truth_table<4u>>( _db );

    _db.foreach_node( [&]( auto n ) {
      if ( std::get<0>( _repr[*sim_res[n].cbegin()] ) == sim_res[n] )
      {
        if ( _repr_to_signal.count( sim_res[n] ) == 0 )
        {
          _repr_to_signal.insert( {sim_res[n], {_db.make_signal( n )}} );
        }
        else
        {
          _repr_to_signal[sim_res[n]].push_back( _db.make_signal( n ) );
        }
      }
      else
      {
        const auto f = ~sim_res[n];
        if ( std::get<0>( _repr[*f.cbegin()] ) == f )
        {
          if ( _repr_to_signal.count( f ) == 0 )
          {
            _repr_to_signal.insert( {f, {!_db.make_signal( n )}} );
          }
          else
          {
            _repr_to_signal[f].push_back( !_db.make_signal( n ) );
          }
        }
      }
    } );

    st.db_size = _db.size();
    st.covered_classes = static_cast<uint32_t>( _repr_to_signal.size() );
  }

  xmg4_npn_resynthesis_params ps;
  xmg4_npn_resynthesis_stats st;
  xmg4_npn_resynthesis_stats* pst{nullptr};

  std::vector<std::tuple<kitty::static_truth_table<4u>, uint32_t, std::vector<uint8_t>>> _repr;
  std::unordered_map<kitty::static_truth_table<4u>, std::vector<signal<DatabaseNtk>>, kitty::hash<kitty::static_truth_table<4u>>> _repr_to_signal;

  DatabaseNtk _db;

  std::vector<uint32_t> const subgraphs;  
}; // namespace mockturtle

} /* namespace mockturtle */
