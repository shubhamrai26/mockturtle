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
  \file xmg3_npn.hpp
  \brief Replace with size-optimum xmg3s and AIGs from NPN (from ABC rewrite)

  \author Mathias Soeken
*/

#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

#include <fmt/format.h>
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/print.hpp>
#include <kitty/static_truth_table.hpp>

#include "../../algorithms/simulation.hpp"
#include "../../io/write_bench.hpp"
#include "../../networks/xmg.hpp"
#include "../../utils/node_map.hpp"
#include "../../utils/stopwatch.hpp"

namespace mockturtle
{

struct xmg3_npn_resynthesis_params
{
  /*! \brief Be verbose. */
  bool verbose{false};
};

struct xmg3_npn_resynthesis_stats
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
 * produce a network based on pre-computed xmg3s with up to at most 4 variables.
 * Consequently, the nodes' fan-in sizes in the input network must not exceed
 * 4.
 *
   \verbatim embed:rst

   Example

   .. code-block:: c++

      const aig_network aig = ...;
      xmg3_npn_resynthesis<aig_network> resyn;
      cut_rewriting( aig, resyn );

   .. note::

      The implementation of this algorithm was heavily inspired by the rewrite
      command in AIG.  It uses the same underlying database of subcircuits.
   \endverbatim
 */
template<class Ntk, class DatabaseNtk = xmg_network>
class xmg3_npn_resynthesis
{
public:
  xmg3_npn_resynthesis( xmg3_npn_resynthesis_params const& ps = {}, xmg3_npn_resynthesis_stats* pst = nullptr )
      : ps( ps ),
        pst( pst ),
        _classes( 1 << 16 )
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

    _repr.reserve( 222u );
    build_classes();
    build_db();
  }

  virtual ~xmg3_npn_resynthesis()
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
  void operator()( Ntk& ntk, kitty::dynamic_truth_table const& function, LeavesIterator begin, LeavesIterator end, Fn&& fn )
  {
    kitty::static_truth_table<4> tt = kitty::extend_to<4>( function );

    /* get representative of function */
    const auto repr = _repr[_classes[*tt.cbegin()]];

    /* check if representative has circuits */
    const auto it = _repr_to_signal.find( repr );
    if ( it == _repr_to_signal.end() )
    {
      return;
    }

    const auto config = kitty::exact_npn_canonization( tt );

    assert( repr == std::get<0>( config ) );

    std::vector<signal<Ntk>> pis( 4, ntk.get_constant( false ) );
    std::copy( begin, end, pis.begin() );

    std::vector<signal<Ntk>> pis_perm;
    auto perm = std::get<2>( config );
    for ( auto i = 0; i < 4; ++i )
    {
      pis_perm.push_back( pis[perm[i]] );
    }

    const auto& phase = std::get<1>( config );
    for ( auto i = 0; i < 4; ++i )
    {
      if ( ( phase >> perm[i] ) & 1 )
      {
        pis_perm[i] = ntk.create_not( pis_perm[i] );
      }
    }

    for ( auto const& cand : it->second )
    {
      std::unordered_map<node<DatabaseNtk>, signal<Ntk>> db_to_ntk;

      db_to_ntk.insert( {0, ntk.get_constant( false )} );
      for ( auto i = 0; i < 4; ++i )
      {
        db_to_ntk.insert( {i + 1, pis_perm[i]} );
      }
      auto f = copy_db_entry( ntk, _db.get_node( cand ), db_to_ntk );
      if ( _db.is_complemented( cand ) != ( ( phase >> 4 ) & 1 ) )
      {
        f = ntk.create_not( f );
      }
      if ( !fn( f ) )
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

    std::vector<signal<Ntk>> fanin;
    //std::array<signal<Ntk>, 2> fanin;
    _db.foreach_fanin( n, [&]( auto const& f ) {
      auto ntk_f = copy_db_entry( ntk, _db.get_node( f ), db_to_ntk );
      if ( _db.is_complemented( f ) )
      {
        ntk_f = ntk.create_not( ntk_f );
      }
      fanin.push_back( ntk_f );
    } );

    const auto f = _db.is_xor3( n ) ? ntk.create_xor3( fanin[0], fanin[1], fanin[2] ) : ntk.create_maj( fanin[0], fanin[1], fanin[2] );
    db_to_ntk.insert( {n, f} );
    return f;
  }

  void build_classes()
  {
    stopwatch t( st.time_classes );

    kitty::dynamic_truth_table map( 16u );
    std::transform( map.cbegin(), map.cend(), map.begin(), []( auto word ) { return ~word; } );

    int64_t index = 0;
    kitty::static_truth_table<4> tt;
    while ( index != -1 )
    {
      kitty::create_from_words( tt, &index, &index + 1 );
      const auto res = kitty::exact_npn_canonization( tt, [&]( const auto& tt ) {
        _classes[*tt.cbegin()] = _repr.size();
        kitty::clear_bit( map, *tt.cbegin() );
      } );
      _repr.push_back( std::get<0>( res ) );

      /* find next non-classified truth table */
      index = find_first_one_bit( map );
    }
  }

  void build_db()
  {
    stopwatch t( st.time_db );

    /* four primary inputs */
    _db.create_pi();
    _db.create_pi();
    _db.create_pi();
    _db.create_pi();

    auto* p = subgraphs;
    while ( true )
    {
      auto entry0 = *p++;
      auto entry1 = *p++;
      auto entry2 = *p++;

      if ( entry0 == 0 && entry1 == 0 && entry2 == 0)
        break;

      auto is_xor = entry0 & 1;
      entry0 >>= 1;

      const auto child0 = _db.make_signal( entry0 >> 1 ) ^ ( entry0 & 1 );
      const auto child1 = _db.make_signal( entry1 >> 1 ) ^ ( entry1 & 1 );
      const auto child2 = _db.make_signal( entry2 >> 1 ) ^ ( entry2 & 1 );

      if ( is_xor )
      {
        _db.create_xor3( child0, child1, child2 );
      }
      else
      {
        _db.create_maj( child0, child1, child2 );
      }
    }

    const auto sim_res = simulate_nodes<kitty::static_truth_table<4>>( _db );

    _db.foreach_node( [&]( auto n ) {
      if ( _repr[_classes[*sim_res[n].cbegin()]] == sim_res[n] )
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
        if ( _repr[_classes[*f.cbegin()]] == f )
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
    st.covered_classes = _repr_to_signal.size();
  }

  xmg3_npn_resynthesis_params ps;
  xmg3_npn_resynthesis_stats st;
  xmg3_npn_resynthesis_stats* pst{nullptr};

  std::vector<kitty::static_truth_table<4>> _repr;
  std::vector<uint8_t> _classes;
  std::unordered_map<kitty::static_truth_table<4>, std::vector<signal<DatabaseNtk>>, kitty::hash<kitty::static_truth_table<4>>> _repr_to_signal;

  DatabaseNtk _db;

  // clang-format off
  inline static const uint16_t subgraphs[]
  {
      0x2,0x2,0x4,0xc,0x8,0xa,0x0,0xa,0xd,0x8,0x7,0x8,0x2,0x6,0x10,0x2,0x6,0x8,0x29,0x4,0x2,0x0,0x15,0x16,0x0,0x3,0x4,0x8,0x14,0x1b,0x8,0x7,0x1a,0x3d,0x6,0x2,0x0,0x9,0x20,0x4,0x4,0x6,0x2,0x8,0x24,0x49,0x6,0x0,0x0,0x9,0x28,0x8,0x7,0xa,0x0,0x2,0x2c,0x10,0xb,0x2e,0x6,0x4,0x8,0x4,0x6,0x8,0x2,0x32,0x34,0x15,0x6,0x0,0x2,0x8,0x38,0x4,0x4,0x9,0x0,0x7,0x3c,0x10,0x3c,0x3e,0x11,0x6,0x4,0x0,0x9,0x42,0x0,0x6,0x9,0x8d,0x4,0x0,0x12,0xa,0x48,0x0,0x4,0x7,0x8,0x8,0x4d,0xd,0x4,0x2,0x2,0x8,0x50,0x0,0x2,0x15,0x10,0x16,0x55,0x11,0x4,0x2,0x0,0x6,0x59,0x0,0x9,0x5a,0x4,0x4,0x7,0xbd,0x8,0x6,0x0,0x9,0x60,0x4,0x4,0x8,0x0,0x6,0x65,0xc,0x8,0x66,0x4,0x6,0x51,0xd5,0x50,0x8,0x0,0x6b,0x6c,0x2,0x4,0x6,0x8,0x8,0x70,0x4,0x15,0x72,0x6,0x4,0x6,0x4,0x8,0x76,0xed,0x8,0x4,0x0,0x79,0x7a,0x69,0x4,0x2,0x2,0x34,0x7e,0x11,0x6,0x2,0x8,0x6,0x8,0x0,0x82,0x85,0x0,0x2,0x35,0x8,0x34,0x88,0x0,0x5,0x82,0xc,0x8,0x8c,0xa,0x8c,0x8e,0xc,0x8,0xb,0x85,0x6,0x2,0xc,0x8,0x95,0x0,0x42,0x97,0x0,0x2,0x5,0x135,0x6,0x0,0x2,0x84,0x9c,0x0,0x2,0x9,0xc,0x64,0xa1,0x4,0x7,0x42,0x4,0x44,0xa5,0x12,0xa,0x42,0x85,0x4,0x2,0x0,0x7,0xaa,0x10,0x43,0xac,0x0,0x2,0x8,0x2,0x84,0xb0,0x0,0x3,0x6,0x169,0x8,0x2,0x0,0x51,0xb6,0x4,0x6,0x25,0x0,0x9,0x24,0x8,0xba,0xbd,0x109,0x8,0x2,0x0,0x8,0x84,0x0,0xc0,0xc2,0x10,0x50,0x71,0x0,0x7,0x42,0x10,0x59,0xc8,0xc,0x8,0x59,0x0,0x3,0x58,0x10,0xcc,0xcf,0x4,0x5,0x8,0x0,0x4,0xd3,0x6,0x6,0xd4,0x1ad,0x8,0x6,0xc,0x8,0x64,0x8,0x7,0x32,0x1b9,0x4,0x0,0x0,0x4,0x8,0x0,0x2,0xe0,0x8,0x7,0xe2,0x1c9,0x4,0x0,0x0,0x8,0xd2,0x4,0x6,0xe9,0x1d5,0xd2,0x0,0x0,0x2,0x59,0x8,0x6,0xee,0x1e1,0xee,0x0,0x0,0x2,0x4,0xc,0x58,0xf5,0x1ed,0x6,0x0,0x0,0x8,0x64,0x0,0x6,0xfb,0x1f9,0x64,0x6,0x2,0x2,0x6,0x8,0xa1,0x100,0x205,0x4,0x0,0x65,0x4,0x0,0x2,0x8,0x32,0xc,0x107,0x108,0xc9,0x8,0x0,0xc,0x8,0x10d,0x11,0x4,0x0,0x221,0x6,0x2,0x4,0x8,0x113,0x10,0x111,0x114,0x0,0x112,0x117,0xb1,0x6,0x0,0x12,0x100,0x11a,0x8,0x6,0x9b,0x23d,0x32,0x0,0x2,0x8,0x100,0x9,0x2,0x0,0x200,0x122,0x125,0x0,0x5,0x6,0x251,0x8,0x2,0x8,0x8,0x128,0x0,0x12a,0x12d,0x10,0x71,0x124,0x1a5,0x6,0x4,0x1a5,0x2,0x0,0x12,0x132,0x134,0x161,0x4,0x2,0xc,0x8,0x139,0x0,0xb0,0x13b,0x4,0x8,0x71,0x27d,0x8,0x0,0xe1,0x2,0x0,0x0,0x3,0x70,0x12,0x142,0x144,0x0,0x6,0x8,0x291,0x8,0x2,0x8,0x8,0x149,0x2,0x14a,0x14c,0x12,0x70,0x142,0x6,0x6,0x8,0x0,0x8,0x153,0x8,0x6,0x154,0x2ad,0x152,0x0,0x0,0x7,0xd2,0x2b5,0x152,0x0,0x11,0x2,0x0,0x10,0xc8,0x15f,0x0,0x6,0x1b,0x2c5,0x1a,0x2,0xc,0x8,0x165,0xc,0x9,0x32,0x2d1,0x4c,0x32,0x4,0x4,0x4d,0x99,0x6,0x2,0x10,0x16c,0x16e,0xc,0x9,0x128,0x4,0x4,0x173,0x2e9,0x128,0x8,0xa,0x6,0x15e,0xc,0x8,0x179,0x0,0x4,0x9,0x8,0x6,0x15e,0x2fd,0x8,0x0,0x2bc,0x17d,0x180,0xc9,0x50,0x0,0x10,0x50,0x184,0x0,0x3,0x84,0x311,0x8,0x2,0x4,0x7,0xa,0x12,0x38,0x18c,0x15,0x8,0x0,0x2,0x6,0x190,0x14,0x190,0x193,0xc,0x8,0x190,0x2,0x6,0xa,0x331,0x8,0x6,0x8,0x6,0x9,0x0,0x8,0x19d,0x33d,0x6,0x4,0x8,0x6,0x43,0x0,0x2,0x43,0x2,0x1a2,0x1a4,0x2,0x2,0x8,0x351,0x128,0x4,0x0,0x4,0x1a8,0xc,0x1a8,0x1ad,0x2,0x110,0x1a8,0x2,0x6,0x1b0,0x365,0x1a8,0x0,0x8,0x8,0x83,0x36d,0x9a,0x0,0x0,0x3,0x8,0x375,0x6,0x2,0x2,0x4,0x1bc,0x10,0x1bc,0x1bf,0x289,0x6,0x2,0x12,0x70,0x1c2,0x8,0x8,0x1a9,0x2,0x6,0x1c6,0x391,0x1a8,0x0,0x0,0x2,0x6,0x10,0x71,0x1cc,0x4,0x6,0x10,0x0,0x9,0x1d0,0x3a5,0x10,0x4,0x141,0x6,0x4,0x0,0xe1,0x1d6,0xa,0x6,0x8,0x3b5,0xa0,0x6,0xa,0x8,0xa0,0x8,0x6,0x1de,0x3c1,0x1de,0xa0,0x399,0x8,0x4,0xc,0x46,0x1e5,0x0,0x7,0xa0,0x12,0x1d6,0x1e8,0x4,0x4,0x85,0x0,0x8,0x1ed,0x3dd,0x1ec,0x84,0x351,0x4,0x0,0xc,0x8,0x1f3,0x10,0x4c,0x65,0x3ed,0x64,0x6,0x84,0xe1,0x1a8,0x8,0x7,0xa0,0x3f9,0x8,0x4,0x21,0x4,0x0,0x4,0x6,0x11,0x20,0x201,0x202,0xd,0x4,0x0,0x0,0xe0,0x207,0x4,0x6,0x9,0x0,0x9,0x20a,0x419,0x84,0x0,0x4,0x4,0x111,0xc,0x8,0x210,0x425,0x210,0x110,0x1c0,0x101,0x206,0xc,0x8,0x111,0x8,0x6,0xa0,0x435,0x8,0x0,0xe,0xa0,0x1da,0x43d,0x8,0x0,0x4,0x7,0x128,0x2,0x8,0x222,0x449,0x128,0x4,0x10,0x24,0x206,0xc,0x64,0x9b,0x455,0x8,0x0,0x6,0x8,0x10,0x0,0xa,0x11,0x461,0x22e,0x4,0x0,0x2,0x207,0x469,0x4,0x2,0x10,0x206,0x237,0x251,0x8,0x4,0x0,0x3,0x23a,0x250,0x23b,0x23c,0x291,0x8,0x4,0x14,0x148,0x241,0x2,0x8,0xd2,0x489,0x70,0x0,0x99,0x8,0x6,0x4,0x4,0x82,0x495,0x2,0x0,0x4,0x4,0x83,0x0,0x5,0x8,0x104,0x24f,0x250,0x4a5,0x4,0x0,0x11,0x6,0x0,0x4ad,0x4,0x2,0xc,0x8,0x258,0x134,0x256,0x25b,0xa,0x6,0x46,0x4bd,0x152,0x0,0x0,0x7,0x8,0x10,0x125,0x262,0x4,0x8,0x46,0xc,0x125,0x266,0x2bd,0x6,0x4,0xc,0x46,0x26a,0x8,0x6,0x15f,0x2,0x6,0x26e,0x10,0x26a,0x270,0xc,0x8,0x125,0x2,0x4,0x8,0x248,0x256,0x277,0x2a5,0x4,0x0,0x8,0x6,0x153,0x10,0x27a,0x27c,0x48,0x46,0x125,0x0,0x8,0x82,0x6,0x4,0x282,0x509,0x82,0x0,0x4,0x7,0x8,0xd,0x2,0x0,0x6,0x288,0x28a,0x8,0x289,0x28c,0x51d,0x28a,0x8,0x0,0x9,0x84,0x4,0x6,0x292,0x529,0x84,0x0,0x8,0x6,0x110,0x10,0x112,0x298,0x351,0x6,0x0,0x10,0x125,0x29c,0x0,0x59,0x256,0x48,0x46,0x58,0x4,0x6,0x263,0x549,0x262,0x4,0x10,0x262,0x2a7,0x8d,0x4,0x2,0xc,0x46,0x2aa,0x8,0x8,0x1cd,0x55d,0x4,0x2,0xc,0x2ae,0x2b0,0xc,0x8,0x58,0xc,0x8,0x124,0x0,0x276,0x2b6,0x571,0x6,0x0,0x6,0xe0,0x206,0x579,0x8,0x4,0xa0,0x256,0x277,0x0,0x4,0x6,0x10,0x152,0x2c2,0x589,0x2c2,0x4,0x10,0x76,0x124,0x591,0x6,0x0,0x8,0x7,0xe0,0x1c1,0x2,0x0,0x0,0x2cc,0x2ce,0x5a1,0x8,0x6,0xa0,0x9a,0x148,0x8,0x6,0x125,0x249,0x6,0x0,0x10,0x2d6,0x2d9,0x4,0x6,0xe1,0x5b9,0x4,0x2,0x0,0x2dc,0x2de,0x5c1,0x2de,0x8,0x351,0x6,0x2,0x4,0x4,0x2e4,0x5cd,0x8,0x2,0x0,0x6,0x125,0x248,0x256,0x2eb,0x2,0x2,0x256,0x5dd,0x4,0x2,0x10,0x257,0x2ee,0x4ae,0x2f0,0x2f2,0x12,0x256,0x278,0x12,0x124,0x2d8,0x0,0x6,0x124,0x5f5,0x124,0x8,0x6a,0x124,0x256,0x4,0x4,0x149,0x0,0x148,0x300,0x605,0x8,0x6,0x68,0x7e,0x257,0x0,0x2,0x1db,0x611,0x8,0x6,0x4,0x4,0x257,0x12,0x256,0x30c,0x61d,0x256,0x0,0x61a,0x30e,0x310,0x515,0x8,0x2,0x0,0x5,0x314,0x514,0x315,0x316,0x66,0x100,0x256,0x6,0x6,0x32,0x0,0x8,0x31d,0x63d,0x32,0x6,0x4,0x6,0x33,0x645,0xb0,0x6,0x611,0x276,0x6,0xa,0x6,0x1a,0x0,0x8,0x1b,0x655,0x328,0x4,0xc,0x8,0x3c,0x65d,0x8,0x0,0xe,0x8,0x3c,0x0,0x6,0x333,0x669,0x3c,0x0,0x0,0x4,0x3c,0xc,0x9,0x338,0x675,0x3c,0x0,0x79,0x6,0x0,0x10,0x19c,0x33e,0xc,0xb,0x58,0x685,0x8,0x0,0x4,0x6,0x9b,0x8,0x9,0x9a,0x691,0x346,0x0,0x8,0x6,0x1b,0xc,0x256,0x34d,0x8,0x8,0x257,0x4,0x257,0x350,0x140,0xf5,0x256,0x0,0x6,0xf5,0x6ad,0xf4,0x8,0x4,0x5,0x6,0x8,0x43,0x35a,0x6b9,0x2,0x0,0x0,0x2,0x7,0x6,0x8,0x360,0x6,0x4,0x362,0x6c9,0x42,0x0,0x0,0x8,0x19c,0x6d1,0x4,0x2,0x0,0x368,0x36b,0x6d9,0x19c,0x6,0x85,0x2,0x0,0xc,0x8,0x370,0x84,0x370,0x372,0xc,0x9,0xf4,0x6ed,0x28a,0x4,0x2,0x8,0xa,0xc,0xf4,0x37b,0x6f9,0xa,0x0,0x65,0x6,0x4,0x201,0x8,0x4,0x0,0x381,0x382,0x6,0x8,0x2c2,0x70d,0x6,0x4,0x1a5,0x6,0x2,0x1a5,0x4,0x0,0x0,0x38b,0x38c,0x0,0x3,0x32,0x64,0x380,0x390,0x0,0x6,0x3d,0x729,0x2,0x0,0x72d,0x8,0x4,0xc,0x8,0x70,0x4,0x4,0x39a,0x739,0x70,0x8,0x2,0x4,0x32,0x65,0x6,0x0,0x4,0x3a1,0x3a2,0x749,0x6,0x4,0x0,0x8,0x3c,0x0,0x6,0x3a8,0x755,0x4,0x2,0x8,0xd2,0x38a,0x639,0x262,0x4,0x4,0x9,0x262,0xa,0x8,0x3b2,0x769,0x262,0x2,0x4,0x8,0x263,0xa,0x6,0x3b8,0x775,0x262,0x2,0x0,0x2,0x58,0xc,0x8,0x3bf,0x781,0x3be,0x58,0x515,0x4,0x0,0x12,0x28a,0x3c4,0xe,0x3c4,0x3c6,0x0,0x4,0x1bb,0x0,0x6,0x3ca,0x799,0x58,0x0,0xc,0x72,0x15f,0x7a1,0x70,0x0,0x2,0x6,0xf4,0xa,0x8,0xf4,0x7ad,0x3d4,0x2,0x585,0x8,0x2,0x0,0x2,0x3db,0x7b9,0x6,0x4,0xa,0x6,0xb4,0x7c1,0x58,0x0,0x0,0x4,0x207,0x6,0x8,0x206,0x7cd,0x3e4,0x2,0xe,0x110,0x15e,0x4,0x7,0x15e,0x7d9,0x6,0x4,0x6,0x3ec,0x3ee,0x8,0x8,0x15,0x4,0x6,0x3f2,0x7e5,0x14,0x0,0x7ed,0x3f4,0x4,0xa,0x6,0x1ba,0x7f5,0xe0,0x2,0x1c1,0x6,0x4,0x4,0x4,0x3ff,0x801,0x8,0x0,0x49,0x8,0x0,0x49,0x8,0x6,0xc,0x14,0xf5,0x811,0xf4,0x2,0x815,0x8,0x4,0x8,0x6,0xb5,0x4,0x8,0xb4,0x821,0x40e,0x6,0x0,0x4,0x1cd,0x829,0x8,0x2,0x4,0x8,0x42,0x831,0x6,0x0,0x0,0x2,0x206,0x839,0x8,0x4,0x15,0x8,0x6,0x515,0x8,0x4,0x0,0x0,0x0
  };
  // clang-format on
}; // namespace mockturtle

} /* namespace mockturtle */
