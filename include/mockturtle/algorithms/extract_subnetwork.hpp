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

#pragma once

#include <vector>
#include <iostream>

namespace mockturtle
{

/*! \brief Fanout-free cone cut
 *
 *  Eagerly compute a fanout-free cut into fanin-direction.
 */
template<typename Ntk>
class ffc_cut
{
public:
  using node = typename Ntk::node;

  struct subnetwork
  {
    std::vector<node> leaves{};
    std::vector<node> roots{}; /* not used */
    std::vector<node> divs{}; /* not used */
  }; /* subnetwork */

public:
  explicit ffc_cut( Ntk const& ntk, int32_t cut_size = -1 )
    : ntk( ntk )
    , cut_size( cut_size )
  {
  }

  std::vector<subnetwork> operator()( node const& root )
  {
    std::vector<node> leaves = { root };
    expand_fanin_cut( leaves );
    std::sort( std::begin( leaves ), std::end( leaves ) );
    return { subnetwork{leaves} };
  }

private:
  void expand_fanin_cut( std::vector<node>& leaves )
  {
    while ( true )
    {
      /* step 1: select a node from the leaves to expand the cut */
      auto it = std::begin( leaves );
      while ( it != std::end( leaves ) )
      {
        /* skip PIs */
        if ( ntk.is_pi( *it ) )
        {
          ++it;
          continue;
        }

        /* skip node if not fanout-free */
        if ( ntk.fanout_size( *it ) > 1 )
        {
          ++it;
          continue;
        }

        if ( cut_size > 0 && leaves.size() - 1 + ntk.fanin_size( *it ) > uint32_t( cut_size ) )
        {
          ++it;
          continue;
        }

        /* found a possible node at which we should expand */
        break;
      }

      /* if we cannot find a leave to expand, we are done with this cut */
      if ( it == std::end( leaves ) )
        return;

      /* step 2: expand the cut, i.e., remove the node from the cut and add its fanin */
      auto const node = *it;
      leaves.erase( it );

      ntk.foreach_fanin( node, [&]( auto const &f ){
          auto const n = ntk.get_node( f );

          /* unique push back */
          if ( std::find( std::begin( leaves ), std::end( leaves ), n ) == std::end( leaves ) )
            leaves.push_back( n );
        });

      /* sort the leaves */
      std::sort( std::begin( leaves ), std::end( leaves ) );
    }
  }

private:
  Ntk const& ntk;
  int32_t cut_size;
}; /* ffc_cut */

/*! \brief Extended FFC cut.
 *
 * Eagerly compute a fanout-free cone cut into fanin-direction and
 * collect additional divisors that depend on the same support.
 */
template<typename Ntk>
class xcut
{
public:
  using node = typename Ntk::node;

  struct subnetwork
  {
    std::vector<node> leaves{};
    std::vector<node> roots{}; /* not used */
    std::vector<node> divs{};
  }; /* subnetwork */

public:
  explicit xcut( Ntk const& ntk, int32_t cut_size = -1 )
    : ntk( ntk )
    , cut_size( cut_size )
  {
  }

  std::vector<subnetwork> operator()( node const& root )
  {
    /* register two traversal ids */
    ntk.incr_trav_id();
    cover_id = ntk.trav_id();

    ntk.incr_trav_id();
    divisor_id = ntk.trav_id();

    std::vector<node> leaves = { root };
    expand_fanin_cut( leaves );
    std::sort( std::begin( leaves ), std::end( leaves ) );

    /* skip all the computations of the divisors if the leave size is too small */
    if ( leaves.size() <= 2u )
    {
      return { subnetwork{leaves,{}} };
    }

    /* mark leaves visited */
    for ( const auto& l : leaves )
    {
      ntk.set_visited( l, cover_id );
    }

    /* TODO: could be replaced with a cheaper depth check */
#if 0
    /* mark tfo of root as visited */
    std::vector<node> roots = { root };
    while ( true )
    {
      auto it = std::begin( roots );

      /* if we cannot find a root to expand, we are done  */
      if ( it == std::end( roots ) )
        break;

      /* expand roots, i.e., remove the node from roots and add its fanout */
      auto const node = *it;
      roots.erase( it );
      ntk.set_visited( node, divisor_id ); /* mark node as visited */

      ntk.foreach_fanout( node, [&]( auto const &n ){
          assert( ntk.visited( n ) != cover_id );
          if ( ntk.visited( n ) == divisor_id )
            return; /* next */

          /* unique push back */
          if ( std::find( std::begin( roots ), std::end( roots ), n ) == std::end( roots ) )
            roots.push_back( n );
        });

      /* sort the leaves */
      std::sort( std::begin( roots ), std::end( roots ) );
    }
#endif

    std::vector<node> divs;
    collect_divisors( divs, root, leaves );

    // print( root, leaves, divs );

    return { subnetwork{leaves,divs} };
  }

  void print( node const& root, std::vector<node> const& leaves, std::vector<node> const& divs, std::ostream& os = std::cout ) const
  {
    os << "[xcut] r:" << root << " l:{ ";
    for ( const auto& l : leaves )
    {
      os << l << ' ';
    }
    os << "} ";

    os << "d:{ ";
    for ( const auto& d : divs )
    {
      os << d << ' ';
    }
    os << "}";
    os << std::endl;
  }

private:
  void expand_fanin_cut( std::vector<node>& leaves )
  {
    while ( true )
    {
      /* step 1: select a node from the leaves to expand the cut */
      auto it = std::begin( leaves );
      while ( it != std::end( leaves ) )
      {
        /* skip PIs */
        if ( ntk.is_pi( *it ) )
        {
          ++it;
          continue;
        }

        /* skip node if not fanout-free */
        if ( ntk.fanout_size( *it ) > 1 )
        {
          ++it;
          continue;
        }

        if ( cut_size > 0 && leaves.size() - 1 + ntk.fanin_size( *it ) > uint32_t( cut_size ) )
        {
          ++it;
          continue;
        }

        /* found a possible node at which we should expand */
        break;
      }

      /* if we cannot find a leave to expand, we are done with this cut */
      if ( it == std::end( leaves ) )
        return;

      /* step 2: expand the cut, i.e., remove the node from the cut and add its fanin */
      auto const node = *it;
      leaves.erase( it );
      ntk.set_visited( *it, cover_id ); /* mark node as part of the current cone */

      ntk.foreach_fanin( node, [&]( auto const &f ){
          auto const n = ntk.get_node( f );

          /* unique push back */
          if ( std::find( std::begin( leaves ), std::end( leaves ), n ) == std::end( leaves ) )
            leaves.push_back( n );
        });

      /* sort the leaves */
      std::sort( std::begin( leaves ), std::end( leaves ) );
    }
  }

  void collect_divisors( std::vector<node>& divs, node const& root, std::vector<node> const& leaves )
  {
    /* traverse in fanout-direction from leaves and merge nodes with fanin in the current cover */
    for ( const auto r : leaves )
    {
      ntk.foreach_fanout( r, [&]( const auto& d ){
          if ( ntk.visited( d ) == cover_id || ntk.visited( d ) == divisor_id )
            return; /* next */

          if ( ntk.level( d ) > ntk.level( root ) )
          {
            ntk.set_visited( d, divisor_id );
            return; /* next */
          }

          /* check if all fanins are part of the cone */
          bool all_fanins_in_cover = true;
          ntk.foreach_fanin( d, [&]( const auto& f ){
              auto const n = ntk.get_node( f );
              if ( ntk.visited( n ) != cover_id )
              {
                all_fanins_in_cover = false;
                return;
              }
            });
          if ( all_fanins_in_cover )
          {
            divs.push_back( d );
            ntk.set_visited( d, cover_id );
          }
          else
          {
            ntk.set_visited( d, divisor_id );
          }
        });
    }

    std::sort( std::begin( divs ), std::end( divs ) );
  }

private:
  Ntk const& ntk;
  int32_t cut_size;

  uint32_t cover_id{0};
  uint32_t divisor_id{0};
}; /* xcut */

} /* mockturtle */
