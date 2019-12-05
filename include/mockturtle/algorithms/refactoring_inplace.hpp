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
  \file refactoring_inplace.hpp
  \brief Refactoring inplace

  \author Heinz Riener
*/

#pragma once

#include "../traits.hpp"
#include "../algorithms/node_resynthesis/exact.hpp"
#include "../utils/progress_bar.hpp"
#include "../utils/stopwatch.hpp"
#include "../views/depth_view.hpp"
#include "../views/fanout_view2.hpp"

#include <percy/percy.hpp>

#include <kitty/static_truth_table.hpp>
#include <kitty/constructors.hpp>

namespace mockturtle
{

/*! \brief Parameters for refactoring_inplace.
 *
 * The data structure `refactoring_inplace_params` holds configurable parameters with
 * default arguments for `refactoring_inplace`.
 */
struct refactoring_inplace_params
{
  /*! \brief Maximum number of PIs of reconvergence-driven cuts. */
  uint32_t max_pis{6};

  /*! \brief Maximum number of divisors to consider. */
  uint32_t max_divisors{150};

  /*! \brief Maximum fanout of a node to be considered as root. */
  uint32_t skip_fanout_limit_for_roots{1000};

  /*! \brief Maximum fanout of a node to be considered as divisor. */
  uint32_t skip_fanout_limit_for_divisors{100};

  /*! \brief Show progress. */
  bool progress{false};

  /*! \brief Be verbose. */
  bool verbose{false};

  /*! \brief Allow zero-gain rewriting. */
  bool allow_zero_gain{false};

  /*! \brief Ignore the limit of cuts per node. */
  bool ignore_num_cut_limit{true}; /* if true (all cuts per node are consider) */

  /*! \brief Consider multiple cuts per node. */
  uint32_t num_cuts_per_node{10u};
};

/*! \brief Statistics for refactoring_inplace.
 *
 * The data structure `refactoring_inplace_stats` provides data collected by running
 * `refactoring_inplace`.
 */
struct refactoring_inplace_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{0};

  /*! \brief Accumulated runtime for cut computation. */
  stopwatch<>::duration time_cuts{0};

  /*! \brief Accumulated runtime for cut evaluation/computing a resubsitution. */
  stopwatch<>::duration time_eval{0};

  /*! \brief Accumulated runtime for mffc computation. */
  stopwatch<>::duration time_mffc{0};

  /*! \brief Accumulated runtime for divisor computation. */
  stopwatch<>::duration time_divs{0};

  /*! \brief Accumulated runtime for updating the network. */
  stopwatch<>::duration time_substitute{0};

  /*! \brief Accumulated runtime for simulation. */
  stopwatch<>::duration time_simulation{0};

  /*! \brief Initial network size (before rewriting) */
  uint64_t initial_size{0};

  /*! \brief Total number of divisors  */
  uint64_t num_total_divisors{0};

  /*! \brief Total number of leaves  */
  uint64_t num_total_leaves{0};

  /*! \brief Total number of gain  */
  uint64_t estimated_gain{0};

  uint64_t num_synthesis_timeouts{0};
  uint64_t num_synthesis_successes{0};
  uint64_t cache_hits{0};
  uint64_t cache_misses{0};

  void report() const
  {
    fmt::print( "[i] synthesis success/timeout = {}/{}\n", num_synthesis_successes, num_synthesis_timeouts );
    fmt::print( "[i] cache hits/misses = {}/{}\n", cache_hits, cache_misses );
    fmt::print( "[i] total time                                                  ({:>5.2f} secs)\n", to_seconds( time_total ) );
    fmt::print( "[i]   cut time                                                  ({:>5.2f} secs)\n", to_seconds( time_cuts ) );
    fmt::print( "[i]   mffc time                                                 ({:>5.2f} secs)\n", to_seconds( time_mffc ) );
    fmt::print( "[i]   divs time                                                 ({:>5.2f} secs)\n", to_seconds( time_divs ) );
    fmt::print( "[i]   simulation time                                           ({:>5.2f} secs)\n", to_seconds( time_simulation ) );
    fmt::print( "[i]   evaluation time                                           ({:>5.2f} secs)\n", to_seconds( time_eval ) );
    fmt::print( "[i]   substitute                                                ({:>5.2f} secs)\n", to_seconds( time_substitute ) );
    fmt::print( "[i] total divisors            = {:8d}\n",         ( num_total_divisors ) );
    fmt::print( "[i] total leaves              = {:8d}\n",         ( num_total_leaves ) );
    fmt::print( "[i] estimated gain            = {:8d} ({:>5.2f}%)\n",
                              estimated_gain, ( (100.0 * estimated_gain) / initial_size ) );
  }
};

namespace detail
{

/* based on abcRefs.c */
template<typename Ntk>
class node_mffc_inside
{
public:
  using node = typename Ntk::node;

public:
  explicit node_mffc_inside( Ntk const& ntk )
    : ntk( ntk )
  {
  }

  int32_t run( node const& n, std::vector<node> const& leaves, std::vector<node>& inside )
  {
    /* increment the fanout counters for the leaves */
    for ( const auto& l : leaves )
      ntk.incr_fanout_size( l );

    /* dereference the node */
    auto count1 = node_deref_rec( n );

    /* collect the nodes inside the MFFC */
    node_mffc_cone( n, inside );

    /* reference it back */
    auto count2 = node_ref_rec( n );
    (void)count2;
    assert( count1 == count2 );

    for ( const auto& l : leaves )
      ntk.decr_fanout_size( l );

    return count1;
  }

private:
  /* ! \brief Dereference the node's MFFC */
  int32_t node_deref_rec( node const& n )
  {
    if ( ntk.is_pi( n ) )
      return 0;

    int32_t counter = 1;
    ntk.foreach_fanin( n, [&]( const auto& f ){
        auto const& p = ntk.get_node( f );

        ntk.decr_fanout_size( p );
        if ( ntk.fanout_size( p ) == 0 )
          counter += node_deref_rec( p );
      });

    return counter;
  }

  /* ! \brief Reference the node's MFFC */
  int32_t node_ref_rec( node const& n )
  {
    if ( ntk.is_pi( n ) )
      return 0;

    int32_t counter = 1;
    ntk.foreach_fanin( n, [&]( const auto& f ){
        auto const& p = ntk.get_node( f );

        auto v = ntk.fanout_size( p );
        ntk.incr_fanout_size( p );
        if ( v == 0 )
          counter += node_ref_rec( p );
      });

    return counter;
  }

  void node_mffc_cone_rec( node const& n, std::vector<node>& cone, bool top_most )
  {
    /* skip visited nodes */
    if ( ntk.visited( n ) == ntk.trav_id() )
      return;
    ntk.set_visited( n, ntk.trav_id() );

    if ( !top_most && ( ntk.is_pi( n ) || ntk.fanout_size( n ) > 0 ) )
      return;

    /* recurse on children */
    ntk.foreach_fanin( n, [&]( const auto& f ){
        node_mffc_cone_rec( ntk.get_node( f ), cone, false );
      });

    /* collect the internal nodes */
    cone.emplace_back( n );
  }

  void node_mffc_cone( node const& n, std::vector<node>& cone )
  {
    cone.clear();
    ntk.incr_trav_id();
    node_mffc_cone_rec( n, cone, true );
  }

private:
  Ntk const& ntk;
};

template<typename Ntk, typename TT>
class simulator
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  using truthtable_t = TT;

  explicit simulator( Ntk const& ntk, uint32_t num_divisors, uint32_t max_pis )
    : ntk( ntk )
    , num_divisors( num_divisors )
    , tts( num_divisors + 1 )
    , node_to_index( ntk.size(), 0u )
    , phase( ntk.size(), false )
  {
    auto tt = kitty::create<truthtable_t>( max_pis );
    tts[0] = tt;

    for ( auto i = 0; i < tt.num_vars(); ++i )
    {
      kitty::create_nth_var( tt, i );
      tts[i+1] = tt;
    }
  }

  void resize()
  {
    if ( ntk.size() > node_to_index.size() )
      node_to_index.resize( ntk.size(), 0u );
    if ( ntk.size() > phase.size() )
      phase.resize( ntk.size(), false );
  }

  void assign( node const& n, uint32_t index )
  {
    assert( n < node_to_index.size() );
    assert( index < num_divisors + 1 );
    node_to_index[n] = index;
  }

  truthtable_t get_tt( signal const& s ) const
  {
    auto const tt = tts.at( node_to_index.at( ntk.get_node( s ) ) );
    return ntk.is_complemented( s ) ? ~tt : tt;
  }

  void set_tt( uint32_t index, truthtable_t const& tt )
  {
    tts[index] = tt;
  }

  void normalize( std::vector<node> const& nodes )
  {
    for ( const auto& n : nodes )
    {
      assert( n < phase.size() );
      assert( n < node_to_index.size() );

      if ( n == 0 )
        return;

      auto& tt = tts[node_to_index.at( n )];
      if ( kitty::get_bit( tt, 0 ) )
      {
        tt = ~tt;
        phase[n] = true;
      }
      else
      {
        phase[n] = false;
      }
    }
  }

  bool get_phase( node const& n ) const
  {
    assert( n < phase.size() );
    return phase.at( n );
  }

private:
  Ntk const& ntk;
  uint32_t num_divisors;

  std::vector<truthtable_t> tts;
  std::vector<uint32_t> node_to_index;
  std::vector<bool> phase;
}; /* simulator */

template<class Ntk, class CutCompFn, class RefactoringFn>
class refactoring_inplace_impl
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

  explicit refactoring_inplace_impl( Ntk& ntk, CutCompFn&& cut_comp_fn, RefactoringFn&& refactoring_fn, refactoring_inplace_params const& ps, refactoring_inplace_stats& st )
    : ntk( ntk )
    , sim( ntk, ps.max_divisors, ps.max_pis )
    , cut_comp_fn( cut_comp_fn )
    , refactoring_fn( refactoring_fn )
    , ps( ps )
    , st( st )
  {
    st.initial_size = ntk.num_gates();

    auto const update_level_of_new_node = [&]( const auto& n ){
      ntk.resize_levels();
      update_node_level( n );
    };

    auto const update_level_of_existing_node = [&]( node const& n, const auto& old_children ){
      (void)old_children;
      update_node_level( n );
    };

    auto const update_level_of_deleted_node = [&]( const auto& n ){
      /* update fanout */
      ntk.set_level( n, -1 );
    };

    ntk._events->on_add.emplace_back( update_level_of_new_node );
    ntk._events->on_modified.emplace_back( update_level_of_existing_node );
    ntk._events->on_delete.emplace_back( update_level_of_deleted_node );
  }

  void run()
  {
    stopwatch t( st.time_total );

    progress_bar pbar{ntk.size(), "rewriting |{0}| node = {1:>4}   cand = {2:>4}   est. gain = {3:>5}", ps.progress};

    /* for cost estimation we use reference counters initialized by the fanout size */
    ntk.clear_values();
    ntk.foreach_node( [&]( auto const& n ){
        ntk.set_value( n, ntk.fanout_size( n ) );
      });

    auto const size = ntk.num_gates();
    ntk.foreach_gate( [&]( auto const& n, auto i ){
        if ( i >= size )
          return false; /* terminate */

        if ( ntk.is_dead( n ) )
          return true; /* next */

        /* skip nodes with many fanouts */
        if ( ntk.fanout_size( n ) > ps.skip_fanout_limit_for_roots )
          return true; /* true */

        pbar( i, size - i, candidates, st.estimated_gain, st.num_synthesis_successes, st.num_synthesis_successes + st.num_synthesis_timeouts );

        /* compute a cut for the current node */
        auto const cuts = call_with_stopwatch( st.time_cuts, [&]() {
            return cut_comp_fn( n );
          });

        auto counter = 0u;
        for ( auto const& leaves : cuts )
        {
          if ( !ps.ignore_num_cut_limit && counter++ > ps.num_cuts_per_node )
            return true; /* next node */

          if ( leaves.size() > ps.max_pis )
            continue;

          if ( leaves.size() < 2u || leaves.size() > 15 )
            continue; /* next cut for this node */

          /* evaluate this cut */
          auto const g = call_with_stopwatch( st.time_eval, [&]() {
              return evaluate( n, leaves, /* known functions = */ {} );
            });

          /* if no replacement or self-replacement */
          if ( !g || n >= ntk.get_node( *g ) )
          {
            continue; /* next cut for this node  */
          }

          /* DAG-aware rewriting */
          int32_t value = recursive_deref( n );
          {
            auto [v, contains] = recursive_ref_contains( ntk.get_node( *g ), n );
            recursive_deref( ntk.get_node( *g ) );

            int32_t gain = contains ? -1 : value - v;
            if ( gain > 0 || ( ps.allow_zero_gain && gain == 0 ) )
            {
              /* update progress bar */
              candidates++;
              st.estimated_gain += gain;

              /* update network */
              call_with_stopwatch( st.time_substitute, [&]() {
                  // std::cout << "[i] substitute " << ' ' << n <<  " with " << ntk.get_node( *g ) << std::endl;
                  ntk.substitute_node( n, *g );
                });
            }
          }
          recursive_ref( n );

          return true; /* next node */
        }

        return true; /* next node */
      });
  }

private:
  void simulate( std::vector<node> const &leaves )
  {
    sim.resize();
    for ( auto i = 0u; i < divs.size(); ++i )
    {
      const auto d = divs.at( i );

      /* skip constant 0 */
      if ( d == 0 )
        continue;

      /* assign leaves to variables */
      if ( i < leaves.size() )
      {
        sim.assign( d, i + 1 );
        continue;
      }

      /* compute truth tables of inner nodes */
      sim.assign( d, i - uint32_t( leaves.size() ) + ps.max_pis + 1 );
      std::vector<kitty::dynamic_truth_table> tts;
      ntk.foreach_fanin( d, [&]( const auto& s, auto i ){
          (void)i;
          tts.emplace_back( sim.get_tt( ntk.make_signal( ntk.get_node( s ) ) ) ); /* ignore sign */
        });

      auto const tt = ntk.compute( d, tts.begin(), tts.end() );
      sim.set_tt( i - uint32_t( leaves.size() ) + ps.max_pis + 1, tt );
    }

    /* normalize truth tables */
    sim.normalize( divs );
  }

  std::optional<signal> evaluate( node const& root, std::vector<node> const &leaves, std::vector<node>const& divs )
  {
    uint32_t const required = std::numeric_limits<uint32_t>::max();

    /* collect the MFFC */
    int32_t num_mffc_nodes = call_with_stopwatch( st.time_mffc, [&]() {
        node_mffc_inside collector( ntk );
        auto num_mffc_nodes = collector.run( root, leaves, mffc );
        assert( num_mffc_nodes > 0 );
        return num_mffc_nodes;
      });

    /* collect the divisor nodes in the cut */
    bool div_comp_success = call_with_stopwatch( st.time_divs, [&]() {
        return collect_divisors( root, leaves, required );
      });
    if ( !div_comp_success )
    {
      return std::nullopt;
    }

    assert( num_mffc_nodes > 0u );
    if ( num_mffc_nodes == 1u )
      return std::nullopt; /* next */

    /* update statistics */
    st.num_total_divisors += num_divs;
    st.num_total_leaves += leaves.size();

    /* simulate the collected divisors */
    call_with_stopwatch( st.time_simulation, [&]() { simulate( leaves ); });

    /* get truth table of root */
    auto const tt_root = sim.get_tt( ntk.make_signal( root ) );
    auto const tt_phase = sim.get_phase( root );

    /* trivial cases */
    if ( kitty::is_const0( tt_root ) )
    {
      return ntk.get_constant( false );
    }
    else if ( kitty::is_const0( ~tt_root ) )
    {
      return ntk.get_constant( true );
    }

    std::vector<signal> signal_leaves;
    for ( const auto& l : leaves )
      signal_leaves.emplace_back( ntk.make_signal( l ) );

    std::optional<signal> result = std::nullopt;
    refactoring_fn( ntk, kitty::shrink_to( tt_root, signal_leaves.size() ), std::begin( signal_leaves ), std::end( signal_leaves ),
           [&result]( signal const& s ){
             result = std::make_optional( s );
           } );

    if ( result && tt_phase )
      *result = !*result;

    return result;
  }

  void mark_cone_visited( node const& n )
  {
    /* skip visited nodes */
    if ( ntk.visited( n ) == ntk.trav_id() )
      return;
    ntk.set_visited( n, ntk.trav_id() );

    ntk.foreach_fanin( n, [&]( const auto& f ){
        mark_cone_visited( ntk.get_node( f ) );
      });
  }

  bool collect_divisors( node const& root, std::vector<node> const& leaves, uint32_t required )
  {
    divs.clear();

    /* add the leaves of the cuts to the divisors */
    ntk.incr_trav_id();
    for ( const auto& l : leaves )
    {
      divs.emplace_back( l );
      ntk.set_visited( l, ntk.trav_id() );
    }

    /* mark cone visited (without MFFC) */
    mark_cone_visited( root );

    /* check if the number of divisors is not exceeded */
    if ( divs.size() - leaves.size() + mffc.size() >= ps.max_divisors - ps.max_pis )
      return false;

    /* get the number of divisors to collect */
    int32_t limit = ps.max_divisors - ps.max_pis - ( uint32_t( divs.size() ) + 1 - uint32_t( leaves.size() ) + uint32_t( mffc.size() ) );

    /* explore the fanouts, which are not in the MFFC */
    int32_t counter = 0;
    bool quit = false;

    /* NOTE: this is tricky and cannot be converted to a range-based loop */
    auto size = divs.size();
    for ( auto i = 0u; i < size; ++i )
    {
      auto const d = divs.at( i );

      if ( ntk.fanout_size( d ) > ps.skip_fanout_limit_for_divisors )
        continue;

      /* if the fanout has all fanins in the set, add it */
      ntk.foreach_fanout( d, [&]( node const& p ){
          if ( ntk.visited( p ) == ntk.trav_id() || ntk.level( p ) > required )
            return true; /* next fanout */

          bool all_fanins_visited = true;
          ntk.foreach_fanin( p, [&]( const auto& g ){
              if ( ntk.visited( ntk.get_node( g ) ) != ntk.trav_id() )
              {
                all_fanins_visited = false;
                return false; /* terminate fanin-loop */
              }
              return true; /* next fanin */
            });

          if ( !all_fanins_visited )
            return true; /* next fanout */

          bool has_root_as_child = false;
          ntk.foreach_fanin( p, [&]( const auto& g ){
              if ( ntk.get_node( g ) == root )
              {
                has_root_as_child = true;
                return false; /* terminate fanin-loop */
              }
              return true; /* next fanin */
            });

          if ( has_root_as_child )
            return true; /* next fanout */

          divs.emplace_back( p );
          ++size;
          ntk.set_visited( p, ntk.trav_id() );

          /* quit computing divisors if there are too many of them */
          if ( ++counter == limit )
          {
            quit = true;
            return false; /* terminate fanout-loop */
          }

          return true; /* next fanout */
        });

      if ( quit )
        break;
    }

    /* get the number of divisors */
    num_divs = uint32_t( divs.size() );

    /* add the nodes in the MFFC */
    for ( const auto& t : mffc )
    {
      divs.emplace_back( t );
    }

    assert( root == divs.at( divs.size()-1u ) );
    assert( divs.size() - leaves.size() <= ps.max_divisors - ps.max_pis );

    return true;
  }

private:
    uint32_t recursive_deref( node const& n )
  {
    /* terminate? */
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      return 0;

    /* recursively collect nodes */
    uint32_t value{1u};
    ntk.foreach_fanin( n, [&]( auto const& s ) {
      if ( ntk.decr_value( ntk.get_node( s ) ) == 0 )
      {
        value += recursive_deref( ntk.get_node( s ) );
      }
    } );
    return value;
  }

  uint32_t recursive_ref( node const& n )
  {
    /* terminate? */
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      return 0;

    /* recursively collect nodes */
    uint32_t value{1u};
    ntk.foreach_fanin( n, [&]( auto const& s ) {
      if ( ntk.incr_value( ntk.get_node( s ) ) == 0 )
      {
        value += recursive_ref( ntk.get_node( s ) );
      }
    } );
    return value;
  }

  std::pair<int32_t, bool> recursive_ref_contains( node const& n, node const& repl )
  {
    /* terminate? */
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      return {0, false};

    /* recursively collect nodes */
    int32_t value{1u};
    bool contains = ( n == repl );
    ntk.foreach_fanin( n, [&]( auto const& s ) {
      contains = contains || ( ntk.get_node( s ) == repl );
      if ( ntk.incr_value( ntk.get_node( s ) ) == 0 )
      {
        const auto [v, c] = recursive_ref_contains( ntk.get_node( s ), repl );
        value += v;
        contains = contains || c;
      }
    } );
    return {value, contains};
  }

  void update_node_level( node const& n, bool top_most = true )
  {
    uint32_t curr_level = ntk.level( n );

    uint32_t max_level = 0;
    ntk.foreach_fanin( n, [&]( const auto& f ){
        auto const p = ntk.get_node( f );
        auto const fanin_level = ntk.level( p );
        if ( fanin_level > max_level )
        {
          max_level = fanin_level;
        }
      });
    ++max_level;

    if ( curr_level != max_level )
    {
      ntk.set_level( n, max_level );

      /* update only one more level */
      if ( top_most )
      {
        ntk.foreach_fanout( n, [&]( const auto& p ){
            update_node_level( p, false );
          });
      }
    }
  }

private:
  Ntk& ntk;
  detail::simulator<Ntk, kitty::dynamic_truth_table> sim;
  CutCompFn&& cut_comp_fn;
  RefactoringFn&& refactoring_fn;

  /*! \brief Refactoring parameters */
  refactoring_inplace_params const& ps;

  /*! \brief Refactoring statistics */
  refactoring_inplace_stats& st;

  /* temporary statistics for progress bar */
  uint32_t candidates{0};
  std::vector<node> mffc;
  std::vector<node> divs;
  uint32_t num_divs{0};
};

} /* namespace detail */

/*! \brief Inplace refactoring.
 *
 * **Required network functions:**
 * - `clear_values`
 * - `fanout_size`
 * - `foreach_fanin`
 * - `foreach_gate`
 * - `foreach_node`
 * - `get_constant`
 * - `get_node`
 * - `is_complemented`
 * - `is_pi`
 * - `level`
 * - `make_signal`
 * - `set_value`
 * - `set_visited`
 * - `size`
 * - `substitute_node`
 * - `value`
 * - `visited`
 *
 * \param ntk Input network (will be changed in-place)
 * \param ps Refactoring params
 * \param pst Refactoring statistics
 */
template<class Ntk, class CutCompFn, class RefactoringFn>
void refactoring_inplace( Ntk& ntk, CutCompFn&& cut_comp_fn, RefactoringFn&& refactoring_fn, refactoring_inplace_params const& ps = {}, refactoring_inplace_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_clear_values_v<Ntk>, "Ntk does not implement the clear_values method" );
  static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement the fanout_size method" );
  static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );
  static_assert( has_foreach_gate_v<Ntk>, "Ntk does not implement the foreach_gate method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_get_constant_v<Ntk>, "Ntk does not implement the get_constant method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_make_signal_v<Ntk>, "Ntk does not implement the make_signal method" );
  static_assert( has_set_value_v<Ntk>, "Ntk does not implement the set_value method" );
  static_assert( has_set_visited_v<Ntk>, "Ntk does not implement the set_visited method" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the has_size method" );
  static_assert( has_substitute_node_v<Ntk>, "Ntk does not implement the has substitute_node method" );
  static_assert( has_value_v<Ntk>, "Ntk does not implement the has_value method" );
  static_assert( has_visited_v<Ntk>, "Ntk does not implement the has_visited method" );

  using ntk_view_t = fanout_view2<depth_view<Ntk>>;
  depth_view<Ntk> depth_view{ntk};
  ntk_view_t ntk_view{depth_view};

  refactoring_inplace_stats st;
  {
    detail::refactoring_inplace_impl<ntk_view_t, CutCompFn, RefactoringFn> p( ntk_view, cut_comp_fn, refactoring_fn, ps, st );
    p.run();
    if ( ps.verbose )
    {
      st.report();
    }
  }

  if ( pst )
  {
    *pst = st;
  }
}

} /* namespace mockturtle */
