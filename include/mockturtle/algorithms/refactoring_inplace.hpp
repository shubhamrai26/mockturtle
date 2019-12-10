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

#include <kitty/dynamic_truth_table.hpp>
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

template<class Ntk, class CutCompFn, class RefactoringFn>
class refactoring_inplace_impl
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

  explicit refactoring_inplace_impl( Ntk& ntk, CutCompFn&& cut_comp_fn, RefactoringFn&& refactoring_fn, refactoring_inplace_params const& ps, refactoring_inplace_stats& st )
    : ntk( ntk )
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
        auto /* const */ subnetworks = call_with_stopwatch( st.time_cuts, [&]() {
            return cut_comp_fn( n );
          });

        auto counter = 0u;
        for ( auto const& subntk : subnetworks )
        {
          if ( !ps.ignore_num_cut_limit && counter++ > ps.num_cuts_per_node )
            return true; /* next node */

          if ( subntk.leaves.size() > ps.max_pis )
            continue;

          if ( subntk.leaves.size() < 2u || subntk.leaves.size() > 15 )
          {
            continue; /* next cut for this node */
          }

          /* evaluate this cut */
          auto const g = call_with_stopwatch( st.time_eval, [&]() {
              return evaluate( n, subntk );
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
                  // auto const tt1 = simulate_subnetwork( ntk.make_signal( n ), subnetworks[0u].leaves, subnetworks[0u].divs );

                  // std::cout << "[i] substitute " << ' ' << n <<  " with " << ntk.get_node( *g ) << ':' << ntk.is_complemented( *g ) << std::endl;
                  // cut_comp_fn.print( n, subnetworks[0u].leaves, subnetworks[0u].divs );

                  ntk.substitute_node( n, *g );
                  // auto const tt2 = simulate_subnetwork( *g, subnetworks[0u].leaves, subnetworks[0u].divs );

                  /* ensure local equivalence */
                  // assert( tt1 == tt2 );
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
  kitty::dynamic_truth_table simulate_subnetwork( signal const& s, std::vector<node> const& leaves, std::vector<node> const& divs ) const
  {
    cut_view<Ntk> cutv( ntk, leaves, s, divs );

    unordered_node_map<kitty::dynamic_truth_table,cut_view<Ntk>> values( cutv );
    default_simulator<kitty::dynamic_truth_table> simulator( leaves.size() );
    simulate_nodes<kitty::dynamic_truth_table,cut_view<Ntk>>( cutv, values, simulator );

    return cutv.is_complemented( s ) ? ~values[s] : values[s];
  }

  template<class SubNtk>
  std::optional<signal> evaluate( node const& root, SubNtk const& subntk )
  {
    /* collect the MFFC */
    int32_t num_mffc_nodes = call_with_stopwatch( st.time_mffc, [&]() {
        node_mffc_inside collector( ntk );
        auto num_mffc_nodes = collector.run( root, subntk.leaves, mffc );
        assert( num_mffc_nodes > 0 );
        return num_mffc_nodes;
      });

    assert( num_mffc_nodes > 0 );
    if ( num_mffc_nodes == 1u )
      return std::nullopt; /* next */

#if 0
    /* update statistics */
    st.num_total_divisors += num_divs;
    st.num_total_leaves += subntk.leaves.size();

    /* simulate the collected divisors */
    call_with_stopwatch( st.time_simulation, [&]() { simulate( subntk.leaves ); });

    /* get truth table of root */
    auto const tt_root = sim.get_tt( ntk.make_signal( root ) );
    auto const tt_phase = sim.get_phase( root );
#endif

    cut_view<Ntk> cutv( ntk, subntk.leaves, ntk.make_signal( root ), subntk.divs );
    unordered_node_map<kitty::dynamic_truth_table,cut_view<Ntk>> values( cutv );
    default_simulator<kitty::dynamic_truth_table> simulator( subntk.leaves.size() );
    call_with_stopwatch( st.time_simulation, [&]() {
        simulate_nodes<kitty::dynamic_truth_table,cut_view<Ntk>>( cutv, values, simulator );
      } );

    auto const tt_root = values[root];

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
    for ( const auto& l : subntk.leaves )
      signal_leaves.emplace_back( ntk.make_signal( l ) );

    /* filter divisors by their functions */
    std::vector<std::pair<signal, kitty::dynamic_truth_table>> filtered_divs;
    for ( const auto& d : subntk.divs )
    {
      /* avoid adding divisors with trivial or pure variable functions */
      auto const tt = kitty::shrink_to( values[d]/*sim.get_tt( ntk.make_signal( d ) )*/, signal_leaves.size() );
      if ( kitty::is_const0( tt ) || kitty::is_const0( ~tt ) )
        continue;
      bool next = false;
      for ( auto i = 0u; i < signal_leaves.size(); ++i )
      {
        kitty::dynamic_truth_table nvar( signal_leaves.size() );
        kitty::create_nth_var( nvar, i );
        if ( tt == nvar || ~tt == nvar )
        {
          next = true;
          break;
        }
      }
      if ( next )
        continue;

      /* ensure that all functions are pairwise different (and pairwise different to their complements) */
      if ( std::find_if( std::begin( filtered_divs ), std::end( filtered_divs ),
                      [&tt]( const auto& p ){ return p.second == tt || p.second == ~tt; } ) == std::end( filtered_divs ) )
      {
        filtered_divs.emplace_back( ntk.make_signal( d ), tt );
      }

      /* TODO: in case of equal functions, we should chose the best one */
    }

    /* add divisor functions to synthesis problem */
    refactoring_fn.clear_functions();
    for ( const auto& d : filtered_divs )
    {
      /* normalize d if necessary */
      if ( !kitty::is_normal( d.second ) )
        refactoring_fn.add_function( !d.first, ~d.second );
      else
        refactoring_fn.add_function( d.first, d.second );
    }

    std::optional<signal> result = std::nullopt;
    refactoring_fn( ntk, kitty::shrink_to( tt_root, signal_leaves.size() ), std::begin( signal_leaves ), std::end( signal_leaves ),
                    [&result]( signal const& s ){
                      result = std::make_optional( s );
                    } );

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

  // using ntk_view_t = fanout_view2<depth_view<Ntk>>;
  // depth_view<Ntk> depth_view{ntk};
  // ntk_view_t ntk_view{depth_view};

  refactoring_inplace_stats st;
  {
    detail::refactoring_inplace_impl<Ntk, CutCompFn, RefactoringFn> p( ntk, cut_comp_fn, refactoring_fn, ps, st );
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
