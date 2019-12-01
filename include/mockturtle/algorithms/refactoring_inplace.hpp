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
  \brief Refactoring

  \author Heinz Riener
*/
#pragma once

#include "../utils/progress_bar.hpp"

#include <kitty/dynamic_truth_table.hpp>
#include <fmt/format.h>

namespace mockturtle
{

/*! \brief Replacement strategies for refactoring inplace.
 *
 * The `refactoring_inplace_strategy` defines when to accept a replacement for a cut.
 */
enum refactoring_inplace_replace_strategy
{
  /* \brief Accept first possible cut replacement with gain. */
  first_gain = 0,
  /* \brief Evaluate all possible cut replacements of a fixed set and accept cut with best gain from this set. */
  best_gain_of_fixed_set = 1,
  /* \brief Evaluate all possible cut replacements and accept cut with best gain. */
  best_gain = 2,
};

/*! \brief Parameters for refactoring.
 *
 * The data structure `refactoring_inplace_params` holds configurable parameters with
 * default arguments for `refactoring_inplace`.
 */
struct refactoring_inplace_params
{
  /*! \brief Maximum number of PIs in MFFCs. */
  uint32_t max_pis{6};

  /*! \brief Maximum fanout for a node to be considered as root. */
  uint32_t skip_fanout_limit_for_roots{1000};

  /*! \brief Cut replacement strategy. */
  refactoring_inplace_replace_strategy strategy{best_gain_of_fixed_set};

  /*! \brief Maximum number of cuts per node to be evaluated. */
  uint32_t skip_cut_limit{10};

  /*! \brief Allow zero-gain substitutions */
  bool allow_zero_gain{false};

  /*! \brief Show progress. */
  bool progress{false};

  /*! \brief Be verbose. */
  bool verbose{false};
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

  /*! \brief Accumulated runtime for cut evaluation. */
  stopwatch<>::duration time_eval{0};

  /*! \brief Accumulated runtime for updating the network. */
  stopwatch<>::duration time_substitute{0};

  /*! \brief Number DAG-aware rewriting attempts. */
  uint64_t candidates{0};

  /*! \brief Total gain computed by DAG-aware rewriting. */
  uint64_t estimated_gain{0};

  void report() const
  {
    fmt::print( "[i] total time = ({:>5.2f} secs)\n", to_seconds( time_total ) );
    fmt::print( "[i]   cut time = ({:>5.2f} secs)\n", to_seconds( time_cuts ) );
    fmt::print( "[i]   evaluation time = ({:>5.2f} secs)\n", to_seconds( time_eval ) );
    fmt::print( "[i]   substitute time = ({:>5.2f} secs)\n", to_seconds( time_substitute ) );
  }
};

namespace detail
{

template<class Ntk, class RefactoringFn, class CutComputeFn>
class refactoring_inplace_impl
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  explicit refactoring_inplace_impl( Ntk& ntk, RefactoringFn&& refactoring_fn, CutComputeFn&& cut_compute_fn, refactoring_inplace_params const& ps, refactoring_inplace_stats& st )
    : ntk( ntk )
    , refactoring_fn( refactoring_fn )
    , cut_compute_fn( cut_compute_fn )
    , ps( ps )
    , st( st )
  {
  }

  void run()
  {
    stopwatch t( st.time_total );

    progress_bar pbar{ntk.size(), "rewriting |{0}| node = {1:>4}   cand = {2:>4}   est. gain = {3:>5}", ps.progress};

    /* for cost estimation we use reference counters initialized by the fanout size */
    ntk.clear_values();
    ntk.foreach_node( [&]( auto const& n ) {
      ntk.set_value( n, ntk.fanout_size( n ) );
    } );

    auto const size = ntk.num_gates();
    ntk.foreach_gate( [&]( auto const& n, auto i ){
        if ( i >= size )
          return false; /* terminate */

        if ( ntk.is_dead( n ) )
          return true; /* next */

        /* skip nodes with many fanouts */
        if ( ntk.fanout_size( n ) > ps.skip_fanout_limit_for_roots )
          return true; /* true */

        pbar( i, size - i, st.candidates, st.estimated_gain );

        /* compute a cut for the current node */
        auto const cuts = call_with_stopwatch( st.time_cuts, [&]() {
            return cut_compute_fn( n );
          });

        /* evaluate each cut */
        auto cut_index = 0u;
        for ( const auto& cut : cuts )
        {
          if ( ps.strategy == best_gain_of_fixed_set && cut_index >= ps.skip_cut_limit )
            break;

          ++cut_index;

          /* evaluate this cut */
          auto const g = call_with_stopwatch( st.time_eval, [&]() {
              return evaluate( n, cut );
            });

          /* if no replacement or self-replacement */
          if ( !g || n >= ntk.get_node( *g ) )
            continue; /* next cut for this node  */

          /* DAG-aware rewriting */
          int32_t value = recursive_deref( n );
          auto [v, contains] = recursive_ref_contains( ntk.get_node( *g ), n );
          recursive_deref( ntk.get_node( *g ) );

          int32_t gain = contains ? -1 : value - v;
          if ( gain > 0 || ( ps.allow_zero_gain && gain == 0 ) )
          {
            /* update progress bar */
            st.candidates++;
            st.estimated_gain += gain;

            /* update network */
            call_with_stopwatch( st.time_substitute, [&]() {
                // std::cout << "[i] substitute " << ' ' << n <<  " with " << ntk.get_node( *g ) << std::endl;
                ntk.substitute_node( n, *g );
              });
          }
          recursive_ref( n );
          return true; /* next node */
        }

        return true; /* next node */
      });
  }

private:
  std::optional<signal> evaluate( node const& root, std::vector<node> const &cut )
  {
    (void)root;
    (void)cut;
    return std::optional<signal>{};
  }

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

private:
  Ntk& ntk;
  RefactoringFn&& refactoring_fn;
  CutComputeFn&& cut_compute_fn;
  refactoring_inplace_params const& ps;
  refactoring_inplace_stats& st;
};

} /* namespace detail */

/*! \brief Boolean refactoring.
 *
 * **Required network functions:**
 * - `get_node`
 * - `size`
 * - `make_signal`
 * - `foreach_gate`
 * - `substitute_node`
 * - `clear_visited`
 * - `clear_values`
 * - `fanout_size`
 * - `set_value`
 * - `foreach_node`
 *
 * \param ntk Input network (will be changed in-place)
 * \param refactoring_fn Refactoring function
 * \param ps Refactoring params
 * \param pst Refactoring statistics
 * \param cost_fn Node cost function (a functor with signature `uint32_t(Ntk const&, node<Ntk> const&)`)
 */
template<class Ntk, class RefactoringFn, class CutComputeFn>
void refactoring_inplace( Ntk& ntk, RefactoringFn&& refactoring_fn, CutComputeFn&& cut_compute_fn, refactoring_inplace_params const& ps = {}, refactoring_inplace_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_make_signal_v<Ntk>, "Ntk does not implement the make_signal method" );
  static_assert( has_foreach_gate_v<Ntk>, "Ntk does not implement the foreach_gate method" );
  static_assert( has_substitute_node_v<Ntk>, "Ntk does not implement the substitute_node method" );
  static_assert( has_clear_visited_v<Ntk>, "Ntk does not implement the clear_visited method" );
  static_assert( has_clear_values_v<Ntk>, "Ntk does not implement the clear_values method" );
  static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement the fanout_size method" );
  static_assert( has_set_value_v<Ntk>, "Ntk does not implement the set_value method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );

  refactoring_inplace_stats st;
  detail::refactoring_inplace_impl<Ntk, RefactoringFn, CutComputeFn> p( ntk, refactoring_fn, cut_compute_fn, ps, st );
  p.run();
  if ( ps.verbose )
  {
    st.report();
  }

  if ( pst )
  {
    *pst = st;
  }
}

} /* namespace mockturtle */
