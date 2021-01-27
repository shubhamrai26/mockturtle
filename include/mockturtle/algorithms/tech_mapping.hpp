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
  \file tech_mapping.hpp
  \brief tech mapping

  \author Shubham Rai  
*/

#pragma once

#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>


#include <fmt/format.h>

#include "../utils/stopwatch.hpp"
#include "../views/topo_view.hpp"
#include "../io/genlib_reader.hpp"
#include "cut_enumeration.hpp"
#include "cut_enumeration/mf_cut.hpp"

#define ABC_ALLOC(type, num)    ((type *) malloc(sizeof(type) * (num)))

namespace mockturtle
{


struct tech_mapping_params
{
    bool verbose{false};
};


struct tech_mapping_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{0};

  void report() const
  {
    std::cout << fmt::format( "[i] total time = {:>5.2f} secs\n", to_seconds( time_total ) );
  }
};

namespace detail
{
   
//template<class Ntk, bool StoreFunction, typename CutData>
template<class Ntk> 
class tech_mapping_impl
{
public:
    //using network_cuts_t = network_cuts<Ntk, StoreFunction, CutData>;
    //using cut_t = typename network_cuts_t::cut_t;

public:
    tech_mapping_impl(Ntk& ntk, tech_mapping_params const& ps, tech_mapping_stats& st, std::string techlib)
        : ntk(ntk), 
          ps(ps),
          st(st),
          techlib(std::move(techlib))
    {}

    void run()
    {
        stopwatch t(st.time_total);
        /* compute and save topological order */
        top_order.reserve( ntk.size() );
        topo_view<Ntk>( ntk ).foreach_node( [this]( auto n ) {
                top_order.push_back( n );
                } );

        //init_nodes();
        ////print_state();

        //set_mapping_refs<false>();
        ////print_state();

        //while ( iteration < ps.rounds )
        //{
        //    compute_mapping<false>();
        //}

        //while ( iteration < ps.rounds + ps.rounds_ela )
        //{
        //    compute_mapping<true>();
        //}

        //gate_library = read_genlib();

        //compute_expression(gate_library);
        
        //cut enumeration;
        auto cuts = cut_enumeration<Ntk, true>( ntk ); /* true enables truth table computation */
        ntk.foreach_node( [&]( auto n ) {
                const auto index = ntk.node_to_index( n );
                for ( auto const& cut : cuts.cuts( index ) )
                {
                std::cout << "Cut " << *cut
                << " with truth table " << kitty::to_hex( cuts.truth_table( *cut ) )
                << "\n";
                }
                } );

        // mapping each cut to standard cells 


        // Carry out mapping

    }

    
private:
  Ntk& ntk;
  tech_mapping_params const& ps;
  tech_mapping_stats& st;

  uint32_t iteration{0}; /* current mapping iteration */
  uint32_t delay{0};     /* current delay of the mapping */
  uint32_t area{0};      /* current area of the mapping */

  std::vector<gate_struct_t> gate_library;
  std::string techlib;   /* technology library */
  std::vector<node<Ntk>> top_order;
  generic_library<Ntk> glib;

  // make the library and include the concept of phases 
  // create cuts for any network 
  // Map those cuts to the list of gates in the library 
  // Iteratively see which gates to pick so that you have the least area and delay of the overall network

}; /* class tech_mapping_impl */

} /* namespace detail */
/*! \brief tech mapping.
 *
 * This function implements a tech mapping algorithm.  It is controlled by two
 * template arguments `StoreFunction` (defaulted to `true`) and `CutData`
 * (defaulted to `cut_enumeration_mf_cut`).  The first argument `StoreFunction`
 * controls whether the given logic function is stored in the mapping.  In that case
 * truth tables are computed during cut enumeration, which requires more
 * runtime.  The second argument is similar to the `CutData` argument in
 * `cut_enumeration`, which can specialize the cost function to select priority
 * cuts and store additional data.  For tech mapping using this function the
 * type passed as `CutData` must implement the following three fields:
 *
 * - `uint32_t delay`
 * - `float flow`
 * - `float costs`
 *
 * See `include/mockturtle/algorithms/cut_enumeration/mf_cut.hpp` for one
 * example of a CutData type that implements the cost function that is used in
 * the LUT mapper `&mf` in ABC.
 *
 * **Required network functions:**
 * - `size`
 * - `is_pi`
 * - `is_constant`
 * - `node_to_index`
 * - `index_to_node`
 * - `get_node`
 * - `foreach_po`
 * - `foreach_node`
 * - `fanout_size`
 * - `clear_mapping`
 * - `add_to_mapping`
 * - `set_lut_funtion` (if `StoreFunction` is true)
 *
   \verbatim embed:rst

   .. note::

      The implementation of this algorithm was heavily inspired by the tech
      mapping command ``map`` in ABC.
   \endverbatim
 */
//template<class Ntk, bool StoreFunction = false, typename CutData = cut_enumeration_mf_cut>
template<class Ntk>
void tech_mapping( Ntk& ntk, tech_mapping_params const& ps = {}, tech_mapping_stats* pst = nullptr, std::string ifname = "rfet.genlib")
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_index_to_node_v<Ntk>, "Ntk does not implement the index_to_node method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement the fanout_size method" );
  ///static_assert( has_clear_mapping_v<Ntk>, "Ntk does not implement the clear_mapping method" );
  ///static_assert( has_add_to_mapping_v<Ntk>, "Ntk does not implement the add_to_mapping method" );
  //static_assert( !StoreFunction || has_set_cell_function_v<Ntk>, "Ntk does not implement the set_cell_function method" );

  tech_mapping_stats st;
  //detail::tech_mapping_impl<Ntk, StoreFunction, CutData> p( ntk, ps, st );
  detail::tech_mapping_impl<Ntk> p( ntk, ps, st, ifname );
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
