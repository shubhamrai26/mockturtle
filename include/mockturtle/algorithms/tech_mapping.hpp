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

//template<typename TT, typename T>
struct cut_enumeration_techmap_cut
{
  uint32_t delay{0};
  float flow{0};
  float cost{0};
  bool phase{false}; // false for normal and true for complementation
  std::vector<gate_struct_t> gates;
};

struct tech_mapping_params
{
    tech_mapping_params()
    {
        cut_enumeration_ps.cut_size = 6;
        cut_enumeration_ps.cut_limit = 8;
    }

    cut_enumeration_params cut_enumeration_ps{};

    bool verbose{false};

  /*! \brief carrying out area optimization */    
    bool area{true};

  /*! \brief carrying out delay optimization */    
    bool delay{false};
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

/* function to update all cuts after cut enumeration */
template<typename CutData>
struct tech_mapping_update_cuts
{
  template<typename NetworkCuts, typename Ntk>
  static void apply( NetworkCuts const& cuts, Ntk const& ntk )
  {
    (void)cuts;
    (void)ntk;
  }
};



namespace detail
{
   
template<class Ntk, bool StoreFunction, typename CutData>
class tech_mapping_impl
{
public:
    using network_cuts_t = network_cuts<Ntk, StoreFunction, CutData>;
    using cut_t = typename network_cuts_t::cut_t;

public:
    tech_mapping_impl( Ntk& ntk,  std::vector<gate_struct_t> glib,  tech_mapping_params const& ps, tech_mapping_stats& st )
        : ntk( ntk ), 
          gl( glib ),
          ps( ps ),
          st( st ),
          flow_refs( ntk.size() ),
          flows( ntk.size() ),
          delays( ntk.size() ),
          cuts( cut_enumeration<Ntk, StoreFunction, CutData>( ntk, ps.cut_enumeration_ps ) )
  {
    tech_mapping_update_cuts<CutData>().apply( cuts, ntk );
  }

    void run()
    {
        stopwatch t(st.time_total);
        /* compute and save topological order */
        top_order.reserve( ntk.size() );
        topo_view<Ntk>( ntk ).foreach_node( [this]( auto n ) {
                top_order.push_back( n );
                } );

        init_nodes();
   
        //cut enumeration;
        ntk.foreach_node( [&]( auto n ) {
                const auto index = ntk.node_to_index( n );
                for ( auto & cut : cuts.cuts( index ) )
                {
                    auto neg_tt = ~(cuts.truth_table( *cut ));
                    std::cout << "Cut " << *cut << std::endl;
                    //<< " with truth table "; 
                    //kitty::print_binary( cuts.truth_table( *cut ) );
                    //std::cout <<std::endl;

                    //auto num_vars = evaluate_numvars()

                    for (auto const& g: gl)
                    {
                      if (cuts.truth_table( *cut ) == g.tt) 
                      {
                        std::cout << "There is a match" << "with gate" << g.name << std::endl;
                        cuts.cuts( index )[0]->data.delay = g.delay; 
                        cuts.cuts( index )[0]->data.gates.emplace_back(g);
                      }
                      if (~(cuts.truth_table( *cut )) == g.tt) 
                      {
                        std::cout << "There is negative match" << "with gate" << g.name << std::endl;
                        cuts.cuts( index )[0]->data.delay = g.delay; 
                        cuts.cuts( index )[0]->data.gates.emplace_back(g);
                        cuts.cuts( index )[0]->data.phase = true;
                      }
                    }
                    // mapping each cut to standard cells 
                    //map_cut_to_gates(*cut, gl);
                }
                } );



        // Carry out mapping

    }
private:

    void init_nodes()
    {
        ntk.foreach_node( [this]( auto n, auto ) {
                const auto index = ntk.node_to_index( n );

                if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
                {
                /* all terminals have flow 1.0 */
                flow_refs[index] = 1.0f;
                }
                else
                {
                flow_refs[index] = static_cast<float>( ntk.fanout_size( n ) );
                }

                flows[index] = cuts.cuts( index )[0]->data.flow;
                delays[index] = cuts.cuts( index )[0]->data.delay;
                //phase[index] = cuts.cuts( index )[0]->data.phase;
                } );

    }

    void map_cut_to_gates (cut_t const& cut, std::vector<gate_struct_t> gl)
    {
    }

    
private:
  Ntk& ntk;
  tech_mapping_params const& ps;
  tech_mapping_stats& st;

  uint32_t iteration{0}; /* current mapping iteration */
  uint32_t delay{0};     /* current delay of the mapping */
  uint32_t area{0};      /* current area of the mapping */

  std::vector<gate_struct_t> const gl; 
  std::vector<node<Ntk>> top_order;
  std::vector<float> flow_refs;
  std::vector<uint32_t> map_refs;
  std::vector<float> flows;
  std::vector<uint32_t> delays;
  network_cuts_t cuts;



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
template<class Ntk, bool StoreFunction = true, typename CutData = cut_enumeration_techmap_cut>
//template<class Ntk>
void tech_mapping( Ntk& ntk,  std::vector<gate_struct_t> const& g, tech_mapping_params const& ps = {}, tech_mapping_stats* pst = nullptr)
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
  static_assert( !StoreFunction || has_set_cell_function_v<Ntk>, "Ntk does not implement the set_cell_function method" );

  tech_mapping_stats st;
  detail::tech_mapping_impl<Ntk, StoreFunction, CutData> p( ntk, g, ps, st );
  //detail::tech_mapping_impl<Ntk> p( ntk, ps, st, g);
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
