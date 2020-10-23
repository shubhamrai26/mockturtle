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
   
    
std::string find_formula(std::string a)
{
    uint32_t lp = a.find( "O=" );
    return a.substr(lp + 2);

}

inline char * chomp( char *s )
{
    char *a, *b, *c;
    // remove leading spaces
    for ( b = s; *b; b++ )
        if ( !isspace(*b) )
            break;
    // strsav the string
    a = strcpy( ABC_ALLOC(char, strlen(b)+1), b );
    // remove trailing spaces
    for ( c = a+strlen(a); c > a; c-- )
        if ( *c == 0 || isspace(*c) )
            *c = 0;
        else
            break;
    return a;
}




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

    struct gate_struct_t
    {
        std::string name;         /* Name of the gate */
        double area;              /* given area of the gate */
        double delay;             /* given delay of the gate */
        std::string formula;      /* the given formula in SOP format */
        //std::list<T> pins;      /* Total number of pins = input + output */
        std::string out_name;     /* Name of the output pin */
        uint8_t n_inputs;         /* number of inputs */ 
        bool gate0;               /* constant 0 gate */ 
        bool gate1;               /* constant 1 gate */
        bool gate_inv;            /* inverter gate */
        bool universal_gate;      /* To see if you have a gate which is universal */ 

        gate_struct_t ()
            :area(0),
            delay(0)
        {
        }

        gate_struct_t(const gate_struct_t &g)
        {
            name = std::move( g.name );
            area = g.area;
            delay = g.delay;
            formula = std::move( g.formula );
            out_name = std::move( g.out_name );
        }

    };

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

        gate_library = read_genlib();
        //derive_mapping();
    }

    std::vector<gate_struct_t> read_genlib()
    {
        std::ifstream inf(techlib);
        if (!inf.is_open())
        {
            std::cout << fmt::format( "[i] Unable to open {} file\n", techlib);
            return {};
        }
        std::vector<gate_struct_t> gate_library;
        std::string line;
        while ( std::getline (inf,line) )
        {
            std::cout << line << '\n';
            gate_library.emplace_back( populate_gate_entry( line ) );
        }
        std::cout << "Total gates read = " << gate_library.size() << std::endl;
        return gate_library;
    }

    /*! \brief Read the genlib 
     *
     * Please note that the genlib file should not have any comments 
     */
    gate_struct_t populate_gate_entry(std::string str)
    {
        //if (str.empty())
        //    return {};

        gate_struct_t g;
        char *token;

        token =  strtok( &str[0], " \t\r\n");
        g.name = strtok( NULL, " \t\r\n");
        g.area = std::stod(strtok( NULL, " \t\r\n" ) );
        g.out_name = chomp( strtok ( NULL, "=" ) );
        g.formula = strtok( NULL, ";" ); 
        token = strtok( NULL, " \t\r\n" );

        while (token && strcmp(token, "PIN") == 0 )
        {
            token = strtok( NULL, "  \t\r\n" );
            token = strtok( NULL, "  \t\r\n" );
            token = strtok( NULL, "  \t\r\n" );
            // Takin a simple model of delay here with 1 as the constant
            g.delay = std::stod( strtok( NULL, "  \t\r\n" ) ) ;
        }
        return g;
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

      The implementation of this algorithm was heavily inspired but the LUT
      mapping command ``&mf`` in ABC.
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
