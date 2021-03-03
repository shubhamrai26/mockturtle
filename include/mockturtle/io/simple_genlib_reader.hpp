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
  \file genlib_reader.hpp
  \brief generic library reader

  \author Shubham Rai  
*/

#pragma once 

#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <regex>

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/print.hpp>
#include <fmt/format.h>

#include "../utils/stopwatch.hpp"

namespace mockturtle
{

struct generic_library_params
{
    bool verbose{false};

};

struct generic_library_stats
{
    stopwatch<>::duration time_total{0};

    void report() const
    {
        std::cout << fmt::format( "[i] total time = {:>5.2f} secs\n", to_seconds( time_total ) );
    }
};

struct gate_struct_t
{
    std::string name;         /* Name of the gate */
    double area;              /* given area of the gate */
    double delay;             /* given delay of the gate */
    std::string formula;      /* the given formula in SOP format */
    //std::list<T> pins;      /* Total number of pins = input + output */
    std::string out_name;     /* Name of the output pin */
    uint8_t n_inputs;         /* number of inputs */ 
    bool gate0{false};               /* constant 0 gate */ 
    bool gate1{false};               /* constant 1 gate */
    bool gate_inv{false};            /* inverter gate */
    bool universal_gate;      /* To see if you have a gate which is universal */ 
    kitty::dynamic_truth_table tt{6};

    gate_struct_t ()
        :area(0),
        delay(0)
    {
    }

    gate_struct_t(const gate_struct_t &g)
        :name( g.name ),
         area( g.area ),
         delay( g.delay ),
         formula( g.formula ),
         out_name( g.out_name ),
         n_inputs( g.n_inputs ),
         gate0( g.gate0 ),
         gate1( g.gate1 ),
         gate_inv( g.gate_inv ),
         universal_gate( g.universal_gate ),
         tt( g.tt )
        
    {
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

    inline int evaluate_numvars( std::string s)
    {
        std::unordered_map<char, uint32_t> m; 

        for (uint32_t i = 0; i < s.length(); i++) { 
            m[s[i]]++; 
        } 

        return m.size(); 
    }




template<class Ntk>
class generic_library
{
public:
    explicit generic_library( Ntk& ntk, generic_library_params const& ps, generic_library_stats& st, std::string techlib )
        : ntk( ntk ),
          ps( ps ),
          st( st ),
          techlib( std::move( techlib ) )
    {
        static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
        static_assert( has_create_pi_v<Ntk>, "Ntk does not implement the create_pi function" );
        static_assert( has_create_po_v<Ntk>, "Ntk does not implement the create_po function" );
        static_assert( has_create_not_v<Ntk>, "Ntk does not implement the create_not function" );
        static_assert( has_create_nary_and_v<Ntk>, "Ntk does not implement the create_nary_and function" );
        static_assert( has_create_nary_or_v<Ntk>, "Ntk does not implement the create_nary_or function" );
        static_assert( has_create_nary_xor_v<Ntk>, "Ntk does not implement the create_nary_xor function" );
    }

    std::vector<gate_struct_t> run(std::string techlib)
    {
        read_genlib(techlib);
        if(genlib.size() > 0)
        {
            std::cout << "Total gates read = " << genlib.size() << std::endl;
            return genlib;
        }
        else
        {
            std::cout << "No gates in the generic library passed " << techlib;
            return{};

        }


    }

private:
    void read_genlib(std::string techlib)
    {
        std::ifstream inf(techlib);
        if (!inf.is_open())
        {
            std::cout << fmt::format( "[i] Unable to open {} file\n", techlib);
            return;
        }
        std::string line;
        while ( std::getline (inf,line) )
        {
            if (line[0] == '#' || line.empty())
                continue;
            std::cout << line << '\n';
            genlib.emplace_back( populate_gate_entry( line ) );
        }

        // enumerate np configurations for each gate entry 
        for (auto g : genlib)
        {
        }
        
    }



    /*! \brief Read the genlib 
     *
     * Please note that the genlib file should not have any comments 
     */
    gate_struct_t populate_gate_entry(std::string str)
    {
        gate_struct_t g;
        char *token;

        token =  strtok( &str[0], " \t\r\n");
        g.name = strtok( NULL, " \t\r\n");
        g.area = std::stod(strtok( NULL, " \t\r\n" ) );
        g.out_name = chomp( strtok ( NULL, "=" ) );
        g.formula = strtok( NULL, ";" ); 
        g.formula = std::regex_replace(g.formula, std::regex("^ +| +$|( ) +"), "$1");
        
        auto num_vars = evaluate_numvars(g.formula);
        kitty::dynamic_truth_table tt1(num_vars);

        auto res = kitty::create_from_expression(tt1, g.formula);
        g.tt = tt1;

        if(kitty::is_const0(tt1))
            g.gate0 = true;

        if(!(kitty::is_const0(tt1)))
            g.gate1 = true;

        token = strtok( NULL, " \t\r\n" );

        while (token && strcmp(token, "PIN") == 0 )
        {
            token = strtok( NULL, "  \t\r\n" );
            token = strtok( NULL, "  \t\r\n" );
            token = strtok( NULL, "  \t\r\n" );
            // Takin a simple model of delay here with 1 as the constant
            g.delay = std::stod( strtok( NULL, "  \t\r\n" ) ) ;
        }

        // Need to add detect_special_gates
        return g;
    }

private:
    Ntk& ntk; 
    generic_library_stats& st;
    generic_library_params const& ps;
    std::string techlib;
    std::vector<gate_struct_t> genlib;

}; /* class generic_library */

} /* namespace detail */

template<typename Ntk>
std::vector<gate_struct_t> reading_genlib( Ntk& ntk, generic_library_params const& ps = {}, generic_library_stats *gst = nullptr, std::string genlib = "" )
{
    generic_library_stats st;
    mockturtle::detail::generic_library<Ntk> g( ntk, ps, st, genlib );

    if ( gst )
    {
        *gst = st;
    }
    
    return g.run( genlib );
}
    
} /* namespace mockturtle */
