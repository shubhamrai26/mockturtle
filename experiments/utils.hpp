#include <fmt/format.h>
#include <lorina/lorina.hpp>
#include<mockturtle/networks/xmg.hpp>
#include <mockturtle/properties/xmgcost.hpp>
#include <mockturtle/algorithms/xmg_resub.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/resubstitution.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/xmg_optimization.hpp>
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/algorithms/node_resynthesis/xmg3_npn.hpp>
#include <mockturtle/algorithms/xmg_algebraic_rewriting.hpp>
#include<iostream>
#include <mockturtle/io/write_verilog.hpp>
#include<time.h>
#include <mockturtle/views/depth_view.hpp>


#include <experiments.hpp>
   
using namespace mockturtle; 
using namespace experiments;

std::string const genlib_path = "/home/shubham/My_work/abc-vlsi-cad-flow/std_libs/date_lib_count_tt_4.genlib";


struct opt_parameters
{
    uint32_t size;
    float opt_time;
};

float area_from_abc ( const xmg_network& xmg)
{
    return abc_map( xmg , genlib_path );
}

void profile( const xmg_network& xmg)
{
    xmg_cost_params xmg_ps;
    xmg_ps.reset();
    num_gate_profile( xmg, xmg_ps );
    xmg_ps.report();
    std::cout << "xmg size " << xmg.num_gates() << std::endl;
}

mockturtle::xmg_network::signal create_xmg_sd_node ( xmg_network& xmg, std::vector<std::pair<mockturtle::xmg_network::signal, bool>> sl, const uint32_t& rand_id, const uint32_t& i1, const uint32_t& i2, const uint32_t& i3)
{
    //std::cout << "xmg create self-dual function " << std::endl;
    return (rand_id%2) ? xmg.create_maj(sl[i1].first, sl[i2].first, xmg.create_not(sl[i3].first)) : xmg.create_xor3(xmg.create_not(sl[i1].first), sl[i2].first, sl[i3].first);
}

mockturtle::xmg_network::signal create_xmg_node ( xmg_network& xmg, std::vector<std::pair<mockturtle::xmg_network::signal, bool>> sl, const uint32_t& rand_id, const uint32_t& i1, const uint32_t& i2, const uint32_t& i3 )
{
    //std::cout << "xmg create and " << std::endl;
    if (rand_id%5 == 0)
        return xmg.create_and(sl[i1].first,xmg.create_xor(sl[i3].first, sl[i2].first));
    if (rand_id%5 == 1 )
        return xmg.create_or(sl[i1].first, xmg.create_not(xmg.create_and(sl[i3].first, sl[i1].first)));
    if (rand_id%5 == 2 )
        return xmg.create_or(sl[i1].first, xmg.create_not(xmg.create_and(sl[i3].first, sl[i2].first)));
    if (rand_id%5 == 3 )
        return xmg.create_or(sl[i1].first, xmg.get_constant(1));
    if (rand_id%5 == 4 )
        return xmg.create_xor(sl[i1].first, xmg.create_not(xmg.create_and(sl[i3].first, sl[i2].first)));
        
}

void create_xmg( xmg_network& xmg, const uint32_t& num_pis, const uint32_t& num_lev, const uint32_t& max_nodes_per_levels, const uint32_t sd_ratio)  
{
    using signal = mockturtle::xmg_network::signal;

    std::vector<std::pair<signal, bool>> sl;
    uint32_t sd_attempts     = 0;
    uint32_t normal_attempts = 0;
    bool sd_or_normal = true;

    for (uint32_t i =0; i < num_pis; i++)
    {
        sl.emplace_back( xmg.create_pi( ), false );
    }

    for (uint32_t i =1; i < num_lev; i++)
    {
        for (uint32_t j = 0; j < max_nodes_per_levels; j++)
        {
            auto size_before = xmg.num_gates( ); 
            signal wire;
            uint32_t i1, i2, i3;
            //std::cout << i1 << " " << i2 << " " << i3 << std::endl;
            do 
            {
                i1 =  rand()%sl.size();
                i2 =  rand()%sl.size();
                i3 =  rand()%sl.size();
                //wire =  (rand()%2) ? xmg.create_maj(sl[i1].first, sl[i2].first, sl[i3].first) : xmg.create_xor3(sl[i1].first, sl[i2].first, sl[i3].first);
                //std::cout << "sd or normal " << (sd_or_normal ? "true" : "false" ) << std::endl;
                wire = ( sd_or_normal) ? create_xmg_sd_node( xmg, sl, rand(), i1, i2, i3 ): create_xmg_node(xmg, sl, rand(), i1, i2, i3 );
            } while ( size_before == xmg.num_gates( ) );
            sl.emplace_back( wire, false  );
            sl[i1].second = true;
            sl[i2].second = true;
            sl[i3].second = (sd_or_normal) ? true: false;
            //if ( (sd_attempts < sd_ratio) && (normal_attempts < (10 -sd_ratio) ) )
            if (sd_attempts < sd_ratio) 
            {
                //std::cout <<"sd attempts " << sd_attempts << std::endl;
                sd_or_normal = true;
                sd_attempts++;
                normal_attempts = 0;
            }
            else if (normal_attempts < (10 -sd_ratio) ) 
            {
                //std::cout <<"normal attempts " << normal_attempts << std::endl;
                normal_attempts++;
                sd_or_normal =  false;
            }
            else 
                sd_attempts = 0;
        }
        //std::cout<< "size of sl " <<  sl.size() << std::endl;
        
    }

    for (uint32_t i = 0 ; i < sl.size( );  i++)
    {
        if (sl[i].second == false )
            xmg.create_po( sl[i].first);
    }
    profile(xmg);
}

opt_parameters call_rs( xmg_network& xmg)
{

    // XMG resubstitution parameter set
    resubstitution_params resub_ps;
    resubstitution_stats resub_st;
    resub_ps.max_pis = 8u;
    //resub_ps.progress = true;
    resub_ps.max_inserts = 1u;  
    resub_ps.use_dont_cares = true; 
    resub_ps.window_size = 12u;  

    xmg_resubstitution( xmg, resub_ps, &resub_st);
    xmg = cleanup_dangling( xmg ); 
    
    opt_parameters oparam;
    oparam.size = xmg.num_gates( );
    oparam.opt_time = to_seconds( resub_st.time_total );

    return oparam;
}

opt_parameters call_rw( xmg_network& xmg)
{
    // XMG rewriting parameter set
    cut_rewriting_params cr_ps;
    cut_rewriting_stats cr_st;
    cr_ps.cut_enumeration_ps.cut_size = 4;
    //cr_ps.progress = true;
    
    /* load database from file */
    mockturtle::xmg_network db;
    if ( read_verilog( "xmg3_without_sd.v", mockturtle::verilog_reader( db ) ) != lorina::return_code::success )
    {
        std::cout << "ERROR" << std::endl;
        std::abort();
        return -1;
    }
    else
    {
        std::cout << "[i] DB loaded" << std::endl;
    }

    mockturtle::xmg4_npn_resynthesis<mockturtle::xmg_network> npn_resyn( mockturtle::detail::to_index_list( db ) );
    //mockturtle::xmg3_npn_resynthesis<xmg_network> resyn;
    cut_rewriting( xmg, npn_resyn, cr_ps, &cr_st );
    xmg = cleanup_dangling( xmg );
    
    opt_parameters oparam;
    oparam.size = xmg.num_gates( );
    oparam.opt_time = to_seconds( cr_st.time_total );

    return oparam;
}

