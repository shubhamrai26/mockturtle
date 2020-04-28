#include <mockturtle/algorithms/resubstitution.hpp>
#include <mockturtle/algorithms/xmg_resub.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/collapse_mapped.hpp>
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/algorithms/detail/database_generator.hpp>
#include <mockturtle/algorithms/lut_mapping.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
#include <mockturtle/algorithms/node_resynthesis/cached.hpp>
#include <mockturtle/algorithms/node_resynthesis/xmg4_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/xmg3_npn.hpp>
#include <mockturtle/properties/xmgcost.hpp>
#include <mockturtle/algorithms/node_resynthesis/xmg_npn.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/index_list.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/views/mapping_view.hpp>
#include <lorina/lorina.hpp>
#include <utils.hpp>

#include <fmt/format.h>

#include <string>
#include <vector>


void experiment_compare()
{
    using namespace experiments;
    using namespace mockturtle;

    /* Winston & Mathias's results from ASP-DAC'17 */
    std::map<std::string, std::tuple<uint32_t, uint32_t>> aspdac17_xmg;
    aspdac17_xmg.emplace( "adder", std::make_tuple( 639, 251 ) );
    aspdac17_xmg.emplace( "bar", std::make_tuple( 3281, 888 ) );
    aspdac17_xmg.emplace( "div", std::make_tuple( 29607, 12094 ) );
    aspdac17_xmg.emplace( "hyp", std::make_tuple( 155349, 50835 ) );
    aspdac17_xmg.emplace( "log2", std::make_tuple( 27936, 8438 ) );
    aspdac17_xmg.emplace( "max", std::make_tuple( 2296, 745 ) );
    aspdac17_xmg.emplace( "multiplier", std::make_tuple( 17508, 5700 ));
    aspdac17_xmg.emplace( "sin", std::make_tuple( 5100, 1655 ) );
    aspdac17_xmg.emplace( "sqrt", std::make_tuple( 20130, 6595 ) );
    aspdac17_xmg.emplace( "square", std::make_tuple( 15070, 3969 ) );
    aspdac17_xmg.emplace( "arbiter", std::make_tuple( 10621, 3752 ) );
    aspdac17_xmg.emplace( "cavlc", std::make_tuple( 706, 139 ) );
    aspdac17_xmg.emplace( "ctrl", std::make_tuple( 116, 29 ) );
    aspdac17_xmg.emplace( "i2c", std::make_tuple( 1264, 372 ) );
    aspdac17_xmg.emplace( "int2float", std::make_tuple( 245, 56 ) );
    aspdac17_xmg.emplace( "mem_ctrl", std::make_tuple( 42019, 12736 ) );
    aspdac17_xmg.emplace( "priority", std::make_tuple( 750, 233 ) );
    aspdac17_xmg.emplace( "router", std::make_tuple( 212, 97 ) );
    aspdac17_xmg.emplace( "voter", std::make_tuple( 6737, 2163 ) );


    uint32_t const size = 6;
    exact_xmg_resynthesis_params xmg2_ps;
    xmg2_ps.use_xor3 = false;
    exact_xmg_resynthesis<mockturtle::xmg_network> xmg2_exact( xmg2_ps );
    cached_resynthesis<mockturtle::xmg_network, decltype( xmg2_exact )> cached_xmg2_exact( xmg2_exact, size, "exact_xmg2_cache6.v" ); 

    experiments::experiment<std::string, uint32_t, uint32_t, std::string, uint32_t, uint32_t, std::string>
        exp( "cut_rewriting", "benchmark", "size_aspdac", "size_ours", "xmg_improv", "klut6_aspdac", "klut6_ours", "klut_improv" );
    for ( auto const& benchmark : epfl_benchmarks( ) )
    {
        fmt::print( "[i] processing {}\n", benchmark );

        /* read aig */
        mockturtle::aig_network aig;
        if ( lorina::read_aiger( experiments::benchmark_path( benchmark ), mockturtle::aiger_reader( aig ) ) != lorina::return_code::success )
        {
            std::cout << "ERROR 2" << std::endl;
            std::abort();
            return;
        }

        /* LUT map AIG into k-LUT network */
        auto klut = lut_map( aig, size );

        /* translate to XMG */
        xmg_network xmg;

        while ( true )
        {
            auto const klut_size_before = klut.size();
            xmg = mockturtle::node_resynthesis<mockturtle::xmg_network>( klut, cached_xmg2_exact );
            //xmg = mockturtle::node_resynthesis<mockturtle::xmg_network>( klut, npn_resyn );

            /* resubstitution */ 
            //mockturtle::resubstitution_params resub_ps;
            //mockturtle::resubstitution_stats resub_st;
            //resub_ps.max_pis = 8u;
            ////resub_ps.progress = true;
            //resub_ps.max_inserts = 1u;  
            //resub_ps.use_dont_cares = true; 
            //resub_ps.window_size = 12u;  
            //mockturtle::xmg_resubstitution(xmg, resub_ps, &resub_st);
            //xmg = mockturtle::cleanup_dangling( xmg );

            /* option 3: ABC if-mapping flow */
            //auto const new_klut = lut_map( xmg );

            /* option 4: mockturtle mf-mapping flow */
            mockturtle::lut_mapping_params ps;
            ps.cut_enumeration_ps.cut_size = 6;
            ps.cut_enumeration_ps.cut_limit = 16;

            mockturtle::mapping_view<mockturtle::xmg_network, true> mapped_xmg{xmg};
            mockturtle::lut_mapping<decltype( mapped_xmg ), true>( mapped_xmg, ps );
            const auto new_klut = *mockturtle::collapse_mapped_network<mockturtle::klut_network>( mapped_xmg );

            if ( new_klut.size() >= klut_size_before )
                break;

            klut = new_klut;
        }
        xmg_cost_params ps1;
        ps1.reset();
        num_gate_profile( xmg, ps1);

        mockturtle::lut_mapping_params ps;
        ps.cut_enumeration_ps.cut_size = 6;
        ps.cut_enumeration_ps.cut_limit = 16;

        mockturtle::mapping_view<mockturtle::xmg_network, true> mapped_xmg{xmg};
        mockturtle::lut_mapping<decltype( mapped_xmg ), true>( mapped_xmg, ps );
        const auto new_klut6 = *mockturtle::collapse_mapped_network<mockturtle::klut_network>( mapped_xmg );


        std::cout << "final XMG size = " << xmg.size() << std::endl;
        std::cout << "final KLUT-6 size = " << new_klut6.size() << std::endl;
        exp( benchmark,
                std::get<0>( aspdac17_xmg[benchmark] ), xmg.size(),
                fmt::format( "{:3.2f}", (double(std::get<0>( aspdac17_xmg[benchmark] ) ) - xmg.size())/std::get<0>( aspdac17_xmg[benchmark] ) ),
                std::get<1>( aspdac17_xmg[benchmark] ), new_klut6.size(),
                fmt::format( "{:3.2f}", (double(std::get<1>( aspdac17_xmg[benchmark] ) ) - new_klut6.size())/std::get<1>( aspdac17_xmg[benchmark] ) ) );
        exp.save();
        exp.table();
    }

    exp.save();
    exp.table();
}


int main()
{
    experiment_compare();
    return 0;
}

