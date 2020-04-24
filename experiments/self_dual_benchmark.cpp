#include <utils.hpp>

int main( int argc, char** argv )
{
    //srand(time(NULL));
    srand(5);
    if (argc != 4)
    {
        std::cout << "[e] Usage executable num pis num_levels nodes per level and sd ratio" << std::endl;
        exit(0);
    }
    std::cout << "num_pis "           << argv[1] << std::endl;
    std::cout << "num_levels "        << argv[2] << std::endl;
    std::cout << "nodes per levels "  << argv[3] << std::endl;
    //std::cout << "sd_ratio "          << argv[4] << std::endl;

    uint32_t num_pis                = std::stoi( std::string( argv[1] ) );
    uint32_t num_levels             = std::stoi( std::string( argv[2] ) );
    uint32_t max_nodes_per_levels   = std::stoi( std::string( argv[3] ) );
    uint32_t sd_ratio               = 0; //std::stoi( std::string( argv[4] ) );

    experiments::experiment<std::string, uint32_t, double, double, double, double, double, uint32_t, uint32_t>
        exp( "RFET_area", "benchmark", "init_size", "init_area", "c2rs_area", "dc2_area", "dch_area", "final_area","final_size", "num_pos" );

    for (int i = 0; i < 10; i++)
    {
        sd_ratio++;
        xmg_network xmg;
        create_xmg( xmg, num_pis, num_levels, max_nodes_per_levels, sd_ratio );
        xmg = cleanup_dangling( xmg );

        auto const init_area  = abc_map( xmg, genlib_path );
        auto const c2rs_area  = abc_map_compress2rs( xmg, genlib_path); 
        auto const dch_area   = abc_map_dch( xmg, genlib_path );
        auto const dc2_area   = abc_map_dc2( xmg, genlib_path );

        auto const init_size = xmg.num_gates();

        uint32_t num_iters = 0;
        uint32_t size_per_iteration = 0;
        float total_opt_time = 0; 
        float total_imp = 0;

        do 
        {
            num_iters++;

            auto rw_params = call_rw( xmg );
            auto rs_params = call_rs( xmg );

            total_opt_time = total_opt_time + rs_params.opt_time + rw_params.opt_time;
            size_per_iteration = rs_params.size; 

            if (size_per_iteration == 0u)
                total_imp = 0;
            else
            {
                int diff = size_per_iteration - xmg.num_gates();
                total_imp = 100 * (double(std::abs(diff))/size_per_iteration);
            }
            std::cout << "Iterations # " << num_iters <<  std::endl;

        } while ( total_imp > 0.5 );

        profile(xmg);
        float area_after = abc_map( xmg, genlib_path );

        std::string ofname = "benchmarks_" + std::string( argv[1] ) + "_" + std::string( argv[2] ) + "_" + std::string( argv[3] ) + "_" + std::to_string( sd_ratio ) + ".v";
        mockturtle::write_verilog( xmg, ofname);

        std::cout << "init area "  << init_area     << std::endl;
        std::cout << "c2rs area "  << c2rs_area     << std::endl;
        std::cout << "dch area "   << dch_area      << std::endl;
        std::cout << "dc2 area "   << dc2_area      << std::endl;
        std::cout << "after area " << area_after    << std::endl;
        std::cout << "num_pos "    << xmg.num_pos() << std::endl;
        exp ( ofname, init_size, init_area, c2rs_area, dc2_area, dch_area, area_after, xmg.num_gates(), xmg.num_pos() );
        exp.save();
        exp.table();
    }
    exp.save();
    exp.table();
    return 0;
}
