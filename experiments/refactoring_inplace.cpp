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

#include <mockturtle/algorithms/extract_subnetwork.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
#include <mockturtle/algorithms/node_resynthesis/dsd.hpp>
#include <mockturtle/algorithms/node_resynthesis/composed.hpp>
#include <mockturtle/algorithms/node_resynthesis/cached.hpp>
#include <mockturtle/algorithms/refactoring_inplace.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/views/fanout_view.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/utils/json_utils.hpp>
#include <lorina/aiger.hpp>

#include <experiments.hpp>

inline bool file_exists( std::string const& name )
{
  std::ifstream f(name.c_str());
  return f.good();
}

namespace nlohmann
{
  template <>
  struct adl_serializer<percy::chain>
  {
    static void to_json( json& j, percy::chain const& c )
    {
      j = nlohmann::json{
        {"nr_in", c.nr_in},
        {"fanin", c.fanin},
        {"op_tt_size", c.op_tt_size},
        {"compiled_functions", c.compiled_functions},
        {"steps", c.steps},
        {"operators", c.operators},
        {"outputs", c.outputs},
      };
    }

    static void from_json( nlohmann::json const& j, percy::chain& c )
    {
      if ( !j.is_null() )
      {
        c.nr_in = j.at( "nr_in" ).get<int>();
        c.fanin = j.at( "fanin" ).get<int>();
        c.op_tt_size = j.at( "op_tt_size" ).get<int>();
        c.compiled_functions = j.at( "compiled_functions" ).get<std::vector<kitty::dynamic_truth_table>>();
        c.steps = j.at( "steps" ).get<std::vector<std::vector<int>>>();
        c.operators = j.at( "operators" ).get<std::vector<kitty::dynamic_truth_table>>();
        c.outputs = j.at( "outputs" ).get<std::vector<int>>();
      }
    }
  };
} // namespace nlohmann

namespace nlohmann
{
  template <>
  struct adl_serializer<mockturtle::exact_resynthesis_params>
  {
    static void to_json( json& j, mockturtle::exact_resynthesis_params const& ps )
    {
      j = nlohmann::json{
        {"cache", *ps.cache},
        {"blacklist_cache", *ps.blacklist_cache},
        {"add_alonce_clauses", ps.add_alonce_clauses},
        {"add_colex_clauses", ps.add_colex_clauses},
        {"add_lex_clauses", ps.add_lex_clauses},
        {"add_lex_func_clauses", ps.add_lex_func_clauses},
        {"add_nontriv_clauses", ps.add_nontriv_clauses},
        {"add_noreapply_clauses", ps.add_noreapply_clauses},
        {"add_symvar_clauses", ps.add_symvar_clauses},
        {"conflict_limit", ps.conflict_limit},
        {"solver_type", ps.solver_type},
        {"encoder_type", ps.encoder_type},
        {"synthesis_method", ps.synthesis_method}
      };
    }

    static void from_json( nlohmann::json const& j, mockturtle::exact_resynthesis_params& ps )
    {
      if ( !j.is_null() )
      {
        ps.cache = std::make_shared<mockturtle::exact_resynthesis_params::cache_map_t>( j.at( "cache" ).get<mockturtle::exact_resynthesis_params::cache_map_t>() );
        ps.blacklist_cache = std::make_shared<mockturtle::exact_resynthesis_params::blacklist_cache_map_t>( j.at( "blacklist_cache" ).get<mockturtle::exact_resynthesis_params::blacklist_cache_map_t>() );
        ps.add_alonce_clauses = j.at( "add_alonce_clauses" ).get<bool>();
        ps.add_colex_clauses = j.at( "add_colex_clauses" ).get<bool>();
        ps.add_lex_clauses = j.at( "add_lex_clauses" ).get<bool>();
        ps.add_lex_func_clauses = j.at( "add_lex_func_clauses" ).get<bool>();
        ps.add_nontriv_clauses = j.at( "add_nontriv_clauses" ).get<bool>();
        ps.add_noreapply_clauses = j.at( "add_noreapply_clauses" ).get<bool>();
        ps.add_symvar_clauses = j.at( "add_symvar_clauses" ).get<bool>();
        ps.conflict_limit = j.at( "conflict_limit" ).get<int>();
        ps.solver_type = j.at( "solver_type" ).get<percy::SolverType>();
        ps.encoder_type = j.at( "encoder_type" ).get<percy::EncoderType>();
        ps.synthesis_method = j.at( "synthesis_method" ).get<percy::SynthMethod>();
      }
    }
  };
} // namespace nlohmann

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  using network_type = aig_network;
  
  experiment<std::string, uint32_t, uint32_t, uint32_t, double, float, bool> exp( "refactoring_inplace", "benchmark", "size_before", "size_after", "diff", "diff[%]", "runtime", "equivalent" );

  /* refactoring parameters */
  refactoring_inplace_params ps;
  // ps.progress = true;
  ps.max_pis = 6;

  /* resynthesis function */
  dsd_resynthesis_params dsd_ps;
  dsd_ps.prime_input_limit = 7u;
  dsd_ps.dsd_ps.with_xor = std::is_same_v<network_type, xag_network>;

  auto cexact_resyn = cached_exact_xag_resynthesis<network_type>( "/Users/riener/exact_cache.json", 10e2, *dsd_ps.prime_input_limit );
  dsd_resynthesis<network_type, decltype( cexact_resyn )> dsd_resyn( cexact_resyn, dsd_ps );
  cached_resynthesis<network_type, decltype( dsd_resyn )> cdsd_resyn( dsd_resyn, ps.max_pis, "/Users/riener/dsd_cache.json" );

  cexact_resyn.report();
  cdsd_resyn.report();
  
  /* print parameters of exact resynthesis */
  // std::cout << "[i] conflict limit = " << exact_ps.conflict_limit << std::endl;
  // std::cout << "[i] cache size = " << ( exact_ps.cache ? exact_ps.cache->size() : 0 ) << std::endl;
  // std::cout << "[i] blacklist cache size = " << ( exact_ps.blacklist_cache ? exact_ps.blacklist_cache->size() : 0 ) << std::endl;

  for ( auto const& benchmark : epfl_benchmarks( ~arbiter ) )
  {
    using ntk_view_t = fanout_view2<depth_view<network_type>>;

    fmt::print( "[i] processing {}\n", benchmark );

    network_type aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      fmt::print( "[i] could not read benchmark {}\n", benchmark );
      continue;
    }

    uint32_t const size_before = aig.num_gates();

    refactoring_inplace_stats st;
    depth_view<network_type> depth_aig{aig};
    ntk_view_t ntk_view{depth_aig};

    /* cut computing function */
    xcut<ntk_view_t> cut_comp( ntk_view, ps.max_pis );
    refactoring_inplace( ntk_view, cut_comp, cexact_resyn, ps, &st );
    aig = cleanup_dangling( aig ); // invalidates ntk_view

    uint32_t const size_after = aig.num_gates();

    auto const cec = abc_cec( aig, benchmark );
    exp( benchmark,
         size_before,
         aig.num_gates(),
         size_before - size_after,
         100.0*(size_before - size_after) / size_before,
         to_seconds( st.time_total ),
         cec );

    st.report();
    std::cout << "cec = " << cec << std::endl;

    /* serialize caches after each processed benchmark */
    // nlohmann::json js = exact_ps;
    //
    // std::ofstream ofs( "exact_ps.json" );
    // ofs << std::setw( 2 ) << js << std::endl;
    // ofs.close();
  }

  exp.save();
  exp.table();
  // exp.compare( {}, {}, {"size_after"});

  cexact_resyn.report();
  cdsd_resyn.report();

  return 0;
}
