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

#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/dsd.hpp>
#include <mockturtle/algorithms/refactoring_inplace.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/io/aiger_reader.hpp>
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
  struct adl_serializer<kitty::dynamic_truth_table>
  {
    static void to_json( json& j, kitty::dynamic_truth_table const& tt )
    {
      j = {
        {"_bits", tt._bits},
        {"_num_vars", tt._num_vars}
      };
    }

    static void from_json( json const& j, kitty::dynamic_truth_table& tt )
    {
      if ( !j.is_null() )
      {
        tt._bits = j.at( "_bits" ).get<std::vector<uint64_t>>();
        tt._num_vars = j.at( "_num_vars" ).get<int>();
      }
    }
  };
} // namespace nlohmann

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

template<typename Ntk>
class cut_compute
{
public:
  using node = typename Ntk::node;

public:
  explicit cut_compute()
  {
  }

  std::vector<std::vector<node>> operator()( node const& n )
  {
    (void)n;
    return {{}};
  }
}; /* cut_compute */

/*! \brief Eagerly compute a fanout-free cut into fanin-direction. */
template<typename Ntk>
class mffc_cut
{
public:
  using node = typename Ntk::node;
  using cut = std::vector<node>;

private:
  using cut_iter = typename cut::iterator;

public:
  explicit mffc_cut( Ntk const& ntk, int32_t cut_size = -1 )
    : ntk( ntk )
    , cut_size( cut_size )
  {
  }

  std::vector<cut> operator()( node const& root )
  {
    cut leaves = { root };
    expand_cut( leaves );
    return { leaves };
  }

private:
  void expand_cut( cut& leaves )
  {
    while ( true )
    {
      /* step 1: select a node from the leaves to expand the cut */
      cut_iter it = std::begin( leaves );
      while ( it != std::end( leaves ) )
      {
        /* skip PIs */
        if ( ntk.is_pi( *it ) )
        {
          ++it;
          continue;
        }

        /* skip node if not fanout-free */
        if ( ntk.fanout_size( *it ) > 1 )
        {
          ++it;
          continue;
        }

        if ( cut_size > 0 && leaves.size() - 1 + ntk.fanin_size( *it ) > uint32_t( cut_size ) )
        {
          ++it;
          continue;
        }

        /* found a possible node at which we should expand */
        break;
      }

      /* if we cannot find a leave to expand, we are done with this cut */
      if ( it == std::end( leaves ) )
        return;

      /* step 2: expand the cut, i.e., remove the node from the cut and add its fanin */
      auto const node = *it;
      leaves.erase( it );

      ntk.foreach_fanin( node, [&]( auto const &f ){
          auto const n = ntk.get_node( f );

          /* unique push back */
          if ( std::find( std::begin( leaves ), std::end( leaves ), n ) == std::end( leaves ) )
            leaves.push_back( n );
        });

      /* sort the leaves */
      std::sort( std::begin( leaves ), std::end( leaves ) );
    }
  }

private:
  Ntk const& ntk;
  int32_t cut_size;
}; /* mffc_cut */

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, double, float, bool> exp( "cut_rewriting", "benchmark", "size_before", "size_after", "diff", "diff[%]", "runtime", "equivalent" );

  /* refactoring parameters */
  refactoring_inplace_params ps;
  ps.progress = true;
  ps.max_pis = 10;

  /* exact resynthesis params */
  exact_resynthesis_params exact_ps;
  if ( file_exists( "exact_ps.json" ) )
  {
    std::ifstream ifs( "exact_ps.json" );
    nlohmann::json js;
    ifs >> js;
    exact_ps = js.get<exact_resynthesis_params>();
    ifs.close();
  }
  else
  {
    exact_ps.conflict_limit = 10000;
    exact_ps.cache = std::make_shared<exact_resynthesis_params::cache_map_t>();
    exact_ps.blacklist_cache = std::make_shared<exact_resynthesis_params::blacklist_cache_map_t>();
  }

  /* resynthesis function */
  exact_aig_resynthesis<aig_network> resyn( false, exact_ps );
  dsd_resynthesis<aig_network, decltype( resyn )> dsd_resyn( resyn );

  std::cout << "cache size = " << exact_ps.cache->size() << std::endl;
  std::cout << "blacklist cache size = " << exact_ps.blacklist_cache->size() << std::endl;
  
  for ( auto const& benchmark : epfl_benchmarks() )
  {
    if ( benchmark == "hyp" )
      continue;

    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      fmt::print( "[i] could not read benchmark {}\n", benchmark );
      continue;
    }

    /* cut computing function */
    mffc_cut<aig_network> cut_comp( aig, ps.max_pis );

    uint32_t size_before = aig.num_gates();

    refactoring_inplace_stats st;
    refactoring_inplace( aig, cut_comp, dsd_resyn, ps, &st );
    aig = cleanup_dangling( aig );

    auto const cec = abc_cec( aig, benchmark );
    exp( benchmark,
         size_before,
         aig.num_gates(),
         size_before -  aig.num_gates(),
         100.0*(size_before -  aig.num_gates()) / size_before,
         to_seconds( st.time_total ),
         cec );
  }

  exp.save();
  exp.table();
  // exp.compare( {}, {}, {"size_after"});

  std::cout << "cache size = " << exact_ps.cache->size() << std::endl;
  std::cout << "blacklist size = " << exact_ps.blacklist_cache->size() << std::endl;

  /* store exact resynthesis params */
  nlohmann::json js = exact_ps;

  std::ofstream ofs( "exact_ps.json" );
  ofs << std::setw( 2 ) << js << std::endl;
  ofs.close();
  
  return 0;
}
