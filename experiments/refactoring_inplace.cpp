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
  explicit mffc_cut( Ntk const& ntk )
    : ntk( ntk )
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
  
  auto cache = std::make_shared<exact_resynthesis_params::cache_map_t>();
  auto blacklist_cache = std::make_shared<exact_resynthesis_params::blacklist_cache_map_t>();
  
  /* resynthesis function */
  exact_resynthesis_params exact_ps;
  exact_ps.conflict_limit = 10000;
  exact_ps.cache = cache;
  exact_ps.blacklist_cache = blacklist_cache;
  exact_aig_resynthesis<aig_network> resyn( false, exact_ps );
  
  dsd_resynthesis<aig_network, decltype( resyn )> dsd_resyn( resyn );
  
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
    mffc_cut<aig_network> cut_comp( aig );

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

  std::cout << "cache size = " << cache->size() << std::endl;
  std::cout << "blacklist size = " << blacklist_cache->size() << std::endl;
  return 0;
}
