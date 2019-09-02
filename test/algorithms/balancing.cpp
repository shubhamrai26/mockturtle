#include <catch.hpp>

#include <mockturtle/algorithms/balancing.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/views/depth_view.hpp>

using namespace mockturtle;

TEST_CASE( "Balancing of path", "[balancing]" )
{
  xag_network xag;
  const auto x1 = xag.create_pi();
  const auto x2 = xag.create_pi();
  const auto x3 = xag.create_pi();
  const auto x4 = xag.create_pi();
  const auto f1 = xag.create_and( x1, x2 );
  const auto f2 = xag.create_and( x3, f1 );
  const auto f3 = xag.create_and( x4, f2 );
  xag.create_po( f3 );

  depth_view dxag{xag};
  balancing( dxag );
}
