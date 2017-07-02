#include <memory>
#include <gtest/gtest.h>

#include "WCECommon.hpp"

using namespace std;
using namespace FenestrationCommon;

class TestCalcPolynom : public testing::Test {

protected:
  virtual void SetUp() {
  }


};

TEST_F( TestCalcPolynom, Test1 ) {
  SCOPED_TRACE( "Begin Test: Calculate polynom test 1." );

  vector< double > input = { -6.75, 8.65, -0.75 };

  auto poly = Polynom( input );

  EXPECT_NEAR( -10.95, poly.value( 12 ), 1e-6 );

}

TEST_F( TestCalcPolynom, Test2 ) {
  SCOPED_TRACE( "Begin Test: Calculate polynom test 2." );

  vector< double > input = { -6.75, 8.65, -0.75 };

  auto poly = Polynom( input );

  EXPECT_NEAR( 1.15, poly.value( 1 ), 1e-6 );

}

TEST_F( TestCalcPolynom, Test3 ) {
  SCOPED_TRACE( "Begin Test: Calculate polynom test 3." );

  vector< double > input = { -9.27348E-06, 2.288300764, 1.646894009, -15.39761441, 26.12276881,
                             -19.1483186, 5.322076488 };

  auto poly = Polynom( input );

  EXPECT_NEAR( 0.807353444, poly.value( 0.7 ), 1e-6 );

}