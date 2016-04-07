#include <memory>
#include <gtest/gtest.h>

#include "UniformDiffuseBSDFLayer.hpp"
#include "WovenCell.hpp"
#include "WovenCellDescription.hpp"
#include "MaterialDescription.hpp"
#include "BSDFDirections.hpp"
#include "SquareMatrix.hpp"
#include "BSDFResults.hpp"
#include "FenestrationCommon.hpp"

using namespace std;
using namespace LayerOptics;
using namespace FenestrationCommon;

class TestWovenShadeUniformMaterial : public testing::Test {

private:
  shared_ptr< CBSDFLayer > m_Shade;

protected:
  virtual void SetUp() {
    // create material
    double Tmat = 0.2;
    double Rfmat = 0.75;
    double Rbmat = 0.75;
    double minLambda = 0.3;
    double maxLambda = 2.5;
    shared_ptr< CMaterialSingleBand > aMaterial = 
      make_shared< CMaterialSingleBand >( Tmat, Tmat, Rfmat, Rbmat, minLambda, maxLambda );

    // make cell geometry
    double diameter = 6.35; // mm
    double spacing = 19.05; // mm
    shared_ptr< CCellDescription > aCellDescription = 
      make_shared< CWovenCellDescription >( diameter, spacing );

    shared_ptr< CBSDFHemisphere > aBSDF = make_shared< CBSDFHemisphere >( BSDFBasis::Quarter );

    shared_ptr< CUniformDiffuseCell > aCell = make_shared< CWovenCell >( aMaterial, aCellDescription );
    
    m_Shade = make_shared< CUniformDiffuseBSDFLayer >( aCell, aBSDF );

  };

public:
  shared_ptr< CBSDFLayer > GetShade() { return m_Shade; };

};

TEST_F( TestWovenShadeUniformMaterial, TestSolarProperties ) {
  SCOPED_TRACE( "Begin Test: Woven shade uniform material." );
  
  shared_ptr< CBSDFLayer > aShade = GetShade();

  shared_ptr< CBSDFResults > aResults = aShade->getResults();

  double tauDiff = aResults->TauDiff( Side::Front );
  EXPECT_NEAR( 0.467578877, tauDiff, 1e-6 );

  double RfDiff = aResults->RhoDiff( Side::Front );
  EXPECT_NEAR( 0.496269069, RfDiff, 1e-6 );

  double RbDiff = aResults->RhoDiff( Side::Back );
  EXPECT_NEAR( 0.496269069, RbDiff, 1e-6 );

  shared_ptr< CSquareMatrix > aT = aResults->Tau( Side::Front );

  size_t size = aT->getSize();

  // Test diagonal
  vector< double > correctResults;
  correctResults.push_back( 5.8208903745997009 );
  correctResults.push_back( 6.1112742537069362 );
  correctResults.push_back( 6.1105232559211569 );
  correctResults.push_back( 6.0948781743621048 );
  correctResults.push_back( 6.1105232559211569 );
  correctResults.push_back( 6.1112742537069362 );
  correctResults.push_back( 6.1105232559211569 );
  correctResults.push_back( 6.0948781743621048 );
  correctResults.push_back( 6.1105232559211569 );
  correctResults.push_back( 5.1378571879714912 );
  correctResults.push_back( 5.1628657464653953 );
  correctResults.push_back( 5.0391653915946177 );
  correctResults.push_back( 4.8274726056748234 );
  correctResults.push_back( 5.0391653915946177 );
  correctResults.push_back( 5.1628657464653953 );
  correctResults.push_back( 5.1378571879714912 );
  correctResults.push_back( 5.1628657464653953 );
  correctResults.push_back( 5.0391653915946186 );
  correctResults.push_back( 4.8274726056748234 );
  correctResults.push_back( 5.0391653915946177 );
  correctResults.push_back( 5.1628657464653962 );
  correctResults.push_back( 3.8105833651825955 );
  correctResults.push_back( 4.0021223385767088 );
  correctResults.push_back( 3.4180170235700418 );
  correctResults.push_back( 1.3822013523549581 );
  correctResults.push_back( 3.4180170235700409 );
  correctResults.push_back( 4.0021223385767088 );
  correctResults.push_back( 3.8105833651825955 );
  correctResults.push_back( 4.0021223385767088 );
  correctResults.push_back( 3.4180170235700431 );
  correctResults.push_back( 1.3822013523549581 );
  correctResults.push_back( 3.4180170235700418 );
  correctResults.push_back( 4.0021223385767106 );
  correctResults.push_back( 0.1035160170829772 );
  correctResults.push_back( 0.0936502628857958 );
  correctResults.push_back( 0.1035160170829772 );
  correctResults.push_back( 0.0936502628857958 );
  correctResults.push_back( 0.1035160170829772 );
  correctResults.push_back( 0.0936502628857958 );
  correctResults.push_back( 0.1035160170829772 );
  correctResults.push_back( 0.0936502628857958 );

  vector< double > calculatedResults;
  for( size_t i = 0; i < size; ++i ) {
    calculatedResults.push_back( (*aT)[i][i] );
  }

  EXPECT_EQ( correctResults.size(), calculatedResults.size() );
  for( size_t i = 0; i < size; ++i ) {
    EXPECT_NEAR( correctResults[i], calculatedResults[i], 1e-5 );
  }

  // Test first row
  correctResults.clear();
  correctResults.push_back( 5.82089   );
  correctResults.push_back( 0.0406216 );
  correctResults.push_back( 0.0406251 );
  correctResults.push_back( 0.0406964 );
  correctResults.push_back( 0.0406251 );
  correctResults.push_back( 0.0406216 );
  correctResults.push_back( 0.0406251 );
  correctResults.push_back( 0.0406964 );
  correctResults.push_back( 0.0406251 );
  correctResults.push_back( 0.0432346 );
  correctResults.push_back( 0.0431121 );
  correctResults.push_back( 0.0437211 );
  correctResults.push_back( 0.0447624 );
  correctResults.push_back( 0.0437211 );
  correctResults.push_back( 0.0431121 );
  correctResults.push_back( 0.0432346 );
  correctResults.push_back( 0.0431121 );
  correctResults.push_back( 0.0437211 );
  correctResults.push_back( 0.0447624 );
  correctResults.push_back( 0.0437211 );
  correctResults.push_back( 0.0431121 );
  correctResults.push_back( 0.0596576 );
  correctResults.push_back( 0.0623422 );
  correctResults.push_back( 0.0678335 );
  correctResults.push_back( 0.0716109 );
  correctResults.push_back( 0.0678335 );
  correctResults.push_back( 0.0623422 );
  correctResults.push_back( 0.0596576 );
  correctResults.push_back( 0.0623422 );
  correctResults.push_back( 0.0678335 );
  correctResults.push_back( 0.0716109 );
  correctResults.push_back( 0.0678335 );
  correctResults.push_back( 0.0623422 );
  correctResults.push_back( 0.103516  );
  correctResults.push_back( 0.0936503 );
  correctResults.push_back( 0.103516  );
  correctResults.push_back( 0.0936503 );
  correctResults.push_back( 0.103516  );
  correctResults.push_back( 0.0936503 );
  correctResults.push_back( 0.103516  );
  correctResults.push_back( 0.0936503 );

  calculatedResults.clear();
  for( size_t i = 0; i < size; ++i ) {
    calculatedResults.push_back( (*aT)[ 0 ][ i ] );
  }

  EXPECT_EQ( correctResults.size(), calculatedResults.size() );
  for( size_t i = 0; i < size; ++i ) {
    EXPECT_NEAR( correctResults[i], calculatedResults[i], 1e-5 );
  }

  // Test first row for reflectance matrix
  shared_ptr< CSquareMatrix > aRf = aResults->Rho( Side::Front );

  correctResults.clear();
  correctResults.push_back( 0.128103 );
  correctResults.push_back( 0.130833 );
  correctResults.push_back( 0.130846 );
  correctResults.push_back( 0.131114 );
  correctResults.push_back( 0.130846 );
  correctResults.push_back( 0.130833 );
  correctResults.push_back( 0.130846 );
  correctResults.push_back( 0.131114 );
  correctResults.push_back( 0.130846 );
  correctResults.push_back( 0.140626 );
  correctResults.push_back( 0.140164 );
  correctResults.push_back( 0.142447 );
  correctResults.push_back( 0.146355 );
  correctResults.push_back( 0.142447 );
  correctResults.push_back( 0.140164 );
  correctResults.push_back( 0.140626 );
  correctResults.push_back( 0.140164 );
  correctResults.push_back( 0.142447 );
  correctResults.push_back( 0.146355 );
  correctResults.push_back( 0.142447 );
  correctResults.push_back( 0.140164 );
  correctResults.push_back( 0.155466 );
  correctResults.push_back( 0.148387 );
  correctResults.push_back( 0.156614 );
  correctResults.push_back( 0.200291 );
  correctResults.push_back( 0.156614 );
  correctResults.push_back( 0.148387 );
  correctResults.push_back( 0.155466 );
  correctResults.push_back( 0.148387 );
  correctResults.push_back( 0.156614 );
  correctResults.push_back( 0.200291 );
  correctResults.push_back( 0.156614 );
  correctResults.push_back( 0.148387 );
  correctResults.push_back( 0.198878 );
  correctResults.push_back( 0.208744 );
  correctResults.push_back( 0.198878 );
  correctResults.push_back( 0.208744 );
  correctResults.push_back( 0.198878 );
  correctResults.push_back( 0.208744 );
  correctResults.push_back( 0.198878 );
  correctResults.push_back( 0.208744 );

  calculatedResults.clear();
  for( size_t i = 0; i < size; ++i ) {
    calculatedResults.push_back( (*aRf)[ 0 ][ i ] );
  }

  EXPECT_EQ( correctResults.size(), calculatedResults.size() );
  for( size_t i = 0; i < size; ++i ) {
    EXPECT_NEAR( correctResults[i], calculatedResults[i], 1e-5 );
  }

  // Test first row for reflectance matrix
  shared_ptr< CSquareMatrix > aRb = aResults->Rho( Side::Back );

  correctResults.clear();
  correctResults.push_back( 0.128103 );
  correctResults.push_back( 0.130833 );
  correctResults.push_back( 0.130846 );
  correctResults.push_back( 0.131114 );
  correctResults.push_back( 0.130846 );
  correctResults.push_back( 0.130833 );
  correctResults.push_back( 0.130846 );
  correctResults.push_back( 0.131114 );
  correctResults.push_back( 0.130846 );
  correctResults.push_back( 0.140626 );
  correctResults.push_back( 0.140164 );
  correctResults.push_back( 0.142447 );
  correctResults.push_back( 0.146355 );
  correctResults.push_back( 0.142447 );
  correctResults.push_back( 0.140164 );
  correctResults.push_back( 0.140626 );
  correctResults.push_back( 0.140164 );
  correctResults.push_back( 0.142447 );
  correctResults.push_back( 0.146355 );
  correctResults.push_back( 0.142447 );
  correctResults.push_back( 0.140164 );
  correctResults.push_back( 0.155466 );
  correctResults.push_back( 0.148387 );
  correctResults.push_back( 0.156614 );
  correctResults.push_back( 0.200291 );
  correctResults.push_back( 0.156614 );
  correctResults.push_back( 0.148387 );
  correctResults.push_back( 0.155466 );
  correctResults.push_back( 0.148387 );
  correctResults.push_back( 0.156614 );
  correctResults.push_back( 0.200291 );
  correctResults.push_back( 0.156614 );
  correctResults.push_back( 0.148387 );
  correctResults.push_back( 0.198878 );
  correctResults.push_back( 0.208744 );
  correctResults.push_back( 0.198878 );
  correctResults.push_back( 0.208744 );
  correctResults.push_back( 0.198878 );
  correctResults.push_back( 0.208744 );
  correctResults.push_back( 0.198878 );
  correctResults.push_back( 0.208744 );

  calculatedResults.clear();
  for( size_t i = 0; i < size; ++i ) {
    calculatedResults.push_back( (*aRb)[ 0 ][ i ] );
  }

  EXPECT_EQ( correctResults.size(), calculatedResults.size() );
  for( size_t i = 0; i < size; ++i ) {
    EXPECT_NEAR( correctResults[i], calculatedResults[i], 1e-5 );
  }

}