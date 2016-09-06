#include <memory>
#include <gtest/gtest.h>

#include "PerforatedCell.hpp"
#include "PerforatedCellDescription.hpp"
#include "MaterialDescription.hpp"
#include "FenestrationCommon.hpp"
#include "BeamDirection.hpp"

using namespace std;
using namespace LayerOptics;
using namespace FenestrationCommon;

class TestRectangularPerforatedCell : public testing::Test {

private:
  shared_ptr< CRectangularCellDescription > m_DescriptionCell;
  shared_ptr< CPerforatedCell > m_PerforatedCell;

protected:
  virtual void SetUp()
  {
    // create material
    double Tmat = 0.1;
    double Rfmat = 0.7;
    double Rbmat = 0.8;
    double minLambda = 0.3;
    double maxLambda = 2.5;
    shared_ptr< CMaterialSingleBand > aMaterial = 
      make_shared< CMaterialSingleBand >( Tmat, Tmat, Rfmat, Rbmat, minLambda, maxLambda );

    // make cell geometry
    double x = 10; // mm
    double y = 10; // mm
    double thickness = 1; // mm
    double xHole = 5; // mm
    double yHole = 5; // mm
    m_DescriptionCell = make_shared< CRectangularCellDescription >( x, y, thickness, xHole, yHole );

    m_PerforatedCell = make_shared< CPerforatedCell >( aMaterial, m_DescriptionCell );
  }

public:
  shared_ptr< CPerforatedCell > GetCell() { return m_PerforatedCell; };
  shared_ptr< CRectangularCellDescription > GetDescription() { return m_DescriptionCell; };

};

TEST_F( TestRectangularPerforatedCell, TestRectangular1 )
{
  SCOPED_TRACE( "Begin Test: Rectangular perforated cell (Theta = 0, Phi = 0)." );
  
  shared_ptr< CPerforatedCell > aCell = GetCell();
  shared_ptr< CCellDescription > aCellDescription = GetDescription();

  double Theta = 0; // deg
  double Phi = 0; // deg
  Side aFrontSide = Side::Front;
  Side aBackSide = Side::Back;

  shared_ptr< CBeamDirection > aDirection = make_shared< CBeamDirection >( Theta, Phi );
  
  double Tdir_dir = aCellDescription->T_dir_dir( aFrontSide, *aDirection );
  EXPECT_NEAR( 0.25, Tdir_dir, 1e-6 );
  
  double Tdir_dif = aCell->T_dir_dif( aFrontSide, aDirection );
  EXPECT_NEAR( 0.075, Tdir_dif, 1e-6 );

  double Rfdir_dif = aCell->R_dir_dif( aFrontSide, aDirection );
  EXPECT_NEAR( 0.525, Rfdir_dif, 1e-6 );

  double Rbdir_dif = aCell->R_dir_dif( aBackSide, aDirection );
  EXPECT_NEAR( 0.6, Rbdir_dif, 1e-6 );

}

TEST_F( TestRectangularPerforatedCell, TestRectangular2 )
{
  SCOPED_TRACE( "Begin Test: Rectangular perforated cell (Theta = 45, Phi = 0)." );
  
  shared_ptr< CPerforatedCell > aCell = GetCell();
  shared_ptr< CCellDescription > aCellDescription = GetDescription();

  double Theta = 45; // deg
  double Phi = 0; // deg
  Side aFrontSide = Side::Front;
  Side aBackSide = Side::Back;

  shared_ptr< CBeamDirection > aDirection = make_shared< CBeamDirection >( Theta, Phi );
  
  double Tdir_dir = aCellDescription->T_dir_dir( aFrontSide, *aDirection );
  EXPECT_NEAR( 0.2, Tdir_dir, 1e-6 );
  
  double Tdir_dif = aCell->T_dir_dif( aFrontSide, aDirection );
  EXPECT_NEAR( 0.08, Tdir_dif, 1e-6 );

  double Rfdir_dif = aCell->R_dir_dif( aFrontSide, aDirection );
  EXPECT_NEAR( 0.56, Rfdir_dif, 1e-6 );

  double Rbdir_dif = aCell->R_dir_dif( aBackSide, aDirection );
  EXPECT_NEAR( 0.64, Rbdir_dif, 1e-6 );

}

TEST_F( TestRectangularPerforatedCell, TestRectangular3 )
{
  SCOPED_TRACE( "Begin Test: Rectangular perforated cell (Theta = 45, Phi = 45)." );
  
  shared_ptr< CPerforatedCell > aCell = GetCell();
  shared_ptr< CCellDescription > aCellDescription = GetDescription();

  double Theta = 45; // deg
  double Phi = 45; // deg
  Side aFrontSide = Side::Front;
  Side aBackSide = Side::Back;

  shared_ptr< CBeamDirection > aDirection = make_shared< CBeamDirection >( Theta, Phi );
  
  double Tdir_dir = aCellDescription->T_dir_dir( aFrontSide, *aDirection );
  EXPECT_NEAR( 0.184289322, Tdir_dir, 1e-6 );
  
  double Tdir_dif = aCell->T_dir_dif( aFrontSide, aDirection );
  EXPECT_NEAR( 0.081571068, Tdir_dif, 1e-6 );

  double Rfdir_dif = aCell->R_dir_dif( aFrontSide, aDirection );
  EXPECT_NEAR( 0.570997475, Rfdir_dif, 1e-6 );

  double Rbdir_dif = aCell->R_dir_dif( aBackSide, aDirection );
  EXPECT_NEAR( 0.652568542, Rbdir_dif, 1e-6 );

}