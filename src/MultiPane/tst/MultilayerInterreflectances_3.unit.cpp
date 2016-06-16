#include <memory>
#include <gtest/gtest.h>

#include "MultiLayerInterRef.hpp"
#include "OpticalSurface.hpp"
#include "OpticalLayer.hpp"
#include "FenestrationCommon.hpp"

using namespace std;
using namespace MultiPane;
using namespace LayerOptics;
using namespace FenestrationCommon;

// Calculation of energies that are incoming to the layers surfaces.
// Layers are added at back side.
class TestMultilayerInterreflectances_3 : public testing::Test {

private:
  // Additional layer added to the back side
  shared_ptr< CInterRef > m_Interref;

protected:
  virtual void SetUp() {

    shared_ptr< CScatteringSurface > aFront = 
      make_shared< CScatteringSurface >( 0.06, 0.04, 0.46, 0.12, 0.46, 0.52 );
    shared_ptr< CScatteringSurface > aBack =
      make_shared< CScatteringSurface >( 0.11, 0.26, 0.34, 0.19, 0.64, 0.22 );
    shared_ptr< CLayer > aLayer1 = make_shared< CLayer >( aFront, aBack );

    aFront = make_shared< CScatteringSurface >( 0.1, 0.05, 0.48, 0.26, 0.56, 0.34 );
    aBack =  make_shared< CScatteringSurface >( 0.15, 0, 0.38, 0.19, 0.49, 0.39 );
    shared_ptr< CLayer > aLayer2 = make_shared< CLayer >( aFront, aBack );

    aFront = make_shared< CScatteringSurface >( 0.08, 0.05, 0.46, 0.23, 0.46, 0.52 );
    aBack = make_shared< CScatteringSurface >( 0.13, 0.25, 0.38, 0.19, 0.64, 0.22 );
    shared_ptr< CLayer > aLayer3 = make_shared< CLayer >( aFront, aBack );
    
    m_Interref = make_shared< CInterRef >( aLayer1 );
    m_Interref->addLayer( aLayer2, Side::Back );
    m_Interref->addLayer( aLayer3, Side::Back );
  
  }

public:
  shared_ptr< CInterRef > getInt() { return m_Interref; };

};

TEST_F( TestMultilayerInterreflectances_3, TestForwardFlowFrontSide ) {
  SCOPED_TRACE( "Begin Test: Triple pane equivalent layer properties (Forward flow - Front Side)." );
  
  shared_ptr< CInterRef > eqLayer = getInt();

  Side aFlow = Side::Front;
  Side aSide = Side::Front;

  // Direct-Direct
  Scattering aScattering = Scattering::DirectDirect;
  double If1 = eqLayer->getEnergyToSurface( 1, aSide, aFlow, aScattering );
  EXPECT_NEAR( 1.0, If1, 1e-6 );

  double If2 = eqLayer->getEnergyToSurface( 2, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.060802286, If2, 1e-6 );

  double If3 = eqLayer->getEnergyToSurface( 3, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.006080229, If3, 1e-6 );

  // Diffuse-Diffuse
  aScattering = Scattering::DiffuseDiffuse;
  If1 = eqLayer->getEnergyToSurface( 1, aSide, aFlow, aScattering );
  EXPECT_NEAR( 1.0, If1, 1e-6 );

  If2 = eqLayer->getEnergyToSurface( 2, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.519291111, If2, 1e-6 );

  If3 = eqLayer->getEnergyToSurface( 3, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.36478051, If3, 1e-6 );

  // Direct-Diffuse
  aScattering = Scattering::DirectDiffuse;
  If1 = eqLayer->getEnergyToSurface( 1, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0, If1, 1e-6 );

  If2 = eqLayer->getEnergyToSurface( 2, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.526442585, If2, 1e-6 );

  If3 = eqLayer->getEnergyToSurface( 3, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.407170225, If3, 1e-6 );

}

TEST_F( TestMultilayerInterreflectances_3, TestForwardFlowBackSide ) {
  SCOPED_TRACE( "Begin Test: Triple pane equivalent layer properties (Forward flow - Back Side)." );

  shared_ptr< CInterRef > eqLayer = getInt();

  Side aFlow = Side::Front;
  Side aSide = Side::Back;

  // Direct-Direct
  Scattering aScattering = Scattering::DirectDirect;
  double If1 = eqLayer->getEnergyToSurface( 1, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.003085716, If1, 1e-6 );

  double If2 = eqLayer->getEnergyToSurface( 2, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.000304011, If2, 1e-6 );

  double If3 = eqLayer->getEnergyToSurface( 3, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0, If3, 1e-6 );

  // Diffuse-Diffuse
  aScattering = Scattering::DiffuseDiffuse;
  If1 = eqLayer->getEnergyToSurface( 1, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.269505052, If1, 1e-6 );

  If2 = eqLayer->getEnergyToSurface( 2, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.189685865, If2, 1e-6 );

  If3 = eqLayer->getEnergyToSurface( 3, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0, If3, 1e-6 );

  // Direct-Diffuse
  aScattering = Scattering::DirectDiffuse;
  If1 = eqLayer->getEnergyToSurface( 1, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.299346813, If1, 1e-6 );

  If2 = eqLayer->getEnergyToSurface( 2, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.21312697, If2, 1e-6 );

  If3 = eqLayer->getEnergyToSurface( 3, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0, If3, 1e-6 );

}

TEST_F( TestMultilayerInterreflectances_3, TestBackwardFlowFrontSide ) {
  SCOPED_TRACE( "Begin Test: Triple pane equivalent layer properties (Backward flow - Front Side)." );

  shared_ptr< CInterRef > eqLayer = getInt();

  Side aFlow = Side::Back;
  Side aSide = Side::Front;

  // Direct-Direct
  Scattering aScattering = Scattering::DirectDirect;
  double If1 = eqLayer->getEnergyToSurface( 1, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0, If1, 1e-6 );

  double If2 = eqLayer->getEnergyToSurface( 2, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.005137793, If2, 1e-6 );

  double If3 = eqLayer->getEnergyToSurface( 3, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.000513779, If3, 1e-6 );

  // Diffuse-Diffuse
  aScattering = Scattering::DiffuseDiffuse;
  If1 = eqLayer->getEnergyToSurface( 1, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0, If1, 1e-6 );

  If2 = eqLayer->getEnergyToSurface( 2, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.097697737, If2, 1e-6 );

  If3 = eqLayer->getEnergyToSurface( 3, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.381724451, If3, 1e-6 );

  // Direct-Diffuse
  aScattering = Scattering::DirectDiffuse;
  If1 = eqLayer->getEnergyToSurface( 1, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0, If1, 1e-6 );

  If2 = eqLayer->getEnergyToSurface( 2, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.07702437, If2, 1e-6 );

  If3 = eqLayer->getEnergyToSurface( 3, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.274147962, If3, 1e-6 );

}

TEST_F( TestMultilayerInterreflectances_3, TestBackwardFlowBackSide ) {
  SCOPED_TRACE( "Begin Test: Triple pane equivalent layer properties (Backward flow - Back Side)." );

  shared_ptr< CInterRef > eqLayer = getInt();

  Side aFlow = Side::Back;
  Side aSide = Side::Back;

  // Direct-Direct
  Scattering aScattering = Scattering::DirectDirect;
  double If1 = eqLayer->getEnergyToSurface( 1, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.019760743, If1, 1e-6 );

  double If2 = eqLayer->getEnergyToSurface( 2, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.130025689, If2, 1e-6 );

  double If3 = eqLayer->getEnergyToSurface( 3, aSide, aFlow, aScattering );
  EXPECT_NEAR( 1, If3, 1e-6 );

  // Diffuse-Diffuse
  aScattering = Scattering::DiffuseDiffuse;
  If1 = eqLayer->getEnergyToSurface( 1, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.444080621, If1, 1e-6 );

  If2 = eqLayer->getEnergyToSurface( 2, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.838496715, If2, 1e-6 );

  If3 = eqLayer->getEnergyToSurface( 3, aSide, aFlow, aScattering );
  EXPECT_NEAR( 1, If3, 1e-6 );

  // Direct-Diffuse
  aScattering = Scattering::DirectDiffuse;
  If1 = eqLayer->getEnergyToSurface( 1, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.333044677, If1, 1e-6 );

  If2 = eqLayer->getEnergyToSurface( 2, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0.522675109, If2, 1e-6 );

  If3 = eqLayer->getEnergyToSurface( 3, aSide, aFlow, aScattering );
  EXPECT_NEAR( 0, If3, 1e-6 );

}