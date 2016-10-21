#include <assert.h>

#include "SpecularBSDFLayer.hpp"
#include "SpecularCell.hpp"
#include "BSDFDirections.hpp"
#include "BSDFPatch.hpp"
#include "BSDFIntegrator.hpp"
#include "FenestrationCommon.hpp"
#include "BeamDirection.hpp"
#include "SquareMatrix.hpp"

using namespace std;
using namespace FenestrationCommon;

namespace SingleLayerOptics {

  CSpecularBSDFLayer::CSpecularBSDFLayer( const shared_ptr< CSpecularCell >& t_Cell, 
    const shared_ptr< const CBSDFHemisphere >& t_Hemisphere ) : CBSDFLayer( t_Cell, t_Hemisphere ) {

  }


  shared_ptr< CSpecularCell > CSpecularBSDFLayer::cellAsSpecular() const {
    shared_ptr< CSpecularCell > aCell = dynamic_pointer_cast< CSpecularCell >( m_Cell );
    assert( aCell != nullptr );
    return aCell;
  }

  void CSpecularBSDFLayer::calcDiffuseDistribution( const Side , const CBeamDirection& , const size_t ) {
    // No diffuse calculations are necessary for specular layer. 
  }

  void CSpecularBSDFLayer::calcDiffuseDistribution_wv( const Side , const CBeamDirection& , const size_t ) {
    // No diffuse calculations are necessary for specular layer.
  }

}