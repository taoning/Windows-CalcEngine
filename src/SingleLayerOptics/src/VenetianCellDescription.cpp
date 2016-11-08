#include <assert.h>
#include <memory>

#include "VenetianCellDescription.hpp"
#include "VenetianSlat.hpp"
#include "ViewSegment2D.hpp"
#include "Geometry2DBeam.hpp"
#include "Geometry2D.hpp"
#include "Segment2D.hpp"
#include "BeamDirection.hpp"
#include "Point2D.hpp"
#include "SquareMatrix.hpp"

using namespace std;
using namespace Viewer;
using namespace FenestrationCommon;

namespace SingleLayerOptics {

  CVenetianCellDescription::CVenetianCellDescription( const double t_SlatWidth, const double t_SlatSpacing, 
    const double t_SlatTiltAngle, const double t_CurvatureRadius, const size_t t_NumOfSlatSegments ) :
    m_Top( make_shared< CVenetianSlat >( t_SlatWidth, t_SlatSpacing, t_SlatTiltAngle, t_CurvatureRadius, 
      t_NumOfSlatSegments, SegmentsDirection::Positive ) ),
    m_Bottom( make_shared< CVenetianSlat >( t_SlatWidth, 0, t_SlatTiltAngle, t_CurvatureRadius, 
      t_NumOfSlatSegments, SegmentsDirection::Negative ) ) {

    shared_ptr< CViewSegment2D > exteriorSegment = 
      make_shared< CViewSegment2D >( m_Bottom->geometry()->lastPoint(), m_Top->geometry()->firstPoint() );

    shared_ptr< CViewSegment2D > interiorSegment = 
      make_shared< CViewSegment2D >( m_Top->geometry()->lastPoint(), m_Bottom->geometry()->firstPoint() );
    
    m_Geometry = make_shared< CGeometry2D >();
    m_Geometry->appendSegment( exteriorSegment );
    m_Geometry->appendGeometry2D( m_Top->geometry() );
    m_Geometry->appendSegment( interiorSegment );
    m_Geometry->appendGeometry2D( m_Bottom->geometry() );

    m_BeamGeometry = make_shared< CGeometry2DBeam >();
    m_BeamGeometry->appendGeometry2D( m_Top->geometry() );
    m_BeamGeometry->appendGeometry2D( m_Bottom->geometry() );
  }

  size_t CVenetianCellDescription::numberOfSegments() const {
    assert( m_Top != nullptr );
    assert( m_Bottom != nullptr );
    // Two additional segments are for interior and exterior openess
    return 2 + m_Top->numberOfSegments() + m_Bottom->numberOfSegments();
  }

  double CVenetianCellDescription::segmentLength( const size_t Index ) const {
    shared_ptr< vector< shared_ptr< CViewSegment2D > > > aSegments =  m_Geometry->segments();
    if( Index > aSegments->size() ) {
      throw runtime_error("Incorrect index for venetian segment.");
    }
    shared_ptr< CViewSegment2D > aSegment = ( *aSegments )[ Index ];
    return aSegment->length(); 
  }

  shared_ptr< CVenetianCellDescription > CVenetianCellDescription::makeBackwardCell() const {
    double slatWidth = m_Top->slatWidth();
    double slatSpacing = m_Top->slatSpacing();
    double slatTiltAngle = -m_Top->slatTiltAngle();
    double curvatureRadius = m_Top->curvatureRadius();
    size_t m_NumOfSlatSegments = m_Top->numberOfSegments();

    shared_ptr< CVenetianCellDescription > aBackwardCell = 
      make_shared< CVenetianCellDescription >( slatWidth, slatSpacing, slatTiltAngle, 
        curvatureRadius, m_NumOfSlatSegments );

    return aBackwardCell;
  }

  shared_ptr< CSquareMatrix > CVenetianCellDescription::viewFactors() {
    return m_Geometry->viewFactors();
  }

  shared_ptr< vector < BeamViewFactor > > CVenetianCellDescription::beamViewFactors( 
    const double t_ProfileAngle, const Side t_Side ) {
    assert( m_BeamGeometry != nullptr );
    return m_BeamGeometry->beamViewFactors( -t_ProfileAngle, t_Side );
  }

  double CVenetianCellDescription::T_dir_dir( const Side t_Side, const CBeamDirection& t_Direction ) {
    assert( m_BeamGeometry != nullptr );
    double aProfileAngle = t_Direction.profileAngle();
    return m_BeamGeometry->directToDirect( -aProfileAngle, t_Side );
  }

  double CVenetianCellDescription::R_dir_dir( const Side, const CBeamDirection& ) {
    return 0;
  }

}