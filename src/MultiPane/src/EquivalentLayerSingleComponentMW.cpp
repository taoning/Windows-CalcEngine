#include <assert.h>

#include "EquivalentLayerSingleComponentMW.hpp"
#include "Series.hpp"
#include "FenestrationCommon.hpp"
#include "EquivalentLayerSingleComponent.hpp"

using namespace FenestrationCommon;
using namespace std;

namespace MultiPane {

  ///////////////////////////////////////////////////////////////////////////
  ///   CSurfaceSeries
  ///////////////////////////////////////////////////////////////////////////

  CSurfaceSeries::CSurfaceSeries( shared_ptr< CSeries > t_T, shared_ptr< CSeries > t_R ) {
    m_Properties[ Property::T ] = t_T;
    m_Properties[ Property::R ] = t_R;
    size_t size = t_T->size();
    shared_ptr< CSeries > aAbs = make_shared< CSeries >();
    for( size_t i = 0; i < size; ++i ) {
      double wl = ( *t_T )[ i ]->x();
      double value = 1 - ( *t_T )[ i ]->value() - ( *t_R )[ i ]->value();
      if( value > 1 || value < 0 ) {
        throw runtime_error("Absorptance value for provided series is out of range.");
      }
      aAbs->addProperty( wl, value );
    }
    m_Properties[ Property::Abs ] = aAbs;
  }

  shared_ptr< CSeries > CSurfaceSeries::getProperties( const Property t_Property ) const {
    return m_Properties.at( t_Property );
  }

  ///////////////////////////////////////////////////////////////////////////
  ///   CLayerSeries
  ///////////////////////////////////////////////////////////////////////////

  CLayerSeries::CLayerSeries( shared_ptr< CSeries > t_Tf, shared_ptr< CSeries > t_Rf, shared_ptr< CSeries > t_Tb, 
    shared_ptr< CSeries > t_Rb ) {
    m_Surfaces[ Side::Front ] = make_shared< CSurfaceSeries >( t_Tf, t_Rf );
    m_Surfaces[ Side::Back ] = make_shared< CSurfaceSeries >( t_Tb, t_Rb );
  }

  shared_ptr< CSeries > CLayerSeries::getProperties( const Side t_Side, const Property t_Property ) const {
    return m_Surfaces.at( t_Side )->getProperties( t_Property );
  }

  ///////////////////////////////////////////////////////////////////////////
  ///   CEquivalentLayerSingleComponentMW
  ///////////////////////////////////////////////////////////////////////////

  CEquivalentLayerSingleComponentMW::CEquivalentLayerSingleComponentMW( shared_ptr< CSeries > t_Tf, shared_ptr< CSeries > t_Tb, 
    shared_ptr< CSeries > t_Rf, shared_ptr< CSeries > t_Rb  ) {
    m_Layer = make_shared< CLayerSeries >( t_Tf, t_Rf, t_Tb, t_Rb );

    size_t size = t_Tf->size();
    for(size_t i = 0; i < size; ++i) {
      shared_ptr< CEquivalentLayerSingleComponent > aLayer = make_shared< CEquivalentLayerSingleComponent >( ( *t_Tf )[i]->value(),
        ( *t_Rf )[i]->value(), ( *t_Tb )[i]->value(), ( *t_Rb )[i]->value() );
      m_EqLayerBySeries.push_back( aLayer );
    }
  }

  void CEquivalentLayerSingleComponentMW::addLayer( shared_ptr< CSeries > t_Tf, shared_ptr< CSeries > t_Tb, 
    shared_ptr< CSeries > t_Rf, shared_ptr< CSeries > t_Rb ) {

    size_t size = t_Tf->size();

    for( size_t i = 0; i < size; ++i ) {
      shared_ptr< CEquivalentLayerSingleComponent > aLayer = m_EqLayerBySeries[ i ];
      aLayer->addLayer( ( *t_Tf )[i]->value(), ( *t_Rf )[i]->value(), ( *t_Tb )[i]->value(), ( *t_Rb )[i]->value() );
    }

    shared_ptr< CSeries > tTotf = make_shared< CSeries >();
    shared_ptr< CSeries > tTotb = make_shared< CSeries >();
    shared_ptr< CSeries > tRfTot = make_shared< CSeries >();
    shared_ptr< CSeries > tRbTot = make_shared< CSeries >();

    for( size_t i = 0; i < size; ++i ) {
      double wl = ( *t_Tf )[ i ]->x();
      
      double Tf = m_EqLayerBySeries[ i ]->getProperty( Property::T, Side::Front );
      tTotf->addProperty( wl, Tf );

      double Rf = m_EqLayerBySeries[ i ]->getProperty( Property::R, Side::Front );
      tRfTot->addProperty( wl, Rf );

      double Tb = m_EqLayerBySeries[ i ]->getProperty( Property::T, Side::Back );
      tTotb->addProperty( wl, Tb );

      double Rb = m_EqLayerBySeries[ i ]->getProperty( Property::R, Side::Back );
      tRbTot->addProperty( wl, Rb );
    }
    
    m_Layer = make_shared< CLayerSeries >( tTotf, tRfTot, tTotb, tRbTot );

  }

  shared_ptr< CSeries > CEquivalentLayerSingleComponentMW::getProperties( const Property t_Property, const Side t_Side ) const {
    return m_Layer->getProperties( t_Side, t_Property );
  }

}