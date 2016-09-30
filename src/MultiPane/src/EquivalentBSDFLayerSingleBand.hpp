#ifndef EQUIVALENTBSDFLAYERSINGLEBAND_H
#define EQUIVALENTBSDFLAYERSINGLEBAND_H

#include <memory>
#include <vector>

namespace FenestrationCommon {

  class CSquareMatrix;
  enum class Side;
  enum class PropertySimple;

}

namespace LayerOptics {

  class CBSDFResults;

}

namespace MultiPane {

  // Class to handle interreflectance calculations
  class CInterReflectance {
  public:
    CInterReflectance( const FenestrationCommon::CSquareMatrix& t_Lambda,
      const FenestrationCommon::CSquareMatrix& t_Rb,
      const FenestrationCommon::CSquareMatrix& t_Rf );

    std::shared_ptr< FenestrationCommon::CSquareMatrix > value() const;

  private:
    std::shared_ptr< FenestrationCommon::CSquareMatrix > m_InterRefl;

  };

  // Class to calculate equivalent BSDF transmittance and reflectances. This will be used by
  // multilayer routines to calculate properties for any number of layers.
  class CBSDFDoubleLayer {
  public:
    CBSDFDoubleLayer( const LayerOptics::CBSDFResults& t_FrontLayer, 
      const LayerOptics::CBSDFResults& t_BackLayer );

    std::shared_ptr< LayerOptics::CBSDFResults > value();

  private:
    std::shared_ptr< FenestrationCommon::CSquareMatrix > equivalentT( 
      const FenestrationCommon::CSquareMatrix& t_Tf2, 
      const FenestrationCommon::CSquareMatrix& t_InterRefl,
      const FenestrationCommon::CSquareMatrix& t_Lambda,
      const FenestrationCommon::CSquareMatrix& t_Tf1 );

    std::shared_ptr< FenestrationCommon::CSquareMatrix > equivalentR( 
      const FenestrationCommon::CSquareMatrix& t_Rf1,
      const FenestrationCommon::CSquareMatrix& t_Tf1,
      const FenestrationCommon::CSquareMatrix& t_Tb1,
      const FenestrationCommon::CSquareMatrix& t_Rf2,
      const FenestrationCommon::CSquareMatrix& t_InterRefl,
      const FenestrationCommon::CSquareMatrix& t_Lambda );

    std::shared_ptr< LayerOptics::CBSDFResults > m_Results;

    std::shared_ptr< FenestrationCommon::CSquareMatrix > m_Tf;
    std::shared_ptr< FenestrationCommon::CSquareMatrix > m_Tb;
    std::shared_ptr< FenestrationCommon::CSquareMatrix > m_Rf;
    std::shared_ptr< FenestrationCommon::CSquareMatrix > m_Rb;
  };

  // Class for equivalent BSDF layer for single material properties (or single wavelength)
  class CEquivalentBSDFLayerSingleBand {
  public:
    explicit CEquivalentBSDFLayerSingleBand( const std::shared_ptr< LayerOptics::CBSDFResults >& t_Layer );
    void addLayer( const std::shared_ptr< LayerOptics::CBSDFResults >& t_Layer );

    std::shared_ptr< FenestrationCommon::CSquareMatrix > Tau( const FenestrationCommon::Side t_Side );
    std::shared_ptr< FenestrationCommon::CSquareMatrix > Rho( const FenestrationCommon::Side t_Side );

    std::shared_ptr< FenestrationCommon::CSquareMatrix > getProperty( const FenestrationCommon::Side t_Side,
      const FenestrationCommon::PropertySimple t_Property );

    std::shared_ptr< std::vector< double > > getLayerAbsorptances( const size_t Index, 
      FenestrationCommon::Side t_Side );

    size_t getNumberOfLayers() const;
  
  private:
    void calcEquivalentProperties();

    std::shared_ptr< std::vector< double > > absTerm1( const std::vector< double >& t_Alpha,
      const FenestrationCommon::CSquareMatrix& t_InterRefl,
      const FenestrationCommon::CSquareMatrix& t_T );

    std::shared_ptr< std::vector< double > > absTerm2( const std::vector< double >& t_Alpha,
      const FenestrationCommon::CSquareMatrix& t_InterRefl,
      const FenestrationCommon::CSquareMatrix& t_R,
      const FenestrationCommon::CSquareMatrix& t_T );

    std::shared_ptr< LayerOptics::CBSDFResults > m_EquivalentLayer;
    std::vector< std::shared_ptr< LayerOptics::CBSDFResults > > m_Layers;
  
    // Forward and backward layers are used for calculation of equivalent absorptances
    std::vector< std::shared_ptr< LayerOptics::CBSDFResults > > m_Forward;
    std::vector< std::shared_ptr< LayerOptics::CBSDFResults > > m_Backward;
  
    std::shared_ptr< std::vector< std::shared_ptr< std::vector< double > > > > m_Af;
    std::shared_ptr< std::vector< std::shared_ptr< std::vector< double > > > > m_Ab;

    bool m_PropertiesCalculated;
    
    std::shared_ptr< const FenestrationCommon::CSquareMatrix > m_Lambda;
  };

}

#endif