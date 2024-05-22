#ifndef DECONSTRUCT_H
#define DECONSTRUCT_H

#include "SpectralSample.hpp"
#include <memory>
#include <WCECommon.hpp>
#include <WCESpectralAveraging.hpp>

namespace SpectralAveraging
{

    struct SurfaceOpticalProperty
    {
        FenestrationCommon::CSeries ts;
        FenestrationCommon::CSeries rfs;
        FenestrationCommon::CSeries rbs;
    };

    struct MonolithicInternalOpticalProperty
    {
        SurfaceOpticalProperty surface;
        FenestrationCommon::CSeries taus;
    };

    MonolithicInternalOpticalProperty
      MonolithicDeconstruct(const std::shared_ptr<CSpectralSampleData> & t_SampleData);

    SurfaceOpticalProperty
      CoatedDeconstruct(const std::shared_ptr<CSpectralSampleData> & sampdat,
                        const std::shared_ptr<CSpectralSampleData> & subdat);

    FenestrationCommon::CSeries
      LaminateDeconstruct(const std::shared_ptr<CSpectralSampleData> & sampdat,
                           const std::shared_ptr<CSpectralSampleData> & subdat1,
                           const std::shared_ptr<CSpectralSampleData> & subdat2);

    FenestrationCommon::CSeries
      LaminateDeconstruct(const std::shared_ptr<CSpectralSampleData> & sampdat,
                           const std::shared_ptr<CSpectralSampleData> & subdat);

    SurfaceOpticalProperty
      EmbeddedCoatingDeconstruct(const std::shared_ptr<CSpectralSampleData> & sampdat,
                                 const std::shared_ptr<CSpectralSampleData> & subdat1,
                                 const std::shared_ptr<CSpectralSampleData> & subdat2,
                                 const std::shared_ptr<CSpectralSampleData> lamdat);

    SurfaceOpticalProperty
      EmbeddedCoatingDeconstruct(const std::shared_ptr<CSpectralSampleData> & sampdat,
                                 const std::shared_ptr<CSpectralSampleData> & subdat,
                                 const std::shared_ptr<CSpectralSampleData> lamdat);

}   // namespace SpectralAveraging


#endif   // DECONSTRUCT_H
