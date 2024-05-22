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
      CoatedDeconstruct(const std::shared_ptr<CSpectralSample> & sampdat,
                        const std::shared_ptr<CSpectralSample> & subdat);

    FenestrationCommon::CSeries
      LaminatedDeconstruct(const std::shared_ptr<CSpectralSample> & sampdat,
                           const std::shared_ptr<CSpectralSample> & subdat1,
                           const std::shared_ptr<CSpectralSample> & subdat2);

    FenestrationCommon::CSeries
      LaminatedDeconstruct(const std::shared_ptr<CSpectralSample> & sampdat,
                           const std::shared_ptr<CSpectralSample> & subdat2);

    SurfaceOpticalProperty
      EmbeddedCoatingDeconstruct(const std::shared_ptr<CSpectralSample> & sampdat,
                                 const std::shared_ptr<CSpectralSample> & subdat1,
                                 const std::shared_ptr<CSpectralSample> & subdat2,
                                 const std::shared_ptr<CSpectralSample> lamdat);

    SurfaceOpticalProperty
      EmbeddedCoatingDeconstruct(const std::shared_ptr<CSpectralSample> & sampdat,
                                 const std::shared_ptr<CSpectralSample> & subdat,
                                 const std::shared_ptr<CSpectralSample> lamdat);

}   // namespace SpectralAveraging


#endif   // DECONSTRUCT_H
