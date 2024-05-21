#ifndef DECONSTRUCT_H
#define DECONSTRUCT_H

#include "SpectralSample.hpp"
#include <memory>
#include <WCECommon.hpp>
#include <WCESpectralAveraging.hpp>

namespace SpectralAveraging
{
    struct MonolithicInternalOpticalProperty
    {
        FenestrationCommon::CSeries ts;
        FenestrationCommon::CSeries rs;
        FenestrationCommon::CSeries taus;
    };

    struct CoatingInternalOpticalProperty
    {
        FenestrationCommon::CSeries tc;
        FenestrationCommon::CSeries rfc;
        FenestrationCommon::CSeries rbc;
    };

    MonolithicInternalOpticalProperty
      MonolithicDeconstruct(const std::shared_ptr<CSpectralSampleData> & t_SampleData);

    CoatingInternalOpticalProperty
      CoatedDeconstruct(const std::shared_ptr<CSpectralSample> & sampdat,
                        const std::shared_ptr<CSpectralSample> & subdat);

    FenestrationCommon::CSeries
      LaminatedDeconstruct(const std::shared_ptr<CSpectralSample> & sampdat,
                           const std::shared_ptr<CSpectralSample> & subdat1,
                           const std::shared_ptr<CSpectralSample> & subdat2);

    FenestrationCommon::CSeries
      LaminatedDeconstruct(const std::shared_ptr<CSpectralSample> & sampdat,
                           const std::shared_ptr<CSpectralSample> & subdat2);

    CoatingInternalOpticalProperty
      EmbeddedCoatingDeconstruct(const std::shared_ptr<CSpectralSample> & sampdat,
                                 const std::shared_ptr<CSpectralSample> & subdat1,
                                 const std::shared_ptr<CSpectralSample> & subdat2,
                                 const std::shared_ptr<CSpectralSample> lamdat);

    CoatingInternalOpticalProperty
      EmbeddedCoatingDeconstruct(const std::shared_ptr<CSpectralSample> & sampdat,
                                 const std::shared_ptr<CSpectralSample> & subdat,
                                 const std::shared_ptr<CSpectralSample> lamdat);

}   // namespace SpectralAveraging


#endif   // DECONSTRUCT_H
