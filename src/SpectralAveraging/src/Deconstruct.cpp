#include <memory>
#include <vector>

#include "SpectralSample.hpp"
#include "MeasuredSampleData.hpp"
#include "Deconstruct.hpp"
#include "WCECommon.hpp"

namespace SpectralAveraging
{

    MonolithicInternalOpticalProperty
      MonolithicDeconstruct(const std::shared_ptr<CSpectralSampleData> & sampdat)
    {
        MonolithicInternalOpticalProperty property;
        FenestrationCommon::CSeries Ts = sampdat->properties(
          FenestrationCommon::Property::T, FenestrationCommon::Side::Front);
        FenestrationCommon::CSeries Rfs = sampdat->properties(
          FenestrationCommon::Property::R, FenestrationCommon::Side::Front);
        FenestrationCommon::CSeries Rbs = sampdat->properties(
          FenestrationCommon::Property::R, FenestrationCommon::Side::Back);

        FenestrationCommon::CSeries Rs = (Rfs + Rbs) * 0.5;

        FenestrationCommon::CSeries beta = Ts * Ts + Rs * 2 - Rs * Rs + 1;

        FenestrationCommon::CSeries gamma = 4 - Rs * 2;

        auto rs = (beta - sqrt(beta * beta - gamma * 2 * Rs)) / gamma;
        auto ts = 1 - rs;
        auto taus = (Rs - rs) / rs / Ts;

        // check if taus has nagative value
        return {{ts, rs,rs}, taus};
    }

    SurfaceOpticalProperty CoatedDeconstruct(const std::shared_ptr<CSpectralSampleData> & sampdat,
                                                    const std::shared_ptr<CSpectralSampleData> & subdat)
    {
        MonolithicInternalOpticalProperty sprop = MonolithicDeconstruct(subdat);
        FenestrationCommon::CSeries Tc = sampdat->properties(
          FenestrationCommon::Property::T, FenestrationCommon::Side::Front);
        FenestrationCommon::CSeries Rfc = sampdat->properties(
          FenestrationCommon::Property::R, FenestrationCommon::Side::Front);
        FenestrationCommon::CSeries Rbc = sampdat->properties(
          FenestrationCommon::Property::R, FenestrationCommon::Side::Back);

        FenestrationCommon::CSeries tau2 = sprop.taus * sprop.taus;

        FenestrationCommon::CSeries rfc = (Rfc - sprop.surface.rfs) / tau2 / (sprop.surface.rfs * (Rfc - 2.) + 1.);

        FenestrationCommon::CSeries tc = Tc / (1 - sprop.surface.rfs * rfc * tau2) / sprop.surface.ts / sprop.taus;

        FenestrationCommon::CSeries rbc =
          Rbc - tc * tc * sprop.surface.rfs * tau2 / (1 - sprop.surface.rfs * rfc * tau2);

        // check if taus has nagative value
        return {tc, rfc, rbc};
    }

    FenestrationCommon::CSeries LaminateDeconstruct(const std::shared_ptr<CSpectralSampleData> & sampdat,
                                                    const std::shared_ptr<CSpectralSampleData> & subdat1,
                                                    const std::shared_ptr<CSpectralSampleData> & subdat2)
    {
        MonolithicInternalOpticalProperty sprop1 = MonolithicDeconstruct(subdat1);
        MonolithicInternalOpticalProperty sprop2 = MonolithicDeconstruct(subdat2);
        FenestrationCommon::CSeries Tl = sampdat->properties(
          FenestrationCommon::Property::T, FenestrationCommon::Side::Front);

        FenestrationCommon::CSeries denoms =
          Tl * sprop1.surface.rfs * sprop2.surface.rfs * sprop1.taus * sprop2.taus * 2;

        FenestrationCommon::CSeries noms = sprop1.surface.ts * sprop2.surface.ts * -1.
                                           + sqrt(sprop1.surface.ts * sprop1.surface.ts * sprop2.surface.ts * sprop2.surface.ts
                                                  + Tl * Tl * sprop1.surface.rfs * sprop2.surface.rfs * 4);

        return noms / denoms;
    }

    // For symmetrical laminates
    FenestrationCommon::CSeries LaminateDeconstruct(const std::shared_ptr<CSpectralSampleData> & sampdat,
                                                    const std::shared_ptr<CSpectralSampleData> & subdat)
    {
        return LaminateDeconstruct(sampdat, subdat, subdat);
    }


    SurfaceOpticalProperty
      EmbeddedCoatingDeconstruct(const std::shared_ptr<CSpectralSampleData> & sampdat,
                                 const std::shared_ptr<CSpectralSampleData> & subdat1,
                                 const std::shared_ptr<CSpectralSampleData> & subdat2,
                                 const std::shared_ptr<CSpectralSampleData> lamdat)
    {
        MonolithicInternalOpticalProperty sprop1 = MonolithicDeconstruct(subdat1);
        MonolithicInternalOpticalProperty sprop2 = MonolithicDeconstruct(subdat2);
        FenestrationCommon::CSeries taupvb = LaminateDeconstruct(lamdat, subdat1, subdat2);

        FenestrationCommon::CSeries Tl = sampdat->properties(
          FenestrationCommon::Property::T, FenestrationCommon::Side::Front);
        FenestrationCommon::CSeries Rfl = sampdat->properties(
          FenestrationCommon::Property::R, FenestrationCommon::Side::Front);
        FenestrationCommon::CSeries Rbl = sampdat->properties(
          FenestrationCommon::Property::R, FenestrationCommon::Side::Back);

        FenestrationCommon::CSeries Rfprim =
          (Rfl - sprop1.surface.rfs) / (sprop1.taus * sprop1.taus)
          / (sprop1.surface.ts * sprop1.surface.ts + sprop1.surface.rfs * (Rfl - sprop1.surface.rfs));
        FenestrationCommon::CSeries Tprim =
          Tl * (1 - sprop1.taus * sprop1.taus * Rfprim * sprop1.surface.rfs) / sprop1.surface.ts / sprop1.taus;
        FenestrationCommon::CSeries Rbprim =
          (Rbl - Tprim * Tprim * sprop1.surface.rfs * sprop1.taus * sprop1.taus)
          / (1 - sprop1.taus * sprop1.taus * Rfprim * sprop1.surface.rfs);

        FenestrationCommon::CSeries rgcdenom =
          (taupvb * taupvb * sprop2.taus * sprop2.taus)
          * ((Rbprim - sprop2.surface.rfs) * sprop2.surface.rfs + sprop2.surface.ts * sprop2.surface.ts);
        FenestrationCommon::CSeries rbc = (Rbprim - sprop2.surface.rfs) / rgcdenom;

        FenestrationCommon::CSeries tcdenom = sprop2.surface.ts * taupvb * sprop2.taus;
        FenestrationCommon::CSeries tcnom =
          Tprim * (1 - rbc * sprop2.surface.rfs * taupvb * taupvb * sprop2.taus * sprop2.taus);
        FenestrationCommon::CSeries tc = tcnom / tcdenom;

        FenestrationCommon::CSeries rfcdenom =
          1 - sprop2.surface.rfs * rbc * taupvb * taupvb * sprop2.taus * sprop2.taus;
        FenestrationCommon::CSeries rfcnom =
          Rfprim - (tc * tc * sprop2.surface.rfs * taupvb * taupvb * sprop2.taus * sprop2.taus);
        FenestrationCommon::CSeries rfc = rfcnom / rfcdenom;

        return {tc, rfc, rbc};
    }

    // Symmetrical laminates
    SurfaceOpticalProperty
      EmbeddedCoatingDeconstruct(const std::shared_ptr<CSpectralSampleData> & sampdat,
                                 const std::shared_ptr<CSpectralSampleData> & subdat,
                                 const std::shared_ptr<CSpectralSampleData> lamdat)
    {
        return EmbeddedCoatingDeconstruct(sampdat, subdat, subdat, lamdat);
    }


}   // namespace SpectralAveraging
