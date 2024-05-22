#include <memory>
#include <gtest/gtest.h>

#include "WCESpectralAveraging.hpp"


using namespace SpectralAveraging;
using namespace FenestrationCommon;

class TestDeconstructors : public testing::Test
{
private:
    std::shared_ptr<CSpectralSampleData> m_SubstrateData;
    std::shared_ptr<CSpectralSampleData> m_CoatedData;
    std::shared_ptr<CSpectralSampleData> m_LaminateData;
    std::shared_ptr<CSpectralSampleData> m_EmbeddedData;

protected:
    std::shared_ptr<CSpectralSampleData> getSubstrateMeasurements() const
    {
        std::shared_ptr<CSpectralSampleData> aMeasurements_sub =
          std::make_shared<CSpectralSampleData>();

        aMeasurements_sub->addRecord(300, 0.202, 0.055, 0.055);
        aMeasurements_sub->addRecord(305, 0.277, 0.059, 0.059);
        aMeasurements_sub->addRecord(310, 0.351, 0.063, 0.063);
        aMeasurements_sub->addRecord(315, 0.426, 0.066, 0.066);
        aMeasurements_sub->addRecord(320, 0.505, 0.070, 0.070);
        return aMeasurements_sub;
    }

    std::shared_ptr<CSpectralSampleData> getLaminateMeasurements() const
    {
        std::shared_ptr<CSpectralSampleData> aMeasurements_lami =
          std::make_shared<CSpectralSampleData>();

        aMeasurements_lami->addRecord(300, 0.003, 0.056, 0.056);
        aMeasurements_lami->addRecord(305, 0.003, 0.057, 0.057);
        aMeasurements_lami->addRecord(310, 0.004, 0.058, 0.058);
        aMeasurements_lami->addRecord(315, 0.004, 0.058, 0.058);
        aMeasurements_lami->addRecord(320, 0.004, 0.058, 0.058);
        return aMeasurements_lami;
    }

    std::shared_ptr<CSpectralSampleData> getCoatedMeasurements() const
    {
        std::shared_ptr<CSpectralSampleData> aMeasurements_coat =
          std::make_shared<CSpectralSampleData>();
        aMeasurements_coat->addRecord(300, 0.001, 0.276, 0.060);
        aMeasurements_coat->addRecord(305, 0.002, 0.280, 0.063);
        aMeasurements_coat->addRecord(310, 0.004, 0.292, 0.068);
        aMeasurements_coat->addRecord(315, 0.006, 0.312, 0.074);
        aMeasurements_coat->addRecord(320, 0.009, 0.338, 0.084);

        return aMeasurements_coat;
    }

    std::shared_ptr<CSpectralSampleData> getEmbeddedMeasurements() const
    {
        std::shared_ptr<CSpectralSampleData> aMeasurements_embed =
          std::make_shared<CSpectralSampleData>();

        aMeasurements_embed->addRecord(300, 0.002, 0.060, 0.056);
        aMeasurements_embed->addRecord(305, 0.002, 0.063, 0.058);
        aMeasurements_embed->addRecord(310, 0.002, 0.068, 0.059);
        aMeasurements_embed->addRecord(315, 0.002, 0.074, 0.059);
        aMeasurements_embed->addRecord(320, 0.002, 0.083, 0.058);

        return aMeasurements_embed;
    }

    void SetUp() override
    {
        m_SubstrateData = getSubstrateMeasurements();
        m_CoatedData = getCoatedMeasurements();
        m_LaminateData = getLaminateMeasurements();
        m_EmbeddedData = getEmbeddedMeasurements();
    }

public:
    std::shared_ptr<CSpectralSampleData> getSubstrateData() const
    {
        return m_SubstrateData;
    };
    std::shared_ptr<CSpectralSampleData> getCoatedData() const
    {
        return m_CoatedData;
    };
    std::shared_ptr<CSpectralSampleData> getLaminateData() const
    {
        return m_LaminateData;
    };
    std::shared_ptr<CSpectralSampleData> getEmbeddedData() const
    {
        return m_EmbeddedData;
    };
};

TEST_F(TestDeconstructors, TestMonoDeconstruct)
{
    const auto aSampleData = getSubstrateData();
    MonolithicInternalOpticalProperty prop = MonolithicDeconstruct(aSampleData);
    const auto ts = prop.surface.ts;
    const auto rs = prop.surface.rfs;
    const auto taus = prop.taus;
    // Answers from Jacob's Matlab code, high tolerance.
    EXPECT_NEAR(0.9474, ts[0].value(), 1e-4);
    EXPECT_NEAR(0.9457, ts[1].value(), 1e-4);
    EXPECT_NEAR(0.9446, ts[2].value(), 1e-4);
    EXPECT_NEAR(0.9451, ts[3].value(), 1e-4);
    EXPECT_NEAR(0.9455, ts[4].value(), 1e-4);

    EXPECT_NEAR(0.0526, rs[0].value(), 1e-4);
    EXPECT_NEAR(0.0543, rs[1].value(), 1e-4);
    EXPECT_NEAR(0.0554, rs[2].value(), 1e-4);
    EXPECT_NEAR(0.0549, rs[3].value(), 1e-4);
    EXPECT_NEAR(0.0545, rs[4].value(), 1e-4);

    EXPECT_NEAR(0.2250, taus[0].value(), 1e-4);
    EXPECT_NEAR(0.3097, taus[1].value(), 1e-4);
    EXPECT_NEAR(0.3932, taus[2].value(), 1e-4);
    EXPECT_NEAR(0.4766, taus[3].value(), 1e-4);
    EXPECT_NEAR(0.5643, taus[4].value(), 1e-4);
}

TEST_F(TestDeconstructors, TestCoatedDeconstruct)
{
    const auto aCoatedData = getCoatedData();
    const auto aSubstrateData = getSubstrateData();
    SurfaceOpticalProperty prop = CoatedDeconstruct(aCoatedData, aSubstrateData);
    const auto ts = prop.ts;
    const auto rfs = prop.rfs;
    const auto rbs = prop.rbs;
    // Answers from Jacob's Matlab code, high tolerance.
    EXPECT_NEAR(0.0046, ts[0].value(), 1e-4);
    EXPECT_NEAR(0.0067, ts[1].value(), 1e-4);
    EXPECT_NEAR(0.0106, ts[2].value(), 1e-4);
    EXPECT_NEAR(0.0131, ts[3].value(), 1e-4);
    EXPECT_NEAR(0.0166, ts[4].value(), 1e-4);

    EXPECT_NEAR(4.8517, rfs[0].value(), 1e-4);
    EXPECT_NEAR(2.5960, rfs[1].value(), 1e-4);
    EXPECT_NEAR(1.6908, rfs[2].value(), 1e-4);
    EXPECT_NEAR(1.2477, rfs[3].value(), 1e-4);
    EXPECT_NEAR(0.9789, rfs[4].value(), 1e-4);

    EXPECT_NEAR(0.0600, rbs[0].value(), 1e-4);
    EXPECT_NEAR(0.0630, rbs[1].value(), 1e-4);
    EXPECT_NEAR(0.0680, rbs[2].value(), 1e-4);
    EXPECT_NEAR(0.0740, rbs[3].value(), 1e-4);
    EXPECT_NEAR(0.0840, rbs[4].value(), 1e-4);
}

TEST_F(TestDeconstructors, TestLaminateDeconstruct)
{
    const auto aLaminateData = getLaminateData();
    const auto aSubstrateData = getSubstrateData();
    FenestrationCommon::CSeries taupvb = LaminateDeconstruct(aLaminateData, aSubstrateData);
    // Answers from Jacob's Matlab code, high tolerance.
    EXPECT_NEAR(0.0660, taupvb[0].value(), 1e-4);
    EXPECT_NEAR(0.0350, taupvb[1].value(), 1e-4);
    EXPECT_NEAR(0.0290, taupvb[2].value(), 1e-4);
    EXPECT_NEAR(0.0197, taupvb[3].value(), 1e-4);
    EXPECT_NEAR(0.0140, taupvb[4].value(), 1e-4);
}

TEST_F(TestDeconstructors, TestEmbeddedDeconstruct)
{
    const auto aEmbeddedData = getEmbeddedData();
    const auto aLaminateData = getLaminateData();
    const auto aSubstrateData = getSubstrateData();
    SurfaceOpticalProperty prop = EmbeddedCoatingDeconstruct(aEmbeddedData, aSubstrateData, aLaminateData);
    const auto ts = prop.ts;
    const auto rfs = prop.rfs;
    const auto rbs = prop.rbs;
    // Answers from Jacob's Matlab code, high tolerance.
    EXPECT_NEAR(0.6662, ts[0].value(), 1e-4);
    EXPECT_NEAR(0.6662, ts[1].value(), 1e-4);
    EXPECT_NEAR(0.4995, ts[2].value(), 1e-4);
    EXPECT_NEAR(0.4993, ts[3].value(), 1e-4);
    EXPECT_NEAR(0.4990, ts[4].value(), 1e-4);

    EXPECT_NEAR(0.1626, rfs[0].value(), 1e-4);
    EXPECT_NEAR(0.1009, rfs[1].value(), 1e-4);
    EXPECT_NEAR(0.0916, rfs[2].value(), 1e-4);
    EXPECT_NEAR(0.0942, rfs[3].value(), 1e-4);
    EXPECT_NEAR(0.1000, rfs[4].value(), 1e-4);

    EXPECT_NEAR(17.1212, rbs[0].value(), 1e-4);
    EXPECT_NEAR(34.8716, rbs[1].value(), 1e-4);
    EXPECT_NEAR(31.3679, rbs[2].value(), 1e-4);
    EXPECT_NEAR(52.4519, rbs[3].value(), 1e-4);
    EXPECT_NEAR(62.7053, rbs[4].value(), 1e-4);
}
