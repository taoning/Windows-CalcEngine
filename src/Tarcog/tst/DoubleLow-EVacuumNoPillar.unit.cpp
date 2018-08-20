#include <memory>
#include <stdexcept>
#include <gtest/gtest.h>

#include "WCETarcog.hpp"
#include "WCECommon.hpp"

using namespace Tarcog;
using namespace FenestrationCommon;

class DoubleLowEVacuumNoPillar : public testing::Test
{
private:
    std::shared_ptr<CSingleSystem> m_TarcogSystem;

protected:
    void SetUp() override
    {
        /////////////////////////////////////////////////////////
        /// Outdoor
        /////////////////////////////////////////////////////////
        auto airTemperature = 255.15;   // Kelvins
        auto pressure = 101325.0;       // Pascals
        auto airSpeed = 5.5;            // meters per second
        auto tSky = 255.15;   // Kelvins
        auto solarRadiation = 0.0;

        auto Outdoor = Environments::outdoor( airTemperature,
											  pressure,
											  airSpeed,
											  solarRadiation,
											  tSky,
											  SkyModel::AllSpecified);
        ASSERT_TRUE(Outdoor != nullptr);
        Outdoor->setHCoeffModel(BoundaryConditionsCoeffModel::CalculateH);

        /////////////////////////////////////////////////////////
        /// Indoor
        /////////////////////////////////////////////////////////

        auto roomTemperature = 294.15;

        auto Indoor = Environments::indoor( roomTemperature, pressure );
        ASSERT_TRUE(Indoor != nullptr);

        /////////////////////////////////////////////////////////
        /// IGU
        /////////////////////////////////////////////////////////
        auto solidLayerThickness = 0.004;   // [m]
        auto solidLayerConductance = 1.0;
        auto TransmittanceIR = 0.0;
        auto emissivityFrontIR = 0.84;
        auto emissivityBackIR = 0.036749500781;

        auto aSolidLayer1 = Layers::solid(solidLayerThickness,
                                          solidLayerConductance,
                                          emissivityFrontIR,
                                          TransmittanceIR,
                                          emissivityBackIR,
                                          TransmittanceIR);

        solidLayerThickness = 0.003962399904;
        emissivityBackIR = 0.84;

        auto aSolidLayer2 = Layers::solid(solidLayerThickness,
                                          solidLayerConductance,
                                          emissivityFrontIR,
                                          TransmittanceIR,
                                          emissivityBackIR,
                                          TransmittanceIR);

        auto gapThickness = 0.0001;
        auto gapPressure = 0.1333;
        auto m_GapLayer = Layers::gap(gapThickness, gapPressure);
        ASSERT_TRUE(m_GapLayer != nullptr);

        auto windowWidth = 1.0;   //[m]
        auto windowHeight = 1.0;
        CIGU aIGU(windowWidth, windowHeight);
        aIGU.addLayer(aSolidLayer1);
        aIGU.addLayer(m_GapLayer);
        aIGU.addLayer(aSolidLayer2);

        /////////////////////////////////////////////////////////
        /// System
        /////////////////////////////////////////////////////////
        m_TarcogSystem = std::make_shared<CSingleSystem>(aIGU, Indoor, Outdoor);
        ASSERT_TRUE(m_TarcogSystem != nullptr);

        m_TarcogSystem->solve();
    }

public:
    std::shared_ptr<CSingleSystem> GetSystem() const
    {
        return m_TarcogSystem;
    };
};

TEST_F(DoubleLowEVacuumNoPillar, Test1)
{
    SCOPED_TRACE("Begin Test: Double Low-E - vacuum with no pillar support");

    auto aSystem = GetSystem();

    ASSERT_TRUE(aSystem != nullptr);

    auto Temperature = *aSystem->getTemperatures();
    std::vector<double> correctTemperature = {255.501938, 255.543003, 292.514948, 292.555627};
    ASSERT_EQ(correctTemperature.size(), Temperature.size());

    for(auto i = 0u; i < correctTemperature.size(); ++i)
    {
        EXPECT_NEAR(correctTemperature[i], Temperature[i], 1e-5);
    }

    auto Radiosity = *aSystem->getRadiosities();
    std::vector<double> correctRadiosity = {241.409657, 407.569595, 413.894817, 416.791085};
    ASSERT_EQ(correctRadiosity.size(), Radiosity.size());

    for(auto i = 0u; i < correctRadiosity.size(); ++i)
    {
        EXPECT_NEAR(correctRadiosity[i], Radiosity[i], 1e-5);
    }

    auto numOfIter = aSystem->getNumberOfIterations();
    EXPECT_EQ(30u, numOfIter);
}
