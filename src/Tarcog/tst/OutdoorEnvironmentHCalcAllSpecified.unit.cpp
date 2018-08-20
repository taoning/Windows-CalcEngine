#include <memory>
#include <stdexcept>
#include <gtest/gtest.h>

#include "WCETarcog.hpp"

using namespace Tarcog;


class TestOutdoorEnvironmentHCalcAllSpecified : public testing::Test
{
private:
    std::shared_ptr<CEnvironment> Outdoor;
    std::shared_ptr<CSingleSystem> m_TarcogSystem;

protected:
    void SetUp() override
    {
        /////////////////////////////////////////////////////////
        /// Outdoor
        /////////////////////////////////////////////////////////
        auto airTemperature = 300.0;   // Kelvins
        auto pressure = 101325.0;      // Pascals
        auto airSpeed = 5.5;           // meters per second
        auto tSky = 270.0;   // Kelvins
        auto solarRadiation = 0.0;

        Outdoor = Environments::outdoor(airTemperature,
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

        auto Indoor = Environments::indoor(roomTemperature, pressure);
        ASSERT_TRUE(Indoor != nullptr);

        /////////////////////////////////////////////////////////
        /// IGU
        /////////////////////////////////////////////////////////
        auto solidLayerThickness = 0.003048;   // [m]
        auto solidLayerConductance = 100.0;

        auto aSolidLayer = Layers::solid(solidLayerThickness, solidLayerConductance);
        ASSERT_TRUE(aSolidLayer != nullptr);

        auto windowWidth = 1.0;
        auto windowHeight = 1.0;
        CIGU aIGU(windowWidth, windowHeight);
        aIGU.addLayer(aSolidLayer);

        /////////////////////////////////////////////////////////
        /// System
        /////////////////////////////////////////////////////////
        m_TarcogSystem = std::make_shared<CSingleSystem>(aIGU, Indoor, Outdoor);
        m_TarcogSystem->solve();
        ASSERT_TRUE(m_TarcogSystem != nullptr);
    }

public:
    std::shared_ptr<CEnvironment> GetOutdoors() const
    {
        return Outdoor;
    };
};

TEST_F(TestOutdoorEnvironmentHCalcAllSpecified, CalculateH_AllSpecified)
{
    SCOPED_TRACE("Begin Test: Outdoors -> H model = Calculate; Sky Model = All Specified");

    auto aOutdoor = GetOutdoors();
    ASSERT_TRUE(aOutdoor != nullptr);

    auto radiosity = aOutdoor->getEnvironmentIR();
    EXPECT_NEAR(380.278401885, radiosity, 1e-6);

    auto hc = aOutdoor->getHc();
    EXPECT_NEAR(26, hc, 1e-6);

    auto outIR = aOutdoor->getRadiationFlow();
    EXPECT_NEAR(52.1067777, outIR, 1e-6);

    auto outConvection = aOutdoor->getConvectionConductionFlow();
    EXPECT_NEAR(-72.9257342, outConvection, 1e-6);

    auto totalHeatFlow = aOutdoor->getHeatFlow();
    EXPECT_NEAR(-20.81895645, totalHeatFlow, 1e-6);
}
