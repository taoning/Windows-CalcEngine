#include <memory>
#include <gtest/gtest.h>

#include "WCETarcog.hpp"

class TestGapLayerSealedForcedVentilation : public testing::Test
{
private:
    std::shared_ptr<Tarcog::ISO15099::CSingleSystem> m_TarcogSystem;

protected:
    void SetUp() override
    {
        /////////////////////////////////////////////////////////
        // Outdoor
        /////////////////////////////////////////////////////////
        auto airTemperature = 255.15;   // Kelvins
        auto airSpeed = 5.5;            // meters per second
        auto tSky = 255.15;             // Kelvins
        auto solarRadiation = 0.0;

        auto Outdoor = Tarcog::ISO15099::Environments::outdoor(
          airTemperature, airSpeed, solarRadiation, tSky, Tarcog::ISO15099::SkyModel::AllSpecified);
        ASSERT_TRUE(Outdoor != nullptr);
        Outdoor->setHCoeffModel(Tarcog::ISO15099::BoundaryConditionsCoeffModel::CalculateH);

        /////////////////////////////////////////////////////////
        /// Indoor
        /////////////////////////////////////////////////////////

        auto roomTemperature = 295.15;

        auto Indoor = Tarcog::ISO15099::Environments::indoor(roomTemperature);
        ASSERT_TRUE(Indoor != nullptr);

        // IGU
        auto solidLayerThickness = 0.005715;   // [m]
        auto solidLayerConductance = 1.0;

        auto solidLayer1 =
          Tarcog::ISO15099::Layers::solid(solidLayerThickness, solidLayerConductance);
        ASSERT_TRUE(solidLayer1 != nullptr);

        auto solidLayer2 =
          Tarcog::ISO15099::Layers::solid(solidLayerThickness, solidLayerConductance);
        ASSERT_TRUE(solidLayer2 != nullptr);

        auto gapThickness = 0.0127;
        auto gapAirSpeed = 0.5;
        Tarcog::ISO15099::ForcedVentilation forcedVentilation = {gapAirSpeed};
        auto gap = Tarcog::ISO15099::Layers::forcedVentilationGap(gapThickness, forcedVentilation);
        ASSERT_TRUE(gap != nullptr);

        double windowWidth = 1;
        double windowHeight = 1;
        Tarcog::ISO15099::CIGU aIGU(windowWidth, windowHeight);
        aIGU.addLayers({solidLayer1, gap, solidLayer2});

        /////////////////////////////////////////////////////////
        /// System
        /////////////////////////////////////////////////////////
        m_TarcogSystem = std::make_shared<Tarcog::ISO15099::CSingleSystem>(aIGU, Indoor, Outdoor);
        ASSERT_TRUE(m_TarcogSystem != nullptr);

        m_TarcogSystem->solve();
    }

public:
    std::shared_ptr<Tarcog::ISO15099::CIGUSolidLayer> GetSolidLayer1() const
    {
        auto solidLayer = m_TarcogSystem->getSolidLayers()[0];
        assert(std::dynamic_pointer_cast<Tarcog::ISO15099::CIGUSolidLayer>(solidLayer) != nullptr);
        return std::dynamic_pointer_cast<Tarcog::ISO15099::CIGUSolidLayer>(solidLayer);
    };

    std::shared_ptr<Tarcog::ISO15099::CIGUVentilatedGapLayer> GetGap() const
    {
        auto gap = m_TarcogSystem->getGapLayers()[0];
        assert(std::dynamic_pointer_cast<Tarcog::ISO15099::CIGUVentilatedGapLayer>(gap) != nullptr);
        return std::dynamic_pointer_cast<Tarcog::ISO15099::CIGUVentilatedGapLayer>(gap);
    };

    std::shared_ptr<Tarcog::ISO15099::CIGUSolidLayer> GetSolidLayer2() const
    {
        auto solidLayer = m_TarcogSystem->getSolidLayers()[1];
        assert(std::dynamic_pointer_cast<Tarcog::ISO15099::CIGUSolidLayer>(solidLayer) != nullptr);
        return std::dynamic_pointer_cast<Tarcog::ISO15099::CIGUSolidLayer>(solidLayer);
    };
};

TEST_F(TestGapLayerSealedForcedVentilation, GainEnergy)
{
    SCOPED_TRACE("Begin Test: Test Sealed Forced Ventilated Gap Layer - Gain Energy");

    auto aLayer = GetGap();

    // Airflow iterations are set to 1e-4 and it cannot exceed that precision

    ASSERT_TRUE(aLayer != nullptr);
    auto gainEnergy = aLayer->getGainFlow();
    EXPECT_NEAR(-57.193291279467125, gainEnergy, 1e-4);
}

TEST_F(TestGapLayerSealedForcedVentilation, Solid1Temperatures)
{
    SCOPED_TRACE("Begin Test: Test Sealed Forced Ventilated Gap Layer - Solid 1 Temperatures");

    auto aLayer = GetSolidLayer1();

    // Airflow iterations are set to 1e-4 and it cannot exceed that precision

    ASSERT_TRUE(aLayer != nullptr);
    auto frontTemperature = aLayer->getTemperature(FenestrationCommon::Side::Front);
    auto backTemperature = aLayer->getTemperature(FenestrationCommon::Side::Back);
    EXPECT_NEAR(257.88610911376719, frontTemperature, 1e-4);
    EXPECT_NEAR(258.34294987245784, backTemperature, 1e-4);
}

TEST_F(TestGapLayerSealedForcedVentilation, GapTemperatures)
{
    SCOPED_TRACE("Begin Test: Test Sealed Forced Ventilated Gap Layer - Gap Temperatures");

    auto aLayer = GetGap();

    // Airflow iterations are set to 1e-4 and it cannot exceed that precision

    ASSERT_TRUE(aLayer != nullptr);
    auto frontTemperature = aLayer->getTemperature(FenestrationCommon::Side::Front);
    auto backTemperature = aLayer->getTemperature(FenestrationCommon::Side::Back);
    auto layerTemperature = aLayer->layerTemperature();
    auto averageTemperature = aLayer->averageTemperature();
    EXPECT_NEAR(258.34294987245784, frontTemperature, 1e-4);
    EXPECT_NEAR(276.01222466583204, backTemperature, 1e-4);
    EXPECT_NEAR(262.42249308270669, layerTemperature, 1e-4);
    EXPECT_NEAR(267.17758726914496, averageTemperature, 1e-4);
}

TEST_F(TestGapLayerSealedForcedVentilation, Solid2Temperatures)
{
    SCOPED_TRACE("Begin Test: Test Sealed Forced Ventilated Gap Layer - Solid 2 Temperatures");

    auto aLayer = GetSolidLayer2();

    // Airflow iterations are set to 1e-4 and it cannot exceed that precision

    ASSERT_TRUE(aLayer != nullptr);
    auto frontTemperature = aLayer->getTemperature(FenestrationCommon::Side::Front);
    auto backTemperature = aLayer->getTemperature(FenestrationCommon::Side::Back);
    EXPECT_NEAR(276.01222466583204, frontTemperature, 1e-4);
    EXPECT_NEAR(276.79592508330529, backTemperature, 1e-4);
}

TEST_F(TestGapLayerSealedForcedVentilation, AirflowReferencePoint)
{
    SCOPED_TRACE("Begin Test: Test Sealed Forced Ventilated Gap Layer - Airflow Reference Point");

    auto aLayer = GetGap();

    // Airflow iterations are set to 1e-4 and it cannot exceed that precision

    ASSERT_TRUE(aLayer != nullptr);
    auto airflowReferencePoint = aLayer->getAirflowReferencePoint(0.5);
    EXPECT_NEAR(6911.4449859804226, airflowReferencePoint, 1e-4);
}
