#include <memory>
#include <stdexcept>
#include <gtest/gtest.h>

#include "WCETarcog.hpp"
#include "WCECommon.hpp"

// Example of double clear window with inital guess for solution
class TestDoubleClear_EN673 : public testing::Test
{
private:
    std::unique_ptr<Tarcog::EN673::IGU> m_IGU;

protected:
    void SetUp() override
    {
        /////////////////////////////////////////////////////////
        /// Outdoor
        /////////////////////////////////////////////////////////
        auto airTemperature = 273.15;   // Kelvins
        auto filmCoefficient = 23;      // [W/m2K]
        const auto outdoor = Tarcog::EN673::Environment(airTemperature, filmCoefficient);

        /////////////////////////////////////////////////////////
        /// Indoor
        /////////////////////////////////////////////////////////

        airTemperature = 293.15;   // Kelvins
        filmCoefficient = 8;      // [W/m2K]
        const auto indoor = Tarcog::EN673::Environment(airTemperature, filmCoefficient);

        //////////////////////////////////////////////////////////
        /// First layer
        //////////////////////////////////////////////////////////
        const auto thickness = 0.003;   // [m]
        const auto conductivity = 1.0;     // [W/m2K]
        const auto emissFront = 0.84;
        const auto emissBack = 0.84;

        const auto layer1 = Tarcog::EN673::Glass(conductivity, thickness, emissFront, emissBack);

        /////////////////////////////////////////////////////////
        /// IGU
        /////////////////////////////////////////////////////////
        m_IGU =
          std::unique_ptr<Tarcog::EN673::IGU>(new Tarcog::EN673::IGU(indoor, outdoor, layer1));

        /////////////////////////////////////////////////////////
        /// gap and layer
        /////////////////////////////////////////////////////////
        Gases::CGas gas;
        /// 100% Air
        gas.addGasItem(1.0, Gases::GasDef::Air);

        const auto gapThickness = 0.0127;   // [mm]
        const auto pressure = 101325;       // [Pa]
        const auto gap = Tarcog::EN673::Gap(gapThickness, pressure, gas);

        m_IGU->addGap(gap);

        const auto layer2 = Tarcog::EN673::Glass(conductivity, thickness, emissFront, emissBack);
        m_IGU->addGlass(layer2);
    }

public:
    Tarcog::EN673::IGU * GetIGU() const
    {
        return m_IGU.get();
    };
};

TEST_F(TestDoubleClear_EN673, Test1)
{
    SCOPED_TRACE("Begin Test: Uvalue");

    auto igu = GetIGU();

    auto Uvalue = igu->Uvalue();

    EXPECT_NEAR(2.82511, Uvalue, 1e-5);
}
