#ifndef GASITEM_H
#define GASITEM_H

#include <memory>
#include <string>

namespace Gases
{
    class CGasData;
    struct GasProperties;

    double const DefaultPressure = 101325;
    double const DefaultTemperature = 273.15;

    class CGasItem
    {
    public:
        CGasItem();
        CGasItem(CGasItem const & t_GasItem);
        CGasItem(double aFraction, CGasData const & t_GasData);
        CGasItem & operator=(CGasItem const & t_GasItem);
        void setTemperature(double t_Temperature);
        void setPressure(double t_Pressure);
        [[nodiscard]] double fraction() const;
        [[nodiscard]] std::shared_ptr<GasProperties> getFractionalGasProperties() const;
        [[nodiscard]] std::shared_ptr<GasProperties> getGasProperties() const;
        bool operator==(CGasItem const & rhs) const;
        bool operator!=(CGasItem const & rhs) const;

        [[nodiscard]] std::string name() const;

    private:
        void fillStandardPressureProperites() const;
        void flllVacuumPressureProperties() const;
        void initialize();
        void resetCalculatedProperties() const;
        double m_Temperature;   // unit in Kelvins
        double m_Pressure;      // unit in Pa
        double m_Fraction;      // value between 0 and 1
        std::shared_ptr<GasProperties> m_GasProperties;
        std::shared_ptr<GasProperties> m_FractionalGasProperties;
        std::unique_ptr<CGasData> m_GasData;
    };

}   // namespace Gases

#endif
