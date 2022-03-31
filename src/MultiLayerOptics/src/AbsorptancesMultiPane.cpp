#include "AbsorptancesMultiPane.hpp"
#include "WCECommon.hpp"

using namespace FenestrationCommon;

namespace MultiLayerOptics
{
    CAbsorptancesMultiPane::CAbsorptancesMultiPane(const CSeries & t_T,
                                                   const CSeries & t_Rf,
                                                   const CSeries & t_Rb) :
        m_R({{Side::Front, std::vector<CSeries>()}, {Side::Back, std::vector<CSeries>()}}),
        m_Abs({{Side::Front, std::vector<CSeries>()}, {Side::Back, std::vector<CSeries>()}}),
        m_rCoeffs({{Side::Front, std::vector<CSeries>()}, {Side::Back, std::vector<CSeries>()}}),
        m_tCoeffs({{Side::Front, std::vector<CSeries>()}, {Side::Back, std::vector<CSeries>()}}),
        Iplus({{Side::Front, std::vector<CSeries>()}, {Side::Back, std::vector<CSeries>()}}),
        Iminus({{Side::Front, std::vector<CSeries>()}, {Side::Back, std::vector<CSeries>()}}),
        m_StateCalculated(false)
    {
        m_T.push_back(t_T);
        m_R.at(Side::Front).push_back(t_Rf);
        m_R.at(Side::Back).push_back(t_Rb);
    }

    void CAbsorptancesMultiPane::addLayer(const CSeries & t_T,
                                          const CSeries & t_Rf,
                                          const CSeries & t_Rb)
    {
        m_T.push_back(t_T);
        m_R.at(Side::Front).push_back(t_Rf);
        m_R.at(Side::Back).push_back(t_Rb);
        m_StateCalculated = false;
    }

    CSeries CAbsorptancesMultiPane::Abs(const size_t Index)
    {
        calculateState();
        return m_Abs.at(Side::Front)[Index];
    }

    size_t CAbsorptancesMultiPane::numOfLayers()
    {
        calculateState();
        return m_Abs.size();
    }

    CSeries CAbsorptancesMultiPane::iplus(size_t Index)
    {
        calculateState();
        return Iplus.at(Side::Front)[Index];
    }

    CSeries CAbsorptancesMultiPane::iminus(size_t Index)
    {
        calculateState();
        return Iminus.at(Side::Front)[Index];
    }

    void CAbsorptancesMultiPane::calculateRTCoefficients()
    {
        const size_t size{m_T.size()};

        for(const auto side: EnumSide())
        {
            const auto oppositeSide{FenestrationCommon::oppositeSide(side)};

            // Calculate r and t coefficients
            CSeries r;
            CSeries t;
            const auto wv{m_T[size - 1].getXArray()};
            r.setConstantValues(wv, 0);
            t.setConstantValues(wv, 0);
            m_rCoeffs.at(side).clear();
            m_tCoeffs.at(side).clear();

            // TODO: Coefficients need to be in other direction for backward calculations
            for(int i = static_cast<int>(size) - 1; i >= 0; --i)
            {
                t = tCoeffs(m_T[i], m_R.at(oppositeSide)[i], r);
                r = rCoeffs(m_T[i], m_R.at(side)[i], m_R.at(oppositeSide)[i], r);

                m_rCoeffs.at(side).insert(m_rCoeffs.at(side).begin(), r);
                m_tCoeffs.at(side).insert(m_tCoeffs.at(side).begin(), t);
            }
        }
    }

    void CAbsorptancesMultiPane::calculateNormalizedRadiances()
    {
        // Calculate normalized radiances
        const auto wv{m_T[0].getXArray()};

        for(const auto side : EnumSide())
        {
            const size_t size{m_rCoeffs.at(side).size()};

            Iplus.at(side).clear();
            Iminus.at(side).clear();

            CSeries Im;
            CSeries Ip;
            Im.setConstantValues(wv, 1);
            Iminus.at(side).push_back(Im);

            for(size_t i = 0; i < size; ++i)
            {
                Ip = m_rCoeffs.at(side)[i] * Im;
                Im = m_tCoeffs.at(side)[i] * Im;
                if(i != 0u)
                {
                    Iplus.at(side).push_back(Ip);
                }
                Iminus.at(side).push_back(Im);
            }
            Ip.setConstantValues(wv, 0);
            Iplus.at(side).push_back(Ip);
        }
    }

    void CAbsorptancesMultiPane::calculateAbsorptances()
    {
        for(const auto side : EnumSide())
        {
            const auto oppositeSide{FenestrationCommon::oppositeSide(side)};
            const size_t size{Iminus.at(side).size()};
            m_Abs.at(side).clear();
            for(size_t i = 0; i < size - 1; ++i)
            {
                const auto Afront{1 - m_T[i] - m_R.at(side)[i]};
                const auto Aback{1 - m_T[i] - m_R.at(oppositeSide)[i]};
                const auto Ifront = Iminus.at(side)[i] * Afront;
                const auto Iback = Iplus.at(side)[i] * Aback;
                m_Abs.at(side).emplace_back(Ifront + Iback);
            }
        }
    }

    void CAbsorptancesMultiPane::calculateState()
    {
        if(!m_StateCalculated)
        {
            calculateRTCoefficients();

            calculateNormalizedRadiances();

            calculateAbsorptances();
        }
    }

    CSeries CAbsorptancesMultiPane::rCoeffs(const CSeries & t_T,
                                            const CSeries & t_Rf,
                                            const CSeries & t_Rb,
                                            const CSeries & t_RCoeffs)
    {
        CSeries rCoeffs;
        size_t size = t_T.size();

        for(size_t i = 0; i < size; ++i)
        {
            double wl = t_T[i].x();
            double rValue = t_Rf[i].value()
                            + t_T[i].value() * t_T[i].value() * t_RCoeffs[i].value()
                                / (1 - t_Rb[i].value() * t_RCoeffs[i].value());
            rCoeffs.addProperty(wl, rValue);
        }

        return rCoeffs;
    }

    CSeries CAbsorptancesMultiPane::tCoeffs(const CSeries & t_T,
                                            const CSeries & t_Rb,
                                            const CSeries & t_RCoeffs)
    {
        CSeries tCoeffs;
        size_t size = t_T.size();

        for(size_t i = 0; i < size; ++i)
        {
            double wl = t_T[i].x();
            double tValue = t_T[i].value() / (1 - t_Rb[i].value() * t_RCoeffs[i].value());
            tCoeffs.addProperty(wl, tValue);
        }

        return tCoeffs;
    }

}   // namespace MultiLayerOptics
