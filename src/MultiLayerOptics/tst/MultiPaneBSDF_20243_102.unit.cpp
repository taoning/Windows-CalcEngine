#include <memory>
#include <gtest/gtest.h>

#include "WCESpectralAveraging.hpp"
#include "WCEMultiLayerOptics.hpp"
#include "WCESingleLayerOptics.hpp"
#include "WCECommon.hpp"


using namespace SingleLayerOptics;
using namespace FenestrationCommon;
using namespace SpectralAveraging;
using namespace MultiLayerOptics;

// Example on how to create multilayer BSDF from specular layers only

class MultiPaneBSDF_20243_103 : public testing::Test
{
private:
    std::shared_ptr<CMultiPaneBSDF> m_Layer;

    std::shared_ptr<CSeries> loadSolarRadiationFile()
    {
        std::shared_ptr<CSeries> aSolarRadiation = std::make_shared<CSeries>();

        // Full ASTM E891-87 Table 1 (Solar radiation)
        aSolarRadiation->addProperty(0.3000, 0.0);
        aSolarRadiation->addProperty(0.3050, 3.4);
        aSolarRadiation->addProperty(0.3100, 15.6);
        aSolarRadiation->addProperty(0.3150, 41.1);
        aSolarRadiation->addProperty(0.3200, 71.2);
        aSolarRadiation->addProperty(0.3250, 100.2);
        aSolarRadiation->addProperty(0.3300, 152.4);
        aSolarRadiation->addProperty(0.3350, 155.6);
        aSolarRadiation->addProperty(0.3400, 179.4);
        aSolarRadiation->addProperty(0.3450, 186.7);
        aSolarRadiation->addProperty(0.3500, 212.0);
        aSolarRadiation->addProperty(0.3600, 240.5);
        aSolarRadiation->addProperty(0.3700, 324.0);
        aSolarRadiation->addProperty(0.3800, 362.4);
        aSolarRadiation->addProperty(0.3900, 381.7);
        aSolarRadiation->addProperty(0.4000, 556.0);
        aSolarRadiation->addProperty(0.4100, 656.3);
        aSolarRadiation->addProperty(0.4200, 690.8);
        aSolarRadiation->addProperty(0.4300, 641.9);
        aSolarRadiation->addProperty(0.4400, 798.5);
        aSolarRadiation->addProperty(0.4500, 956.6);
        aSolarRadiation->addProperty(0.4600, 990.0);
        aSolarRadiation->addProperty(0.4700, 998.0);
        aSolarRadiation->addProperty(0.4800, 1046.1);
        aSolarRadiation->addProperty(0.4900, 1005.1);
        aSolarRadiation->addProperty(0.5000, 1026.7);
        aSolarRadiation->addProperty(0.5100, 1066.7);
        aSolarRadiation->addProperty(0.5200, 1011.5);
        aSolarRadiation->addProperty(0.5300, 1084.9);
        aSolarRadiation->addProperty(0.5400, 1082.4);
        aSolarRadiation->addProperty(0.5500, 1102.2);
        aSolarRadiation->addProperty(0.5700, 1087.4);
        aSolarRadiation->addProperty(0.5900, 1024.3);
        aSolarRadiation->addProperty(0.6100, 1088.8);
        aSolarRadiation->addProperty(0.6300, 1062.1);
        aSolarRadiation->addProperty(0.6500, 1061.7);
        aSolarRadiation->addProperty(0.6700, 1046.2);
        aSolarRadiation->addProperty(0.6900, 859.2);
        aSolarRadiation->addProperty(0.7100, 1002.4);
        aSolarRadiation->addProperty(0.7180, 816.9);
        aSolarRadiation->addProperty(0.7244, 842.8);
        aSolarRadiation->addProperty(0.7400, 971.0);
        aSolarRadiation->addProperty(0.7525, 956.3);
        aSolarRadiation->addProperty(0.7575, 942.2);
        aSolarRadiation->addProperty(0.7625, 524.8);
        aSolarRadiation->addProperty(0.7675, 830.7);
        aSolarRadiation->addProperty(0.7800, 908.9);
        aSolarRadiation->addProperty(0.8000, 873.4);
        aSolarRadiation->addProperty(0.8160, 712.0);
        aSolarRadiation->addProperty(0.8237, 660.2);
        aSolarRadiation->addProperty(0.8315, 765.5);
        aSolarRadiation->addProperty(0.8400, 799.8);
        aSolarRadiation->addProperty(0.8600, 815.2);
        aSolarRadiation->addProperty(0.8800, 778.3);
        aSolarRadiation->addProperty(0.9050, 630.4);
        aSolarRadiation->addProperty(0.9150, 565.2);
        aSolarRadiation->addProperty(0.9250, 586.4);
        aSolarRadiation->addProperty(0.9300, 348.1);
        aSolarRadiation->addProperty(0.9370, 224.2);
        aSolarRadiation->addProperty(0.9480, 271.4);
        aSolarRadiation->addProperty(0.9650, 451.2);
        aSolarRadiation->addProperty(0.9800, 549.7);
        aSolarRadiation->addProperty(0.9935, 630.1);
        aSolarRadiation->addProperty(1.0400, 582.9);
        aSolarRadiation->addProperty(1.0700, 539.7);
        aSolarRadiation->addProperty(1.1000, 366.2);
        aSolarRadiation->addProperty(1.1200, 98.1);
        aSolarRadiation->addProperty(1.1300, 169.5);
        aSolarRadiation->addProperty(1.1370, 118.7);
        aSolarRadiation->addProperty(1.1610, 301.9);
        aSolarRadiation->addProperty(1.1800, 406.8);
        aSolarRadiation->addProperty(1.2000, 375.2);
        aSolarRadiation->addProperty(1.2350, 423.6);
        aSolarRadiation->addProperty(1.2900, 365.7);
        aSolarRadiation->addProperty(1.3200, 223.4);
        aSolarRadiation->addProperty(1.3500, 30.1);
        aSolarRadiation->addProperty(1.3950, 1.4);
        aSolarRadiation->addProperty(1.4425, 51.6);
        aSolarRadiation->addProperty(1.4625, 97.0);
        aSolarRadiation->addProperty(1.4770, 97.3);
        aSolarRadiation->addProperty(1.4970, 167.1);
        aSolarRadiation->addProperty(1.5200, 239.3);
        aSolarRadiation->addProperty(1.5390, 248.8);
        aSolarRadiation->addProperty(1.5580, 249.3);
        aSolarRadiation->addProperty(1.5780, 222.3);
        aSolarRadiation->addProperty(1.5920, 227.3);
        aSolarRadiation->addProperty(1.6100, 210.5);
        aSolarRadiation->addProperty(1.6300, 224.7);
        aSolarRadiation->addProperty(1.6460, 215.9);
        aSolarRadiation->addProperty(1.6780, 202.8);
        aSolarRadiation->addProperty(1.7400, 158.2);
        aSolarRadiation->addProperty(1.8000, 28.6);
        aSolarRadiation->addProperty(1.8600, 1.8);
        aSolarRadiation->addProperty(1.9200, 1.1);
        aSolarRadiation->addProperty(1.9600, 19.7);
        aSolarRadiation->addProperty(1.9850, 84.9);
        aSolarRadiation->addProperty(2.0050, 25.0);
        aSolarRadiation->addProperty(2.0350, 92.5);
        aSolarRadiation->addProperty(2.0650, 56.3);
        aSolarRadiation->addProperty(2.1000, 82.7);
        aSolarRadiation->addProperty(2.1480, 76.2);
        aSolarRadiation->addProperty(2.1980, 66.4);
        aSolarRadiation->addProperty(2.2700, 65.0);
        aSolarRadiation->addProperty(2.3600, 57.6);
        aSolarRadiation->addProperty(2.4500, 19.8);
        aSolarRadiation->addProperty(2.4940, 17.0);
        aSolarRadiation->addProperty(2.5370, 3.0);
        aSolarRadiation->addProperty(2.9410, 4.0);
        aSolarRadiation->addProperty(2.9730, 7.0);
        aSolarRadiation->addProperty(3.0050, 6.0);
        aSolarRadiation->addProperty(3.0560, 3.0);
        aSolarRadiation->addProperty(3.1320, 5.0);
        aSolarRadiation->addProperty(3.1560, 18.0);
        aSolarRadiation->addProperty(3.2040, 1.2);
        aSolarRadiation->addProperty(3.2450, 3.0);
        aSolarRadiation->addProperty(3.3170, 12.0);
        aSolarRadiation->addProperty(3.3440, 3.0);
        aSolarRadiation->addProperty(3.4500, 12.2);
        aSolarRadiation->addProperty(3.5730, 11.0);
        aSolarRadiation->addProperty(3.7650, 9.0);
        aSolarRadiation->addProperty(4.0450, 6.9);

        return aSolarRadiation;
    }

    std::shared_ptr<CSpectralSampleData> loadSampleData_NFRC_102()
    {
        std::shared_ptr<CSpectralSampleData> aMeasurements_102 =
          std::make_shared<CSpectralSampleData>();

        aMeasurements_102->addRecord(0.300, 0.0020, 0.0470, 0.0480);
        aMeasurements_102->addRecord(0.305, 0.0030, 0.0470, 0.0480);
        aMeasurements_102->addRecord(0.310, 0.0090, 0.0470, 0.0480);
        aMeasurements_102->addRecord(0.315, 0.0350, 0.0470, 0.0480);
        aMeasurements_102->addRecord(0.320, 0.1000, 0.0470, 0.0480);
        aMeasurements_102->addRecord(0.325, 0.2180, 0.0490, 0.0500);
        aMeasurements_102->addRecord(0.330, 0.3560, 0.0530, 0.0540);
        aMeasurements_102->addRecord(0.335, 0.4980, 0.0600, 0.0610);
        aMeasurements_102->addRecord(0.340, 0.6160, 0.0670, 0.0670);
        aMeasurements_102->addRecord(0.345, 0.7090, 0.0730, 0.0740);
        aMeasurements_102->addRecord(0.350, 0.7740, 0.0780, 0.0790);
        aMeasurements_102->addRecord(0.355, 0.8180, 0.0820, 0.0820);
        aMeasurements_102->addRecord(0.360, 0.8470, 0.0840, 0.0840);
        aMeasurements_102->addRecord(0.365, 0.8630, 0.0850, 0.0850);
        aMeasurements_102->addRecord(0.370, 0.8690, 0.0850, 0.0860);
        aMeasurements_102->addRecord(0.375, 0.8610, 0.0850, 0.0850);
        aMeasurements_102->addRecord(0.380, 0.8560, 0.0840, 0.0840);
        aMeasurements_102->addRecord(0.385, 0.8660, 0.0850, 0.0850);
        aMeasurements_102->addRecord(0.390, 0.8810, 0.0860, 0.0860);
        aMeasurements_102->addRecord(0.395, 0.8890, 0.0860, 0.0860);
        aMeasurements_102->addRecord(0.400, 0.8930, 0.0860, 0.0860);
        aMeasurements_102->addRecord(0.410, 0.8930, 0.0860, 0.0860);
        aMeasurements_102->addRecord(0.420, 0.8920, 0.0860, 0.0860);
        aMeasurements_102->addRecord(0.430, 0.8920, 0.0850, 0.0850);
        aMeasurements_102->addRecord(0.440, 0.8920, 0.0850, 0.0850);
        aMeasurements_102->addRecord(0.450, 0.8960, 0.0850, 0.0850);
        aMeasurements_102->addRecord(0.460, 0.9000, 0.0850, 0.0850);
        aMeasurements_102->addRecord(0.470, 0.9020, 0.0840, 0.0840);
        aMeasurements_102->addRecord(0.480, 0.9030, 0.0840, 0.0840);
        aMeasurements_102->addRecord(0.490, 0.9040, 0.0850, 0.0850);
        aMeasurements_102->addRecord(0.500, 0.9050, 0.0840, 0.0840);
        aMeasurements_102->addRecord(0.510, 0.9050, 0.0840, 0.0840);
        aMeasurements_102->addRecord(0.520, 0.9050, 0.0840, 0.0840);
        aMeasurements_102->addRecord(0.530, 0.9040, 0.0840, 0.0840);
        aMeasurements_102->addRecord(0.540, 0.9040, 0.0830, 0.0830);
        aMeasurements_102->addRecord(0.550, 0.9030, 0.0830, 0.0830);
        aMeasurements_102->addRecord(0.560, 0.9020, 0.0830, 0.0830);
        aMeasurements_102->addRecord(0.570, 0.9000, 0.0820, 0.0820);
        aMeasurements_102->addRecord(0.580, 0.8980, 0.0820, 0.0820);
        aMeasurements_102->addRecord(0.590, 0.8960, 0.0810, 0.0810);
        aMeasurements_102->addRecord(0.600, 0.8930, 0.0810, 0.0810);
        aMeasurements_102->addRecord(0.610, 0.8900, 0.0810, 0.0810);
        aMeasurements_102->addRecord(0.620, 0.8860, 0.0800, 0.0800);
        aMeasurements_102->addRecord(0.630, 0.8830, 0.0800, 0.0800);
        aMeasurements_102->addRecord(0.640, 0.8790, 0.0790, 0.0790);
        aMeasurements_102->addRecord(0.650, 0.8750, 0.0790, 0.0790);
        aMeasurements_102->addRecord(0.660, 0.8720, 0.0790, 0.0790);
        aMeasurements_102->addRecord(0.670, 0.8680, 0.0780, 0.0780);
        aMeasurements_102->addRecord(0.680, 0.8630, 0.0780, 0.0780);
        aMeasurements_102->addRecord(0.690, 0.8590, 0.0770, 0.0770);
        aMeasurements_102->addRecord(0.700, 0.8540, 0.0760, 0.0770);
        aMeasurements_102->addRecord(0.710, 0.8500, 0.0760, 0.0760);
        aMeasurements_102->addRecord(0.720, 0.8450, 0.0750, 0.0760);
        aMeasurements_102->addRecord(0.730, 0.8400, 0.0750, 0.0750);
        aMeasurements_102->addRecord(0.740, 0.8350, 0.0750, 0.0750);
        aMeasurements_102->addRecord(0.750, 0.8310, 0.0740, 0.0740);
        aMeasurements_102->addRecord(0.760, 0.8260, 0.0740, 0.0740);
        aMeasurements_102->addRecord(0.770, 0.8210, 0.0740, 0.0740);
        aMeasurements_102->addRecord(0.780, 0.8160, 0.0730, 0.0730);
        aMeasurements_102->addRecord(0.790, 0.8120, 0.0730, 0.0730);
        aMeasurements_102->addRecord(0.800, 0.8080, 0.0720, 0.0720);
        aMeasurements_102->addRecord(0.810, 0.8030, 0.0720, 0.0720);
        aMeasurements_102->addRecord(0.820, 0.8000, 0.0720, 0.0720);
        aMeasurements_102->addRecord(0.830, 0.7960, 0.0710, 0.0710);
        aMeasurements_102->addRecord(0.840, 0.7930, 0.0700, 0.0710);
        aMeasurements_102->addRecord(0.850, 0.7880, 0.0700, 0.0710);
        aMeasurements_102->addRecord(0.860, 0.7860, 0.0700, 0.0700);
        aMeasurements_102->addRecord(0.870, 0.7820, 0.0740, 0.0740);
        aMeasurements_102->addRecord(0.880, 0.7800, 0.0720, 0.0720);
        aMeasurements_102->addRecord(0.890, 0.7770, 0.0730, 0.0740);
        aMeasurements_102->addRecord(0.900, 0.7760, 0.0720, 0.0720);
        aMeasurements_102->addRecord(0.910, 0.7730, 0.0720, 0.0720);
        aMeasurements_102->addRecord(0.920, 0.7710, 0.0710, 0.0710);
        aMeasurements_102->addRecord(0.930, 0.7700, 0.0700, 0.0700);
        aMeasurements_102->addRecord(0.940, 0.7680, 0.0690, 0.0690);
        aMeasurements_102->addRecord(0.950, 0.7660, 0.0680, 0.0680);
        aMeasurements_102->addRecord(0.960, 0.7660, 0.0670, 0.0680);
        aMeasurements_102->addRecord(0.970, 0.7640, 0.0680, 0.0680);
        aMeasurements_102->addRecord(0.980, 0.7630, 0.0680, 0.0680);
        aMeasurements_102->addRecord(0.990, 0.7620, 0.0670, 0.0670);
        aMeasurements_102->addRecord(1.000, 0.7620, 0.0660, 0.0670);
        aMeasurements_102->addRecord(1.050, 0.7600, 0.0660, 0.0660);
        aMeasurements_102->addRecord(1.100, 0.7590, 0.0660, 0.0660);
        aMeasurements_102->addRecord(1.150, 0.7610, 0.0660, 0.0660);
        aMeasurements_102->addRecord(1.200, 0.7650, 0.0660, 0.0660);
        aMeasurements_102->addRecord(1.250, 0.7700, 0.0650, 0.0650);
        aMeasurements_102->addRecord(1.300, 0.7770, 0.0670, 0.0670);
        aMeasurements_102->addRecord(1.350, 0.7860, 0.0660, 0.0670);
        aMeasurements_102->addRecord(1.400, 0.7950, 0.0670, 0.0680);
        aMeasurements_102->addRecord(1.450, 0.8080, 0.0670, 0.0670);
        aMeasurements_102->addRecord(1.500, 0.8190, 0.0690, 0.0690);
        aMeasurements_102->addRecord(1.550, 0.8290, 0.0690, 0.0690);
        aMeasurements_102->addRecord(1.600, 0.8360, 0.0700, 0.0700);
        aMeasurements_102->addRecord(1.650, 0.8400, 0.0700, 0.0700);
        aMeasurements_102->addRecord(1.700, 0.8420, 0.0690, 0.0700);
        aMeasurements_102->addRecord(1.750, 0.8420, 0.0690, 0.0700);
        aMeasurements_102->addRecord(1.800, 0.8410, 0.0700, 0.0700);
        aMeasurements_102->addRecord(1.850, 0.8400, 0.0690, 0.0690);
        aMeasurements_102->addRecord(1.900, 0.8390, 0.0680, 0.0680);
        aMeasurements_102->addRecord(1.950, 0.8390, 0.0710, 0.0710);
        aMeasurements_102->addRecord(2.000, 0.8390, 0.0690, 0.0690);
        aMeasurements_102->addRecord(2.050, 0.8400, 0.0680, 0.0680);
        aMeasurements_102->addRecord(2.100, 0.8410, 0.0680, 0.0680);
        aMeasurements_102->addRecord(2.150, 0.8390, 0.0690, 0.0690);
        aMeasurements_102->addRecord(2.200, 0.8300, 0.0700, 0.0700);
        aMeasurements_102->addRecord(2.250, 0.8300, 0.0700, 0.0700);
        aMeasurements_102->addRecord(2.300, 0.8320, 0.0690, 0.0690);
        aMeasurements_102->addRecord(2.350, 0.8320, 0.0690, 0.0700);
        aMeasurements_102->addRecord(2.400, 0.8320, 0.0700, 0.0700);
        aMeasurements_102->addRecord(2.450, 0.8260, 0.0690, 0.0690);
        aMeasurements_102->addRecord(2.500, 0.8220, 0.0680, 0.0680);

        return aMeasurements_102;
    }

    std::shared_ptr<CSpectralSampleData> loadSampleData_NFRC_20243()
    {
        std::shared_ptr<CSpectralSampleData> aMeasurements_20243 =
          std::make_shared<CSpectralSampleData>();
        aMeasurements_20243->addRecord(0.300, 0.0001, 0.0467, 0.0473);
        aMeasurements_20243->addRecord(0.305, 0.0001, 0.0460, 0.0464);
        aMeasurements_20243->addRecord(0.310, 0.0005, 0.0458, 0.0464);
        aMeasurements_20243->addRecord(0.315, 0.0000, 0.0455, 0.0463);
        aMeasurements_20243->addRecord(0.320, 0.0000, 0.0451, 0.0457);
        aMeasurements_20243->addRecord(0.325, 0.0001, 0.0448, 0.0454);
        aMeasurements_20243->addRecord(0.330, 0.0003, 0.0445, 0.0454);
        aMeasurements_20243->addRecord(0.335, 0.0001, 0.0437, 0.0447);
        aMeasurements_20243->addRecord(0.340, 0.0000, 0.0438, 0.0443);
        aMeasurements_20243->addRecord(0.345, 0.0000, 0.0434, 0.0419);
        aMeasurements_20243->addRecord(0.350, 0.0000, 0.0433, 0.0443);
        aMeasurements_20243->addRecord(0.355, 0.0000, 0.0430, 0.0441);
        aMeasurements_20243->addRecord(0.360, 0.0000, 0.0427, 0.0441);
        aMeasurements_20243->addRecord(0.365, 0.0001, 0.0426, 0.0441);
        aMeasurements_20243->addRecord(0.370, 0.0000, 0.0424, 0.0441);
        aMeasurements_20243->addRecord(0.375, 0.0000, 0.0425, 0.0440);
        aMeasurements_20243->addRecord(0.380, 0.0000, 0.0424, 0.0440);
        aMeasurements_20243->addRecord(0.385, 0.0000, 0.0432, 0.0448);
        aMeasurements_20243->addRecord(0.390, 0.0002, 0.0449, 0.0474);
        aMeasurements_20243->addRecord(0.395, 0.0002, 0.0459, 0.0487);
        aMeasurements_20243->addRecord(0.400, 0.0005, 0.0454, 0.0481);
        aMeasurements_20243->addRecord(0.405, 0.0007, 0.0443, 0.0468);
        aMeasurements_20243->addRecord(0.410, 0.0008, 0.0442, 0.0461);
        aMeasurements_20243->addRecord(0.415, 0.0009, 0.0452, 0.0466);
        aMeasurements_20243->addRecord(0.420, 0.0011, 0.0470, 0.0479);
        aMeasurements_20243->addRecord(0.425, 0.0010, 0.0492, 0.0499);
        aMeasurements_20243->addRecord(0.430, 0.0012, 0.0518, 0.0523);
        aMeasurements_20243->addRecord(0.435, 0.0011, 0.0538, 0.0544);
        aMeasurements_20243->addRecord(0.440, 0.0012, 0.0557, 0.0562);
        aMeasurements_20243->addRecord(0.445, 0.0012, 0.0572, 0.0577);
        aMeasurements_20243->addRecord(0.450, 0.0012, 0.0582, 0.0588);
        aMeasurements_20243->addRecord(0.455, 0.0012, 0.0589, 0.0595);
        aMeasurements_20243->addRecord(0.460, 0.0013, 0.0594, 0.0599);
        aMeasurements_20243->addRecord(0.465, 0.0013, 0.0599, 0.0604);
        aMeasurements_20243->addRecord(0.470, 0.0013, 0.0605, 0.0609);
        aMeasurements_20243->addRecord(0.475, 0.0013, 0.0613, 0.0617);
        aMeasurements_20243->addRecord(0.480, 0.0013, 0.0624, 0.0626);
        aMeasurements_20243->addRecord(0.485, 0.0012, 0.0637, 0.0640);
        aMeasurements_20243->addRecord(0.490, 0.0012, 0.0651, 0.0653);
        aMeasurements_20243->addRecord(0.495, 0.0012, 0.0667, 0.0670);
        aMeasurements_20243->addRecord(0.500, 0.0011, 0.0682, 0.0686);
        aMeasurements_20243->addRecord(0.505, 0.0011, 0.0695, 0.0702);
        aMeasurements_20243->addRecord(0.510, 0.0010, 0.0704, 0.0714);
        aMeasurements_20243->addRecord(0.515, 0.0010, 0.0710, 0.0721);
        aMeasurements_20243->addRecord(0.520, 0.0010, 0.0710, 0.0725);
        aMeasurements_20243->addRecord(0.525, 0.0009, 0.0705, 0.0723);
        aMeasurements_20243->addRecord(0.530, 0.0009, 0.0696, 0.0717);
        aMeasurements_20243->addRecord(0.535, 0.0009, 0.0680, 0.0705);
        aMeasurements_20243->addRecord(0.540, 0.0008, 0.0662, 0.0687);
        aMeasurements_20243->addRecord(0.545, 0.0008, 0.0639, 0.0666);
        aMeasurements_20243->addRecord(0.550, 0.0007, 0.0614, 0.0642);
        aMeasurements_20243->addRecord(0.555, 0.0007, 0.0586, 0.0615);
        aMeasurements_20243->addRecord(0.560, 0.0006, 0.0560, 0.0588);
        aMeasurements_20243->addRecord(0.565, 0.0006, 0.0532, 0.0562);
        aMeasurements_20243->addRecord(0.570, 0.0006, 0.0508, 0.0535);
        aMeasurements_20243->addRecord(0.575, 0.0005, 0.0485, 0.0512);
        aMeasurements_20243->addRecord(0.580, 0.0005, 0.0466, 0.0492);
        aMeasurements_20243->addRecord(0.585, 0.0004, 0.0451, 0.0475);
        aMeasurements_20243->addRecord(0.590, 0.0005, 0.0440, 0.0464);
        aMeasurements_20243->addRecord(0.595, 0.0004, 0.0433, 0.0455);
        aMeasurements_20243->addRecord(0.600, 0.0004, 0.0431, 0.0452);
        aMeasurements_20243->addRecord(0.605, 0.0003, 0.0432, 0.0452);
        aMeasurements_20243->addRecord(0.610, 0.0003, 0.0437, 0.0456);
        aMeasurements_20243->addRecord(0.615, 0.0003, 0.0443, 0.0463);
        aMeasurements_20243->addRecord(0.620, 0.0002, 0.0452, 0.0472);
        aMeasurements_20243->addRecord(0.625, 0.0002, 0.0462, 0.0482);
        aMeasurements_20243->addRecord(0.630, 0.0002, 0.0473, 0.0493);
        aMeasurements_20243->addRecord(0.635, 0.0002, 0.0483, 0.0504);
        aMeasurements_20243->addRecord(0.640, 0.0002, 0.0492, 0.0515);
        aMeasurements_20243->addRecord(0.645, 0.0002, 0.0501, 0.0525);
        aMeasurements_20243->addRecord(0.650, 0.0001, 0.0508, 0.0533);
        aMeasurements_20243->addRecord(0.655, 0.0001, 0.0515, 0.0541);
        aMeasurements_20243->addRecord(0.660, 0.0001, 0.0519, 0.0546);
        aMeasurements_20243->addRecord(0.665, 0.0001, 0.0523, 0.0549);
        aMeasurements_20243->addRecord(0.670, 0.0001, 0.0524, 0.0552);
        aMeasurements_20243->addRecord(0.675, 0.0001, 0.0525, 0.0553);
        aMeasurements_20243->addRecord(0.680, 0.0001, 0.0524, 0.0551);
        aMeasurements_20243->addRecord(0.685, 0.0001, 0.0522, 0.0549);
        aMeasurements_20243->addRecord(0.690, 0.0001, 0.0518, 0.0545);
        aMeasurements_20243->addRecord(0.695, 0.0001, 0.0513, 0.0542);
        aMeasurements_20243->addRecord(0.700, 0.0000, 0.0510, 0.0536);
        aMeasurements_20243->addRecord(0.705, 0.0000, 0.0505, 0.0531);
        aMeasurements_20243->addRecord(0.710, 0.0001, 0.0498, 0.0523);
        aMeasurements_20243->addRecord(0.715, 0.0000, 0.0493, 0.0518);
        aMeasurements_20243->addRecord(0.720, 0.0001, 0.0489, 0.0512);
        aMeasurements_20243->addRecord(0.725, 0.0000, 0.0485, 0.0507);
        aMeasurements_20243->addRecord(0.730, 0.0001, 0.0481, 0.0502);
        aMeasurements_20243->addRecord(0.735, 0.0000, 0.0477, 0.0498);
        aMeasurements_20243->addRecord(0.740, 0.0000, 0.0474, 0.0495);
        aMeasurements_20243->addRecord(0.745, 0.0000, 0.0471, 0.0491);
        aMeasurements_20243->addRecord(0.750, 0.0000, 0.0470, 0.0485);
        aMeasurements_20243->addRecord(0.755, 0.0001, 0.0469, 0.0486);
        aMeasurements_20243->addRecord(0.760, 0.0000, 0.0470, 0.0484);
        aMeasurements_20243->addRecord(0.765, 0.0000, 0.0471, 0.0483);
        aMeasurements_20243->addRecord(0.770, 0.0000, 0.0471, 0.0484);
        aMeasurements_20243->addRecord(0.775, 0.0000, 0.0474, 0.0484);
        aMeasurements_20243->addRecord(0.780, 0.0000, 0.0473, 0.0486);
        aMeasurements_20243->addRecord(0.785, 0.0000, 0.0477, 0.0490);
        aMeasurements_20243->addRecord(0.790, 0.0000, 0.0479, 0.0492);
        aMeasurements_20243->addRecord(0.795, 0.0001, 0.0481, 0.0496);
        aMeasurements_20243->addRecord(0.800, 0.0072, 0.0492, 0.0553);
        aMeasurements_20243->addRecord(0.805, 0.0018, 0.0517, 0.0470);
        aMeasurements_20243->addRecord(0.810, 0.0000, 0.0507, 0.0514);
        aMeasurements_20243->addRecord(0.815, 0.0000, 0.0499, 0.0517);
        aMeasurements_20243->addRecord(0.820, 0.0000, 0.0536, 0.0543);
        aMeasurements_20243->addRecord(0.825, 0.0000, 0.0499, 0.0523);
        aMeasurements_20243->addRecord(0.830, 0.0000, 0.0538, 0.0555);
        aMeasurements_20243->addRecord(0.835, 0.0000, 0.0496, 0.0521);
        aMeasurements_20243->addRecord(0.840, 0.0011, 0.0518, 0.0536);
        aMeasurements_20243->addRecord(0.845, 0.0013, 0.0533, 0.0569);
        aMeasurements_20243->addRecord(0.850, 0.0000, 0.0517, 0.0563);
        aMeasurements_20243->addRecord(0.855, 0.0000, 0.0525, 0.0534);
        aMeasurements_20243->addRecord(0.860, 0.0008, 0.0542, 0.0555);
        aMeasurements_20243->addRecord(0.865, 0.0001, 0.0540, 0.0563);
        aMeasurements_20243->addRecord(0.870, 0.0000, 0.0552, 0.0582);
        aMeasurements_20243->addRecord(0.875, 0.0008, 0.0576, 0.0576);
        aMeasurements_20243->addRecord(0.880, 0.0007, 0.0577, 0.0581);
        aMeasurements_20243->addRecord(0.885, 0.0000, 0.0577, 0.0587);
        aMeasurements_20243->addRecord(0.890, 0.0000, 0.0587, 0.0601);
        aMeasurements_20243->addRecord(0.895, 0.0000, 0.0608, 0.0606);
        aMeasurements_20243->addRecord(0.900, 0.0000, 0.0627, 0.0630);
        aMeasurements_20243->addRecord(0.905, 0.0002, 0.0641, 0.0645);
        aMeasurements_20243->addRecord(0.910, 0.0000, 0.0663, 0.0656);
        aMeasurements_20243->addRecord(0.915, 0.0000, 0.0685, 0.0684);
        aMeasurements_20243->addRecord(0.920, 0.0000, 0.0715, 0.0704);
        aMeasurements_20243->addRecord(0.925, 0.0000, 0.0741, 0.0727);
        aMeasurements_20243->addRecord(0.930, 0.0002, 0.0778, 0.0758);
        aMeasurements_20243->addRecord(0.935, 0.0003, 0.0808, 0.0791);
        aMeasurements_20243->addRecord(0.940, 0.0000, 0.0844, 0.0822);
        aMeasurements_20243->addRecord(0.945, 0.0000, 0.0888, 0.0860);
        aMeasurements_20243->addRecord(0.950, 0.0000, 0.0925, 0.0892);
        aMeasurements_20243->addRecord(0.955, 0.0001, 0.0969, 0.0936);
        aMeasurements_20243->addRecord(0.960, 0.0000, 0.1017, 0.0978);
        aMeasurements_20243->addRecord(0.965, 0.0000, 0.1063, 0.1020);
        aMeasurements_20243->addRecord(0.970, 0.0000, 0.1112, 0.1068);
        aMeasurements_20243->addRecord(0.975, 0.0002, 0.1167, 0.1118);
        aMeasurements_20243->addRecord(0.980, 0.0000, 0.1219, 0.1167);
        aMeasurements_20243->addRecord(0.985, 0.0000, 0.1273, 0.1217);
        aMeasurements_20243->addRecord(0.990, 0.0001, 0.1334, 0.1274);
        aMeasurements_20243->addRecord(0.995, 0.0000, 0.1392, 0.1329);
        aMeasurements_20243->addRecord(1.000, 0.0001, 0.1449, 0.1387);
        aMeasurements_20243->addRecord(1.005, 0.0000, 0.1510, 0.1444);
        aMeasurements_20243->addRecord(1.010, 0.0000, 0.1568, 0.1499);
        aMeasurements_20243->addRecord(1.015, 0.0000, 0.1629, 0.1559);
        aMeasurements_20243->addRecord(1.020, 0.0000, 0.1690, 0.1619);
        aMeasurements_20243->addRecord(1.025, 0.0000, 0.1750, 0.1676);
        aMeasurements_20243->addRecord(1.030, 0.0000, 0.1811, 0.1738);
        aMeasurements_20243->addRecord(1.035, 0.0002, 0.1871, 0.1798);
        aMeasurements_20243->addRecord(1.040, 0.0000, 0.1932, 0.1857);
        aMeasurements_20243->addRecord(1.045, 0.0000, 0.1993, 0.1916);
        aMeasurements_20243->addRecord(1.050, 0.0000, 0.2054, 0.1976);
        aMeasurements_20243->addRecord(1.055, 0.0000, 0.2111, 0.2033);
        aMeasurements_20243->addRecord(1.060, 0.0000, 0.2167, 0.2090);
        aMeasurements_20243->addRecord(1.065, 0.0000, 0.2224, 0.2147);
        aMeasurements_20243->addRecord(1.070, 0.0000, 0.2280, 0.2202);
        aMeasurements_20243->addRecord(1.075, 0.0000, 0.2333, 0.2257);
        aMeasurements_20243->addRecord(1.080, 0.0000, 0.2387, 0.2310);
        aMeasurements_20243->addRecord(1.085, 0.0000, 0.2440, 0.2363);
        aMeasurements_20243->addRecord(1.090, 0.0000, 0.2489, 0.2415);
        aMeasurements_20243->addRecord(1.095, 0.0000, 0.2539, 0.2465);
        aMeasurements_20243->addRecord(1.100, 0.0000, 0.2587, 0.2514);
        aMeasurements_20243->addRecord(1.105, 0.0000, 0.2633, 0.2561);
        aMeasurements_20243->addRecord(1.110, 0.0001, 0.2677, 0.2607);
        aMeasurements_20243->addRecord(1.115, 0.0001, 0.2719, 0.2650);
        aMeasurements_20243->addRecord(1.120, 0.0000, 0.2758, 0.2691);
        aMeasurements_20243->addRecord(1.125, 0.0000, 0.2795, 0.2729);
        aMeasurements_20243->addRecord(1.130, 0.0000, 0.2828, 0.2765);
        aMeasurements_20243->addRecord(1.135, 0.0001, 0.2852, 0.2793);
        aMeasurements_20243->addRecord(1.140, 0.0000, 0.2872, 0.2814);
        aMeasurements_20243->addRecord(1.145, 0.0000, 0.2887, 0.2833);
        aMeasurements_20243->addRecord(1.150, 0.0002, 0.2897, 0.2848);
        aMeasurements_20243->addRecord(1.155, 0.0000, 0.2906, 0.2862);
        aMeasurements_20243->addRecord(1.160, 0.0000, 0.2910, 0.2867);
        aMeasurements_20243->addRecord(1.165, 0.0000, 0.2904, 0.2867);
        aMeasurements_20243->addRecord(1.170, 0.0000, 0.2888, 0.2855);
        aMeasurements_20243->addRecord(1.175, 0.0000, 0.2858, 0.2833);
        aMeasurements_20243->addRecord(1.180, 0.0000, 0.2823, 0.2803);
        aMeasurements_20243->addRecord(1.185, 0.0000, 0.2788, 0.2775);
        aMeasurements_20243->addRecord(1.190, 0.0000, 0.2768, 0.2761);
        aMeasurements_20243->addRecord(1.195, 0.0000, 0.2771, 0.2766);
        aMeasurements_20243->addRecord(1.200, 0.0000, 0.2802, 0.2800);
        aMeasurements_20243->addRecord(1.205, 0.0001, 0.2846, 0.2843);
        aMeasurements_20243->addRecord(1.210, 0.0000, 0.2902, 0.2898);
        aMeasurements_20243->addRecord(1.215, 0.0000, 0.2959, 0.2951);
        aMeasurements_20243->addRecord(1.220, 0.0000, 0.3010, 0.3000);
        aMeasurements_20243->addRecord(1.225, 0.0000, 0.3057, 0.3041);
        aMeasurements_20243->addRecord(1.230, 0.0000, 0.3097, 0.3074);
        aMeasurements_20243->addRecord(1.235, 0.0000, 0.3133, 0.3103);
        aMeasurements_20243->addRecord(1.240, 0.0001, 0.3166, 0.3128);
        aMeasurements_20243->addRecord(1.245, 0.0000, 0.3196, 0.3150);
        aMeasurements_20243->addRecord(1.250, 0.0000, 0.3229, 0.3173);
        aMeasurements_20243->addRecord(1.255, 0.0000, 0.3261, 0.3195);
        aMeasurements_20243->addRecord(1.260, 0.0001, 0.3293, 0.3218);
        aMeasurements_20243->addRecord(1.265, 0.0002, 0.3326, 0.3245);
        aMeasurements_20243->addRecord(1.270, 0.0001, 0.3359, 0.3272);
        aMeasurements_20243->addRecord(1.275, 0.0000, 0.3390, 0.3301);
        aMeasurements_20243->addRecord(1.280, 0.0000, 0.3424, 0.3332);
        aMeasurements_20243->addRecord(1.285, 0.0000, 0.3456, 0.3366);
        aMeasurements_20243->addRecord(1.290, 0.0000, 0.3484, 0.3395);
        aMeasurements_20243->addRecord(1.295, 0.0000, 0.3510, 0.3424);
        aMeasurements_20243->addRecord(1.300, 0.0000, 0.3534, 0.3450);
        aMeasurements_20243->addRecord(1.305, 0.0000, 0.3555, 0.3475);
        aMeasurements_20243->addRecord(1.310, 0.0000, 0.3573, 0.3498);
        aMeasurements_20243->addRecord(1.315, 0.0000, 0.3591, 0.3517);
        aMeasurements_20243->addRecord(1.320, 0.0000, 0.3605, 0.3534);
        aMeasurements_20243->addRecord(1.325, 0.0000, 0.3617, 0.3548);
        aMeasurements_20243->addRecord(1.330, 0.0000, 0.3626, 0.3559);
        aMeasurements_20243->addRecord(1.335, 0.0000, 0.3631, 0.3566);
        aMeasurements_20243->addRecord(1.340, 0.0000, 0.3629, 0.3567);
        aMeasurements_20243->addRecord(1.345, 0.0000, 0.3627, 0.3568);
        aMeasurements_20243->addRecord(1.350, 0.0000, 0.3622, 0.3563);
        aMeasurements_20243->addRecord(1.355, 0.0000, 0.3600, 0.3544);
        aMeasurements_20243->addRecord(1.360, 0.0002, 0.3570, 0.3519);
        aMeasurements_20243->addRecord(1.365, 0.0001, 0.3534, 0.3484);
        aMeasurements_20243->addRecord(1.370, 0.0000, 0.3493, 0.3447);
        aMeasurements_20243->addRecord(1.375, 0.0000, 0.3456, 0.3416);
        aMeasurements_20243->addRecord(1.380, 0.0000, 0.3425, 0.3389);
        aMeasurements_20243->addRecord(1.385, 0.0000, 0.3401, 0.3369);
        aMeasurements_20243->addRecord(1.390, 0.0000, 0.3385, 0.3356);
        aMeasurements_20243->addRecord(1.395, 0.0000, 0.3368, 0.3344);
        aMeasurements_20243->addRecord(1.400, 0.0000, 0.3347, 0.3327);
        aMeasurements_20243->addRecord(1.405, 0.0000, 0.3321, 0.3304);
        aMeasurements_20243->addRecord(1.410, 0.0000, 0.3291, 0.3279);
        aMeasurements_20243->addRecord(1.415, 0.0000, 0.3264, 0.3255);
        aMeasurements_20243->addRecord(1.420, 0.0000, 0.3238, 0.3233);
        aMeasurements_20243->addRecord(1.425, 0.0000, 0.3220, 0.3218);
        aMeasurements_20243->addRecord(1.430, 0.0000, 0.3208, 0.3208);
        aMeasurements_20243->addRecord(1.435, 0.0001, 0.3204, 0.3203);
        aMeasurements_20243->addRecord(1.440, 0.0000, 0.3208, 0.3208);
        aMeasurements_20243->addRecord(1.445, 0.0000, 0.3226, 0.3227);
        aMeasurements_20243->addRecord(1.450, 0.0000, 0.3257, 0.3256);
        aMeasurements_20243->addRecord(1.455, 0.0000, 0.3301, 0.3300);
        aMeasurements_20243->addRecord(1.460, 0.0000, 0.3356, 0.3353);
        aMeasurements_20243->addRecord(1.465, 0.0000, 0.3420, 0.3413);
        aMeasurements_20243->addRecord(1.470, 0.0000, 0.3487, 0.3479);
        aMeasurements_20243->addRecord(1.475, 0.0000, 0.3557, 0.3546);
        aMeasurements_20243->addRecord(1.480, 0.0001, 0.3627, 0.3612);
        aMeasurements_20243->addRecord(1.485, 0.0000, 0.3694, 0.3676);
        aMeasurements_20243->addRecord(1.490, 0.0000, 0.3758, 0.3737);
        aMeasurements_20243->addRecord(1.495, 0.0000, 0.3820, 0.3796);
        aMeasurements_20243->addRecord(1.500, 0.0000, 0.3880, 0.3853);
        aMeasurements_20243->addRecord(1.505, 0.0000, 0.3938, 0.3910);
        aMeasurements_20243->addRecord(1.510, 0.0000, 0.3994, 0.3962);
        aMeasurements_20243->addRecord(1.515, 0.0000, 0.4045, 0.4012);
        aMeasurements_20243->addRecord(1.520, 0.0000, 0.4092, 0.4059);
        aMeasurements_20243->addRecord(1.525, 0.0000, 0.4135, 0.4102);
        aMeasurements_20243->addRecord(1.530, 0.0000, 0.4174, 0.4139);
        aMeasurements_20243->addRecord(1.535, 0.0000, 0.4208, 0.4172);
        aMeasurements_20243->addRecord(1.540, 0.0000, 0.4238, 0.4203);
        aMeasurements_20243->addRecord(1.545, 0.0000, 0.4266, 0.4230);
        aMeasurements_20243->addRecord(1.550, 0.0001, 0.4291, 0.4258);
        aMeasurements_20243->addRecord(1.555, 0.0000, 0.4320, 0.4287);
        aMeasurements_20243->addRecord(1.560, 0.0000, 0.4352, 0.4321);
        aMeasurements_20243->addRecord(1.565, 0.0000, 0.4387, 0.4355);
        aMeasurements_20243->addRecord(1.570, 0.0001, 0.4428, 0.4394);
        aMeasurements_20243->addRecord(1.575, 0.0001, 0.4469, 0.4437);
        aMeasurements_20243->addRecord(1.580, 0.0000, 0.4516, 0.4483);
        aMeasurements_20243->addRecord(1.585, 0.0001, 0.4560, 0.4527);
        aMeasurements_20243->addRecord(1.590, 0.0002, 0.4604, 0.4569);
        aMeasurements_20243->addRecord(1.595, 0.0000, 0.4645, 0.4610);
        aMeasurements_20243->addRecord(1.600, 0.0001, 0.4685, 0.4651);
        aMeasurements_20243->addRecord(1.605, 0.0000, 0.4720, 0.4687);
        aMeasurements_20243->addRecord(1.610, 0.0001, 0.4754, 0.4721);
        aMeasurements_20243->addRecord(1.615, 0.0000, 0.4788, 0.4753);
        aMeasurements_20243->addRecord(1.620, 0.0000, 0.4815, 0.4783);
        aMeasurements_20243->addRecord(1.625, 0.0000, 0.4842, 0.4809);
        aMeasurements_20243->addRecord(1.630, 0.0000, 0.4869, 0.4837);
        aMeasurements_20243->addRecord(1.635, 0.0002, 0.4896, 0.4863);
        aMeasurements_20243->addRecord(1.640, 0.0000, 0.4920, 0.4889);
        aMeasurements_20243->addRecord(1.645, 0.0000, 0.4940, 0.4912);
        aMeasurements_20243->addRecord(1.650, 0.0000, 0.4955, 0.4925);
        aMeasurements_20243->addRecord(1.655, 0.0000, 0.4958, 0.4932);
        aMeasurements_20243->addRecord(1.660, 0.0001, 0.4940, 0.4920);
        aMeasurements_20243->addRecord(1.665, 0.0000, 0.4894, 0.4878);
        aMeasurements_20243->addRecord(1.670, 0.0000, 0.4798, 0.4789);
        aMeasurements_20243->addRecord(1.675, 0.0000, 0.4619, 0.4623);
        aMeasurements_20243->addRecord(1.680, 0.0003, 0.4352, 0.4374);
        aMeasurements_20243->addRecord(1.685, 0.0002, 0.3990, 0.4033);
        aMeasurements_20243->addRecord(1.690, 0.0001, 0.3573, 0.3634);
        aMeasurements_20243->addRecord(1.695, 0.0000, 0.3175, 0.3256);
        aMeasurements_20243->addRecord(1.700, 0.0004, 0.2840, 0.2936);
        aMeasurements_20243->addRecord(1.705, 0.0000, 0.2626, 0.2731);
        aMeasurements_20243->addRecord(1.710, 0.0001, 0.2547, 0.2655);
        aMeasurements_20243->addRecord(1.715, 0.0000, 0.2560, 0.2669);
        aMeasurements_20243->addRecord(1.720, 0.0000, 0.2625, 0.2734);
        aMeasurements_20243->addRecord(1.725, 0.0002, 0.2696, 0.2804);
        aMeasurements_20243->addRecord(1.730, 0.0000, 0.2735, 0.2845);
        aMeasurements_20243->addRecord(1.735, 0.0000, 0.2727, 0.2840);
        aMeasurements_20243->addRecord(1.740, 0.0002, 0.2709, 0.2819);
        aMeasurements_20243->addRecord(1.745, 0.0000, 0.2697, 0.2812);
        aMeasurements_20243->addRecord(1.750, 0.0001, 0.2732, 0.2851);
        aMeasurements_20243->addRecord(1.755, 0.0000, 0.2836, 0.2955);
        aMeasurements_20243->addRecord(1.760, 0.0000, 0.3001, 0.3114);
        aMeasurements_20243->addRecord(1.765, 0.0000, 0.3190, 0.3302);
        aMeasurements_20243->addRecord(1.770, 0.0000, 0.3380, 0.3490);
        aMeasurements_20243->addRecord(1.775, 0.0001, 0.3551, 0.3657);
        aMeasurements_20243->addRecord(1.780, 0.0000, 0.3684, 0.3787);
        aMeasurements_20243->addRecord(1.785, 0.0002, 0.3783, 0.3888);
        aMeasurements_20243->addRecord(1.790, 0.0001, 0.3859, 0.3961);
        aMeasurements_20243->addRecord(1.795, 0.0000, 0.3911, 0.4015);
        aMeasurements_20243->addRecord(1.800, 0.0006, 0.3951, 0.4056);
        aMeasurements_20243->addRecord(1.805, 0.0000, 0.3977, 0.4080);
        aMeasurements_20243->addRecord(1.810, 0.0000, 0.3996, 0.4101);
        aMeasurements_20243->addRecord(1.815, 0.0001, 0.4012, 0.4118);
        aMeasurements_20243->addRecord(1.820, 0.0000, 0.4031, 0.4134);
        aMeasurements_20243->addRecord(1.825, 0.0000, 0.4057, 0.4162);
        aMeasurements_20243->addRecord(1.830, 0.0001, 0.4106, 0.4208);
        aMeasurements_20243->addRecord(1.835, 0.0002, 0.4164, 0.4266);
        aMeasurements_20243->addRecord(1.840, 0.0010, 0.4247, 0.4347);
        aMeasurements_20243->addRecord(1.845, 0.0003, 0.4337, 0.4433);
        aMeasurements_20243->addRecord(1.850, 0.0001, 0.4423, 0.4520);
        aMeasurements_20243->addRecord(1.855, 0.0005, 0.4506, 0.4604);
        aMeasurements_20243->addRecord(1.860, 0.0000, 0.4584, 0.4677);
        aMeasurements_20243->addRecord(1.865, 0.0002, 0.4639, 0.4737);
        aMeasurements_20243->addRecord(1.870, 0.0003, 0.4674, 0.4783);
        aMeasurements_20243->addRecord(1.875, 0.0000, 0.4690, 0.4808);
        aMeasurements_20243->addRecord(1.880, 0.0001, 0.4672, 0.4814);
        aMeasurements_20243->addRecord(1.885, 0.0000, 0.4635, 0.4803);
        aMeasurements_20243->addRecord(1.890, 0.0000, 0.4576, 0.4776);
        aMeasurements_20243->addRecord(1.895, 0.0000, 0.4508, 0.4743);
        aMeasurements_20243->addRecord(1.900, 0.0000, 0.4436, 0.4710);
        aMeasurements_20243->addRecord(1.905, 0.0008, 0.4372, 0.4675);
        aMeasurements_20243->addRecord(1.910, 0.0000, 0.4323, 0.4651);
        aMeasurements_20243->addRecord(1.915, 0.0000, 0.4283, 0.4632);
        aMeasurements_20243->addRecord(1.920, 0.0000, 0.4266, 0.4626);
        aMeasurements_20243->addRecord(1.925, 0.0000, 0.4262, 0.4625);
        aMeasurements_20243->addRecord(1.930, 0.0000, 0.4278, 0.4632);
        aMeasurements_20243->addRecord(1.935, 0.0004, 0.4299, 0.4646);
        aMeasurements_20243->addRecord(1.940, 0.0000, 0.4340, 0.4669);
        aMeasurements_20243->addRecord(1.945, 0.0000, 0.4397, 0.4710);
        aMeasurements_20243->addRecord(1.950, 0.0000, 0.4450, 0.4755);
        aMeasurements_20243->addRecord(1.955, 0.0001, 0.4527, 0.4810);
        aMeasurements_20243->addRecord(1.960, 0.0002, 0.4586, 0.4854);
        aMeasurements_20243->addRecord(1.965, 0.0006, 0.4639, 0.4895);
        aMeasurements_20243->addRecord(1.970, 0.0000, 0.4686, 0.4926);
        aMeasurements_20243->addRecord(1.975, 0.0000, 0.4722, 0.4945);
        aMeasurements_20243->addRecord(1.980, 0.0002, 0.4744, 0.4959);
        aMeasurements_20243->addRecord(1.985, 0.0003, 0.4766, 0.4961);
        aMeasurements_20243->addRecord(1.990, 0.0007, 0.4761, 0.4958);
        aMeasurements_20243->addRecord(1.995, 0.0000, 0.4748, 0.4922);
        aMeasurements_20243->addRecord(2.000, 0.0000, 0.4692, 0.4865);
        aMeasurements_20243->addRecord(2.005, 0.0000, 0.4605, 0.4770);
        aMeasurements_20243->addRecord(2.010, 0.0018, 0.4496, 0.4665);
        aMeasurements_20243->addRecord(2.015, 0.0000, 0.4372, 0.4530);
        aMeasurements_20243->addRecord(2.020, 0.0000, 0.4250, 0.4412);
        aMeasurements_20243->addRecord(2.025, 0.0007, 0.4130, 0.4289);
        aMeasurements_20243->addRecord(2.030, 0.0000, 0.4027, 0.4187);
        aMeasurements_20243->addRecord(2.035, 0.0014, 0.3935, 0.4089);
        aMeasurements_20243->addRecord(2.040, 0.0000, 0.3837, 0.3989);
        aMeasurements_20243->addRecord(2.045, 0.0000, 0.3749, 0.3907);
        aMeasurements_20243->addRecord(2.050, 0.0000, 0.3649, 0.3816);
        aMeasurements_20243->addRecord(2.055, 0.0027, 0.3564, 0.3728);
        aMeasurements_20243->addRecord(2.060, 0.0000, 0.3481, 0.3658);
        aMeasurements_20243->addRecord(2.065, 0.0000, 0.3436, 0.3598);
        aMeasurements_20243->addRecord(2.070, 0.0000, 0.3388, 0.3567);
        aMeasurements_20243->addRecord(2.075, 0.0000, 0.3389, 0.3565);
        aMeasurements_20243->addRecord(2.080, 0.0004, 0.3400, 0.3578);
        aMeasurements_20243->addRecord(2.085, 0.0006, 0.3430, 0.3611);
        aMeasurements_20243->addRecord(2.090, 0.0021, 0.3488, 0.3668);
        aMeasurements_20243->addRecord(2.095, 0.0000, 0.3562, 0.3729);
        aMeasurements_20243->addRecord(2.100, 0.0012, 0.3616, 0.3790);
        aMeasurements_20243->addRecord(2.105, 0.0008, 0.3647, 0.3851);
        aMeasurements_20243->addRecord(2.110, 0.0011, 0.3711, 0.3900);
        aMeasurements_20243->addRecord(2.115, 0.0000, 0.3767, 0.3945);
        aMeasurements_20243->addRecord(2.120, 0.0014, 0.3805, 0.3975);
        aMeasurements_20243->addRecord(2.125, 0.0000, 0.3839, 0.4005);
        aMeasurements_20243->addRecord(2.130, 0.0000, 0.3854, 0.4013);
        aMeasurements_20243->addRecord(2.135, 0.0030, 0.3880, 0.4049);
        aMeasurements_20243->addRecord(2.140, 0.0000, 0.3917, 0.4075);
        aMeasurements_20243->addRecord(2.145, 0.0000, 0.3923, 0.4098);
        aMeasurements_20243->addRecord(2.150, 0.0000, 0.3941, 0.4101);
        aMeasurements_20243->addRecord(2.155, 0.0041, 0.3952, 0.4147);
        aMeasurements_20243->addRecord(2.160, 0.0000, 0.3979, 0.4149);
        aMeasurements_20243->addRecord(2.165, 0.0008, 0.3977, 0.4148);
        aMeasurements_20243->addRecord(2.170, 0.0000, 0.3960, 0.4146);
        aMeasurements_20243->addRecord(2.175, 0.0003, 0.3970, 0.4149);
        aMeasurements_20243->addRecord(2.180, 0.0039, 0.3970, 0.4140);
        aMeasurements_20243->addRecord(2.185, 0.0000, 0.3950, 0.4121);
        aMeasurements_20243->addRecord(2.190, 0.0000, 0.3917, 0.4091);
        aMeasurements_20243->addRecord(2.195, 0.0000, 0.3873, 0.4050);
        aMeasurements_20243->addRecord(2.200, 0.0008, 0.3843, 0.4011);
        aMeasurements_20243->addRecord(2.205, 0.0008, 0.3791, 0.3958);
        aMeasurements_20243->addRecord(2.210, 0.0007, 0.3744, 0.3917);
        aMeasurements_20243->addRecord(2.215, 0.0009, 0.3678, 0.3841);
        aMeasurements_20243->addRecord(2.220, 0.0006, 0.3591, 0.3771);
        aMeasurements_20243->addRecord(2.225, 0.0007, 0.3498, 0.3671);
        aMeasurements_20243->addRecord(2.230, 0.0013, 0.3356, 0.3535);
        aMeasurements_20243->addRecord(2.235, 0.0000, 0.3164, 0.3349);
        aMeasurements_20243->addRecord(2.240, 0.0003, 0.2935, 0.3115);
        aMeasurements_20243->addRecord(2.245, 0.0000, 0.2623, 0.2799);
        aMeasurements_20243->addRecord(2.250, 0.0010, 0.2264, 0.2423);
        aMeasurements_20243->addRecord(2.255, 0.0009, 0.1851, 0.2018);
        aMeasurements_20243->addRecord(2.260, 0.0000, 0.1463, 0.1599);
        aMeasurements_20243->addRecord(2.265, 0.0018, 0.1121, 0.1239);
        aMeasurements_20243->addRecord(2.270, 0.0014, 0.0843, 0.0948);
        aMeasurements_20243->addRecord(2.275, 0.0016, 0.0659, 0.0732);
        aMeasurements_20243->addRecord(2.280, 0.0004, 0.0551, 0.0612);
        aMeasurements_20243->addRecord(2.285, 0.0000, 0.0479, 0.0533);
        aMeasurements_20243->addRecord(2.290, 0.0000, 0.0446, 0.0476);
        aMeasurements_20243->addRecord(2.295, 0.0000, 0.0439, 0.0465);
        aMeasurements_20243->addRecord(2.300, 0.0006, 0.0434, 0.0457);
        aMeasurements_20243->addRecord(2.305, 0.0000, 0.0439, 0.0452);
        aMeasurements_20243->addRecord(2.310, 0.0011, 0.0427, 0.0468);
        aMeasurements_20243->addRecord(2.315, 0.0000, 0.0451, 0.0463);
        aMeasurements_20243->addRecord(2.320, 0.0000, 0.0454, 0.0485);
        aMeasurements_20243->addRecord(2.325, 0.0019, 0.0473, 0.0510);
        aMeasurements_20243->addRecord(2.330, 0.0000, 0.0478, 0.0513);
        aMeasurements_20243->addRecord(2.335, 0.0000, 0.0484, 0.0526);
        aMeasurements_20243->addRecord(2.340, 0.0000, 0.0483, 0.0529);
        aMeasurements_20243->addRecord(2.345, 0.0000, 0.0477, 0.0530);
        aMeasurements_20243->addRecord(2.350, 0.0000, 0.0503, 0.0533);
        aMeasurements_20243->addRecord(2.355, 0.0000, 0.0506, 0.0548);
        aMeasurements_20243->addRecord(2.360, 0.0000, 0.0527, 0.0570);
        aMeasurements_20243->addRecord(2.365, 0.0002, 0.0526, 0.0568);
        aMeasurements_20243->addRecord(2.370, 0.0000, 0.0521, 0.0572);
        aMeasurements_20243->addRecord(2.375, 0.0003, 0.0523, 0.0576);
        aMeasurements_20243->addRecord(2.380, 0.0000, 0.0521, 0.0582);
        aMeasurements_20243->addRecord(2.385, 0.0000, 0.0518, 0.0570);
        aMeasurements_20243->addRecord(2.390, 0.0011, 0.0548, 0.0607);
        aMeasurements_20243->addRecord(2.395, 0.0020, 0.0546, 0.0615);
        aMeasurements_20243->addRecord(2.400, 0.0000, 0.0559, 0.0607);
        aMeasurements_20243->addRecord(2.405, 0.0003, 0.0566, 0.0637);
        aMeasurements_20243->addRecord(2.410, 0.0000, 0.0578, 0.0628);
        aMeasurements_20243->addRecord(2.415, 0.0000, 0.0574, 0.0625);
        aMeasurements_20243->addRecord(2.420, 0.0012, 0.0547, 0.0610);
        aMeasurements_20243->addRecord(2.425, 0.0000, 0.0527, 0.0575);
        aMeasurements_20243->addRecord(2.430, 0.0000, 0.0538, 0.0584);
        aMeasurements_20243->addRecord(2.435, 0.0000, 0.0507, 0.0567);
        aMeasurements_20243->addRecord(2.440, 0.0017, 0.0492, 0.0539);
        aMeasurements_20243->addRecord(2.445, 0.0013, 0.0498, 0.0538);
        aMeasurements_20243->addRecord(2.450, 0.0000, 0.0463, 0.0504);
        aMeasurements_20243->addRecord(2.455, 0.0000, 0.0466, 0.0507);
        aMeasurements_20243->addRecord(2.460, 0.0011, 0.0466, 0.0495);
        aMeasurements_20243->addRecord(2.465, 0.0000, 0.0457, 0.0496);
        aMeasurements_20243->addRecord(2.470, 0.0010, 0.0456, 0.0486);
        aMeasurements_20243->addRecord(2.475, 0.0013, 0.0466, 0.0526);
        aMeasurements_20243->addRecord(2.480, 0.0000, 0.0438, 0.0497);
        aMeasurements_20243->addRecord(2.485, 0.0000, 0.0475, 0.0522);
        aMeasurements_20243->addRecord(2.490, 0.0000, 0.0516, 0.0528);
        aMeasurements_20243->addRecord(2.495, 0.0039, 0.0531, 0.0553);
        aMeasurements_20243->addRecord(2.500, 0.0000, 0.0505, 0.0572);

        return aMeasurements_20243;
    }

protected:
    virtual void SetUp()
    {
        // Create material from samples
        double thickness = 3.048e-3;   // [m]
        auto aMaterial_102 = SingleLayerOptics::Material::nBandMaterial(
          loadSampleData_NFRC_102(), thickness, MaterialType::Monolithic, WavelengthRange::Solar);
        thickness = 5.715e-3;   // [m]
        auto aMaterial_20243 = SingleLayerOptics::Material::nBandMaterial(
          loadSampleData_NFRC_20243(), thickness, MaterialType::Monolithic, WavelengthRange::Solar);

        // BSDF definition is needed as well as its material representation
        auto aBSDF = std::make_shared<CBSDFHemisphere>(BSDFBasis::Quarter);
        auto Layer_102 = CBSDFLayerMaker::getSpecularLayer(aMaterial_102, aBSDF);
        auto Layer_20243 = CBSDFLayerMaker::getSpecularLayer(aMaterial_20243, aBSDF);

        // To assure interpolation to common wavelengths. MultiBSDF will NOT work with different
        // wavelengths
        CCommonWavelengths aCommonWL;
        aCommonWL.addWavelength(Layer_102->getBandWavelengths());
        aCommonWL.addWavelength(Layer_20243->getBandWavelengths());

        std::vector<double> commonWavelengths =
          aCommonWL.getCombinedWavelengths(Combine::Interpolate);

        // Equivalent BSDF layer
        std::shared_ptr<CEquivalentBSDFLayer> aEqLayer =
          std::make_shared<CEquivalentBSDFLayer>(commonWavelengths, Layer_20243);
        aEqLayer->addLayer(Layer_102);

        m_Layer = std::make_shared<CMultiPaneBSDF>(
          Layer_20243, commonWavelengths, loadSolarRadiationFile());
        m_Layer->addLayer(Layer_102);
    }

public:
    std::shared_ptr<CMultiPaneBSDF> getLayer()
    {
        return m_Layer;
    };
};

TEST_F(MultiPaneBSDF_20243_103, TestSpecular1)
{
    SCOPED_TRACE("Begin Test: Specular layer - BSDF.");

    const double minLambda = 0.3;
    const double maxLambda = 2.5;

    CMultiPaneBSDF & aLayer = *getLayer();

    double tauDiff = aLayer.DiffDiff(minLambda, maxLambda, Side::Front, PropertySimple::T);
    EXPECT_NEAR(0.000164, tauDiff, 1e-6);

    double rhoDiff = aLayer.DiffDiff(minLambda, maxLambda, Side::Front, PropertySimple::R);
    EXPECT_NEAR(0.179067, rhoDiff, 1e-6);

    double absDiff1 = aLayer.AbsDiff(minLambda, maxLambda, Side::Front, 1);
    EXPECT_NEAR(0.820757, absDiff1, 1e-6);

    double absDiff2 = aLayer.AbsDiff(minLambda, maxLambda, Side::Front, 2);
    EXPECT_NEAR(0.0000116, absDiff2, 1e-6);

    double theta = 0;
    double phi = 0;

    double tauHem = aLayer.DirHem(minLambda, maxLambda, Side::Front, PropertySimple::T, theta, phi);
    EXPECT_NEAR(0.000327, tauHem, 1e-6);

    double tauDir = aLayer.DirDir(minLambda, maxLambda, Side::Front, PropertySimple::T, theta, phi);
    EXPECT_NEAR(0.000327, tauDir, 1e-6);

    double rhoHem = aLayer.DirHem(minLambda, maxLambda, Side::Front, PropertySimple::R, theta, phi);
    EXPECT_NEAR(0.131407, rhoHem, 1e-6);

    double rhoDir = aLayer.DirDir(minLambda, maxLambda, Side::Front, PropertySimple::R, theta, phi);
    EXPECT_NEAR(0.131407, rhoDir, 1e-6);

    double abs1 = aLayer.Abs(minLambda, maxLambda, Side::Front, 1, theta, phi);
    EXPECT_NEAR(0.868248, abs1, 1e-6);

    double abs2 = aLayer.Abs(minLambda, maxLambda, Side::Front, 2, theta, phi);
    EXPECT_NEAR(0.000019, abs2, 1e-6);

    theta = 45;
    phi = 78;

    tauHem = aLayer.DirHem(minLambda, maxLambda, Side::Front, PropertySimple::T, theta, phi);
    EXPECT_NEAR(0.000205, tauHem, 1e-6);

    rhoHem = aLayer.DirHem(minLambda, maxLambda, Side::Front, PropertySimple::R, theta, phi);
    EXPECT_NEAR(0.134284, rhoHem, 1e-6);

    abs1 = aLayer.Abs(minLambda, maxLambda, Side::Front, 1, theta, phi);
    EXPECT_NEAR(0.865498, abs1, 1e-6);

    abs2 = aLayer.Abs(minLambda, maxLambda, Side::Front, 2, theta, phi);
    EXPECT_NEAR(0.000013, abs2, 1e-6);
}
