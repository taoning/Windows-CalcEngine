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

// Example on how to create multilayer BSDF with combination of specular and perforated layer

class MultiPaneBSDF_102_PerforatedCircular_NFRC18000 : public testing::Test
{
private:
    std::unique_ptr<CMultiPaneBSDF> m_Layer;

    CSeries loadSolarRadiationFile()
    {
        // Full ASTM E891-87 Table 1 (Solar radiation)
        CSeries aSolarRadiation(
          {{0.3000, 0.0},    {0.3050, 3.4},    {0.3100, 15.6},   {0.3150, 41.1},   {0.3200, 71.2},
           {0.3250, 100.2},  {0.3300, 152.4},  {0.3350, 155.6},  {0.3400, 179.4},  {0.3450, 186.7},
           {0.3500, 212.0},  {0.3600, 240.5},  {0.3700, 324.0},  {0.3800, 362.4},  {0.3900, 381.7},
           {0.4000, 556.0},  {0.4100, 656.3},  {0.4200, 690.8},  {0.4300, 641.9},  {0.4400, 798.5},
           {0.4500, 956.6},  {0.4600, 990.0},  {0.4700, 998.0},  {0.4800, 1046.1}, {0.4900, 1005.1},
           {0.5000, 1026.7}, {0.5100, 1066.7}, {0.5200, 1011.5}, {0.5300, 1084.9}, {0.5400, 1082.4},
           {0.5500, 1102.2}, {0.5700, 1087.4}, {0.5900, 1024.3}, {0.6100, 1088.8}, {0.6300, 1062.1},
           {0.6500, 1061.7}, {0.6700, 1046.2}, {0.6900, 859.2},  {0.7100, 1002.4}, {0.7180, 816.9},
           {0.7244, 842.8},  {0.7400, 971.0},  {0.7525, 956.3},  {0.7575, 942.2},  {0.7625, 524.8},
           {0.7675, 830.7},  {0.7800, 908.9},  {0.8000, 873.4},  {0.8160, 712.0},  {0.8237, 660.2},
           {0.8315, 765.5},  {0.8400, 799.8},  {0.8600, 815.2},  {0.8800, 778.3},  {0.9050, 630.4},
           {0.9150, 565.2},  {0.9250, 586.4},  {0.9300, 348.1},  {0.9370, 224.2},  {0.9480, 271.4},
           {0.9650, 451.2},  {0.9800, 549.7},  {0.9935, 630.1},  {1.0400, 582.9},  {1.0700, 539.7},
           {1.1000, 366.2},  {1.1200, 98.1},   {1.1300, 169.5},  {1.1370, 118.7},  {1.1610, 301.9},
           {1.1800, 406.8},  {1.2000, 375.2},  {1.2350, 423.6},  {1.2900, 365.7},  {1.3200, 223.4},
           {1.3500, 30.1},   {1.3950, 1.4},    {1.4425, 51.6},   {1.4625, 97.0},   {1.4770, 97.3},
           {1.4970, 167.1},  {1.5200, 239.3},  {1.5390, 248.8},  {1.5580, 249.3},  {1.5780, 222.3},
           {1.5920, 227.3},  {1.6100, 210.5},  {1.6300, 224.7},  {1.6460, 215.9},  {1.6780, 202.8},
           {1.7400, 158.2},  {1.8000, 28.6},   {1.8600, 1.8},    {1.9200, 1.1},    {1.9600, 19.7},
           {1.9850, 84.9},   {2.0050, 25.0},   {2.0350, 92.5},   {2.0650, 56.3},   {2.1000, 82.7},
           {2.1480, 76.2},   {2.1980, 66.4},   {2.2700, 65.0},   {2.3600, 57.6},   {2.4500, 19.8},
           {2.4940, 17.0},   {2.5370, 3.0},    {2.9410, 4.0},    {2.9730, 7.0},    {3.0050, 6.0},
           {3.0560, 3.0},    {3.1320, 5.0},    {3.1560, 18.0},   {3.2040, 1.2},    {3.2450, 3.0},
           {3.3170, 12.0},   {3.3440, 3.0},    {3.4500, 12.2},   {3.5730, 11.0},   {3.7650, 9.0},
           {4.0450, 6.9}

          });

        return aSolarRadiation;
    }

    std::shared_ptr<CSpectralSampleData> loadSampleData_NFRC_31111()
    {
        auto aMeasurements = CSpectralSampleData::create(
          {{0.3, 0.0, 0.0794, 0.0674},   {0.305, 0.0, 0.0953, 0.0669}, {0.31, 0.0, 0.1577, 0.0663},
           {0.315, 0.0, 0.262, 0.0651},  {0.32, 0.0, 0.3688, 0.0638},  {0.325, 0.0, 0.445, 0.0625},
           {0.33, 0.0, 0.5006, 0.0614},  {0.335, 0.0, 0.5515, 0.0605}, {0.34, 0.0, 0.6044, 0.0603},
           {0.345, 0.0, 0.6534, 0.0612}, {0.35, 0.0, 0.6938, 0.0635},  {0.355, 0.0, 0.7303, 0.0674},
           {0.36, 0.0, 0.7649, 0.0728},  {0.365, 0.0, 0.7881, 0.0801}, {0.37, 0.0, 0.8012, 0.0896},
           {0.375, 0.0, 0.8125, 0.1026}, {0.38, 0.0, 0.825, 0.1206},   {0.385, 0.0, 0.8299, 0.1485},
           {0.39, 0.0, 0.8232, 0.1934},  {0.395, 0.0, 0.8162, 0.2594}, {0.4, 0.0, 0.813, 0.3358},
           {0.405, 0.0, 0.8089, 0.4036}, {0.41, 0.0, 0.8041, 0.4459},  {0.415, 0.0, 0.802, 0.4644},
           {0.42, 0.0, 0.8044, 0.4695},  {0.425, 0.0, 0.8113, 0.4697}, {0.43, 0.0, 0.8192, 0.4689},
           {0.435, 0.0, 0.8266, 0.468},  {0.44, 0.0, 0.8329, 0.4672},  {0.445, 0.0, 0.8382, 0.4665},
           {0.45, 0.0, 0.843, 0.4661},   {0.455, 0.0, 0.8466, 0.4652}, {0.46, 0.0, 0.8498, 0.4648},
           {0.465, 0.0, 0.852, 0.464},   {0.47, 0.0, 0.8545, 0.4632},  {0.475, 0.0, 0.8563, 0.4622},
           {0.48, 0.0, 0.8583, 0.4617},  {0.485, 0.0, 0.8597, 0.4607}, {0.49, 0.0, 0.8608, 0.4599},
           {0.495, 0.0, 0.8628, 0.4599}, {0.5, 0.0, 0.8644, 0.461},    {0.505, 0.0, 0.8658, 0.4642},
           {0.51, 0.0, 0.8668, 0.4694},  {0.515, 0.0, 0.8677, 0.475},  {0.52, 0.0, 0.8685, 0.4789},
           {0.525, 0.0, 0.8695, 0.481},  {0.53, 0.0, 0.8703, 0.4823},  {0.535, 0.0, 0.8713, 0.4825},
           {0.54, 0.0, 0.872, 0.4821},   {0.545, 0.0, 0.8718, 0.4814}, {0.55, 0.0, 0.8725, 0.4808},
           {0.555, 0.0, 0.8728, 0.4801}, {0.56, 0.0, 0.8728, 0.4791},  {0.565, 0.0, 0.8731, 0.4779},
           {0.57, 0.0, 0.8733, 0.4769},  {0.575, 0.0, 0.8726, 0.4758}, {0.58, 0.0, 0.8737, 0.4756},
           {0.585, 0.0, 0.8733, 0.4749}, {0.59, 0.0, 0.8731, 0.4743},  {0.595, 0.0, 0.8721, 0.4733},
           {0.6, 0.0, 0.8716, 0.4722},   {0.605, 0.0, 0.8715, 0.4715}, {0.61, 0.0, 0.8711, 0.4712},
           {0.615, 0.0, 0.87, 0.4707},   {0.62, 0.0, 0.8696, 0.4706},  {0.625, 0.0, 0.8689, 0.4705},
           {0.63, 0.0, 0.8684, 0.4708},  {0.635, 0.0, 0.8676, 0.4708}, {0.64, 0.0, 0.8669, 0.4708},
           {0.645, 0.0, 0.8656, 0.4704}, {0.65, 0.0, 0.8648, 0.4705},  {0.655, 0.0, 0.8643, 0.4704},
           {0.66, 0.0, 0.8633, 0.4697},  {0.665, 0.0, 0.8626, 0.4688}, {0.67, 0.0, 0.8618, 0.4675},
           {0.675, 0.0, 0.8607, 0.466},  {0.68, 0.0, 0.8597, 0.4651},  {0.685, 0.0, 0.8591, 0.4646},
           {0.69, 0.0, 0.8581, 0.4641},  {0.695, 0.0, 0.8568, 0.4636}, {0.7, 0.0, 0.8558, 0.4632},
           {0.705, 0.0, 0.8552, 0.4633}, {0.71, 0.0, 0.8543, 0.4638},  {0.715, 0.0, 0.8528, 0.4643},
           {0.72, 0.0, 0.8522, 0.4652},  {0.725, 0.0, 0.8506, 0.4654}, {0.73, 0.0, 0.8496, 0.4657},
           {0.735, 0.0, 0.8485, 0.4659}, {0.74, 0.0, 0.8471, 0.4662},  {0.745, 0.0, 0.8456, 0.4662},
           {0.75, 0.0, 0.8449, 0.4663},  {0.755, 0.0, 0.8431, 0.4658}, {0.76, 0.0, 0.8422, 0.4655},
           {0.765, 0.0, 0.8408, 0.4653}, {0.77, 0.0, 0.8397, 0.4648},  {0.775, 0.0, 0.8382, 0.4642},
           {0.78, 0.0, 0.8366, 0.4632},  {0.785, 0.0, 0.8354, 0.4628}, {0.79, 0.0, 0.8341, 0.4619},
           {0.795, 0.0, 0.8319, 0.461},  {0.8, 0.0, 0.8311, 0.4603},   {0.805, 0.0, 0.8293, 0.4602},
           {0.81, 0.0, 0.8287, 0.4596},  {0.815, 0.0, 0.828, 0.4591},  {0.82, 0.0, 0.8271, 0.4591},
           {0.825, 0.0, 0.8261, 0.4589}, {0.83, 0.0, 0.8259, 0.4593},  {0.835, 0.0, 0.8245, 0.4593},
           {0.84, 0.0, 0.8244, 0.4594},  {0.845, 0.0, 0.8264, 0.4617}, {0.85, 0.0, 0.8255, 0.4611},
           {0.855, 0.0, 0.8247, 0.4612}, {0.86, 0.0, 0.8227, 0.4592},  {0.865, 0.0, 0.828, 0.4651},
           {0.87, 0.0, 0.8261, 0.4674},  {0.875, 0.0, 0.8273, 0.4661}, {0.88, 0.0, 0.8287, 0.469},
           {0.885, 0.0, 0.8295, 0.4682}, {0.89, 0.0, 0.832, 0.4723},   {0.895, 0.0, 0.8291, 0.4711},
           {0.9, 0.0, 0.8326, 0.4727},   {0.905, 0.0, 0.8323, 0.4743}, {0.91, 0.0, 0.8338, 0.4756},
           {0.915, 0.0, 0.8356, 0.4791}, {0.92, 0.0, 0.8391, 0.4816},  {0.925, 0.0, 0.8391, 0.4806},
           {0.93, 0.0, 0.8414, 0.4841},  {0.935, 0.0, 0.8423, 0.4858}, {0.94, 0.0, 0.8435, 0.4877},
           {0.945, 0.0, 0.8442, 0.4882}, {0.95, 0.0, 0.8456, 0.4908},  {0.955, 0.0, 0.8478, 0.4923},
           {0.96, 0.0, 0.8483, 0.4948},  {0.965, 0.0, 0.8499, 0.4967}, {0.97, 0.0, 0.8505, 0.4993},
           {0.975, 0.0, 0.8534, 0.5},    {0.98, 0.0, 0.8528, 0.5004},  {0.985, 0.0, 0.8542, 0.5008},
           {0.99, 0.0, 0.8557, 0.5038},  {0.995, 0.0, 0.8546, 0.5051}, {1.0, 0.0, 0.8569, 0.5056},
           {1.005, 0.0, 0.8574, 0.5078}, {1.01, 0.0, 0.8586, 0.5083},  {1.015, 0.0, 0.8601, 0.5101},
           {1.02, 0.0, 0.8609, 0.5117},  {1.025, 0.0, 0.8614, 0.5125}, {1.03, 0.0, 0.8614, 0.5129},
           {1.035, 0.0, 0.8637, 0.516},  {1.04, 0.0, 0.864, 0.5171},   {1.045, 0.0, 0.8646, 0.5177},
           {1.05, 0.0, 0.8659, 0.5203},  {1.055, 0.0, 0.8668, 0.5206}, {1.06, 0.0, 0.867, 0.5211},
           {1.065, 0.0, 0.8678, 0.5233}, {1.07, 0.0, 0.8686, 0.524},   {1.075, 0.0, 0.8696, 0.5254},
           {1.08, 0.0, 0.8699, 0.5253},  {1.085, 0.0, 0.8703, 0.5273}, {1.09, 0.0, 0.8707, 0.5281},
           {1.095, 0.0, 0.8711, 0.5287}, {1.1, 0.0, 0.8724, 0.5299},   {1.105, 0.0, 0.8716, 0.5314},
           {1.11, 0.0, 0.8713, 0.5316},  {1.115, 0.0, 0.8697, 0.5308}, {1.12, 0.0, 0.8669, 0.5307},
           {1.125, 0.0, 0.8637, 0.5297}, {1.13, 0.0, 0.8627, 0.5297},  {1.135, 0.0, 0.8625, 0.5298},
           {1.14, 0.0, 0.8655, 0.5319},  {1.145, 0.0, 0.8677, 0.5332}, {1.15, 0.0, 0.8695, 0.5344},
           {1.155, 0.0, 0.8695, 0.5342}, {1.16, 0.0, 0.8697, 0.5337},  {1.165, 0.0, 0.869, 0.5327},
           {1.17, 0.0, 0.8682, 0.5316},  {1.175, 0.0, 0.8689, 0.5317}, {1.18, 0.0, 0.8682, 0.5322},
           {1.185, 0.0, 0.8702, 0.534},  {1.19, 0.0, 0.8722, 0.5354},  {1.195, 0.0, 0.8737, 0.5357},
           {1.2, 0.0, 0.8755, 0.5377},   {1.205, 0.0, 0.877, 0.5404},  {1.21, 0.0, 0.8777, 0.5429},
           {1.215, 0.0, 0.8788, 0.5454}, {1.22, 0.0, 0.8794, 0.5481},  {1.225, 0.0, 0.8797, 0.55},
           {1.23, 0.0, 0.8803, 0.5523},  {1.235, 0.0, 0.8812, 0.5547}, {1.24, 0.0, 0.8813, 0.5555},
           {1.245, 0.0, 0.882, 0.5568},  {1.25, 0.0, 0.8831, 0.5581},  {1.255, 0.0, 0.8836, 0.559},
           {1.26, 0.0, 0.8832, 0.5596},  {1.265, 0.0, 0.884, 0.5616},  {1.27, 0.0, 0.8839, 0.5626},
           {1.275, 0.0, 0.8845, 0.5635}, {1.28, 0.0, 0.8845, 0.5641},  {1.285, 0.0, 0.884, 0.5653},
           {1.29, 0.0, 0.8842, 0.5655},  {1.295, 0.0, 0.8844, 0.5666}, {1.3, 0.0, 0.8846, 0.5678},
           {1.305, 0.0, 0.8842, 0.568},  {1.31, 0.0, 0.8842, 0.5689},  {1.315, 0.0, 0.8844, 0.5696},
           {1.32, 0.0, 0.8842, 0.57},    {1.325, 0.0, 0.8835, 0.5706}, {1.33, 0.0, 0.8823, 0.5716},
           {1.335, 0.0, 0.8825, 0.5715}, {1.34, 0.0, 0.882, 0.5716},   {1.345, 0.0, 0.8804, 0.5712},
           {1.35, 0.0, 0.878, 0.5715},   {1.355, 0.0, 0.8762, 0.5705}, {1.36, 0.0, 0.8745, 0.5704},
           {1.365, 0.0, 0.873, 0.57},    {1.37, 0.0, 0.8722, 0.5693},  {1.375, 0.0, 0.8712, 0.5686},
           {1.38, 0.0, 0.8695, 0.567},   {1.385, 0.0, 0.8682, 0.5671}, {1.39, 0.0, 0.8663, 0.5661},
           {1.395, 0.0, 0.8635, 0.5656}, {1.4, 0.0, 0.8608, 0.5642},   {1.405, 0.0, 0.8604, 0.5647},
           {1.41, 0.0, 0.8587, 0.5639},  {1.415, 0.0, 0.8575, 0.5632}, {1.42, 0.0, 0.8605, 0.5643},
           {1.425, 0.0, 0.8632, 0.5655}, {1.43, 0.0, 0.8666, 0.5673},  {1.435, 0.0, 0.8695, 0.5706},
           {1.44, 0.0, 0.871, 0.5736},   {1.445, 0.0, 0.8733, 0.5767}, {1.45, 0.0, 0.8749, 0.5796},
           {1.455, 0.0, 0.876, 0.581},   {1.46, 0.0, 0.8771, 0.5836},  {1.465, 0.0, 0.8779, 0.5852},
           {1.47, 0.0, 0.8792, 0.5865},  {1.475, 0.0, 0.8796, 0.5875}, {1.48, 0.0, 0.8807, 0.5897},
           {1.485, 0.0, 0.8821, 0.5915}, {1.49, 0.0, 0.8831, 0.5929},  {1.495, 0.0, 0.8835, 0.5943},
           {1.5, 0.0, 0.8838, 0.5952},   {1.505, 0.0, 0.8836, 0.5958}, {1.51, 0.0, 0.8848, 0.5976},
           {1.515, 0.0, 0.8862, 0.599},  {1.52, 0.0, 0.8871, 0.6001},  {1.525, 0.0, 0.8866, 0.6003},
           {1.53, 0.0, 0.887, 0.6019},   {1.535, 0.0, 0.8872, 0.603},  {1.54, 0.0, 0.8865, 0.6034},
           {1.545, 0.0, 0.8869, 0.6046}, {1.55, 0.0, 0.8872, 0.6055},  {1.555, 0.0, 0.8877, 0.6067},
           {1.56, 0.0, 0.8887, 0.6076},  {1.565, 0.0, 0.8888, 0.6076}, {1.57, 0.0, 0.8889, 0.6093},
           {1.575, 0.0, 0.8889, 0.6109}, {1.58, 0.0, 0.8884, 0.6112},  {1.585, 0.0, 0.8873, 0.6117},
           {1.59, 0.0, 0.888, 0.6116},   {1.595, 0.0, 0.8872, 0.6116}, {1.6, 0.0, 0.8856, 0.6108},
           {1.605, 0.0, 0.8842, 0.6105}, {1.61, 0.0, 0.8813, 0.6097},  {1.615, 0.0, 0.8776, 0.6081},
           {1.62, 0.0, 0.8739, 0.6059},  {1.625, 0.0, 0.8703, 0.6045}, {1.63, 0.0, 0.8647, 0.6017},
           {1.635, 0.0, 0.8534, 0.5946}, {1.64, 0.0, 0.834, 0.5827},   {1.645, 0.0, 0.8061, 0.5656},
           {1.65, 0.0, 0.7724, 0.545},   {1.655, 0.0, 0.7401, 0.5253}, {1.66, 0.0, 0.7178, 0.5112},
           {1.665, 0.0, 0.7155, 0.5074}, {1.67, 0.0, 0.73, 0.5152},    {1.675, 0.0, 0.7506, 0.5258},
           {1.68, 0.0, 0.7697, 0.5361},  {1.685, 0.0, 0.7853, 0.5414}, {1.69, 0.0, 0.7992, 0.5431},
           {1.695, 0.0, 0.8083, 0.5403}, {1.7, 0.0, 0.81, 0.5277},     {1.705, 0.0, 0.8094, 0.5091},
           {1.71, 0.0, 0.8092, 0.4933},  {1.715, 0.0, 0.8097, 0.4835}, {1.72, 0.0, 0.81, 0.4831},
           {1.725, 0.0, 0.8114, 0.4936}, {1.73, 0.0, 0.8144, 0.5065},  {1.735, 0.0, 0.8201, 0.5199},
           {1.74, 0.0, 0.8286, 0.5294},  {1.745, 0.0, 0.8344, 0.5313}, {1.75, 0.0, 0.839, 0.5354},
           {1.755, 0.0, 0.8423, 0.5404}, {1.76, 0.0, 0.8427, 0.5451},  {1.765, 0.0, 0.8443, 0.5515},
           {1.77, 0.0, 0.8471, 0.5576},  {1.775, 0.0, 0.8493, 0.5626}, {1.78, 0.0, 0.8519, 0.5672},
           {1.785, 0.0, 0.8537, 0.5716}, {1.79, 0.0, 0.853, 0.5746},   {1.795, 0.0, 0.8518, 0.5779},
           {1.8, 0.0, 0.8511, 0.5813},   {1.805, 0.0, 0.8507, 0.5838}, {1.81, 0.0, 0.8516, 0.5861},
           {1.815, 0.0, 0.8528, 0.589},  {1.82, 0.0, 0.8541, 0.5898},  {1.825, 0.0, 0.8543, 0.5901},
           {1.83, 0.0, 0.8555, 0.5909},  {1.835, 0.0, 0.8582, 0.5936}, {1.84, 0.0, 0.8624, 0.5994},
           {1.845, 0.0, 0.865, 0.6022},  {1.85, 0.0, 0.8663, 0.603},   {1.855, 0.0, 0.8672, 0.606},
           {1.86, 0.0, 0.8669, 0.608},   {1.865, 0.0, 0.8646, 0.6083}, {1.87, 0.0, 0.864, 0.6102},
           {1.875, 0.0, 0.8655, 0.6142}, {1.88, 0.0, 0.8644, 0.6141},  {1.885, 0.0, 0.8585, 0.614},
           {1.89, 0.0, 0.8522, 0.6111},  {1.895, 0.0, 0.8435, 0.6054}, {1.9, 0.0, 0.829, 0.5963},
           {1.905, 0.0, 0.816, 0.5899},  {1.91, 0.0, 0.8125, 0.5888},  {1.915, 0.0, 0.8151, 0.5923},
           {1.92, 0.0, 0.8224, 0.6005},  {1.925, 0.0, 0.8297, 0.6071}, {1.93, 0.0, 0.8348, 0.6096},
           {1.935, 0.0, 0.8385, 0.6116}, {1.94, 0.0, 0.8368, 0.6117},  {1.945, 0.0, 0.8339, 0.6104},
           {1.95, 0.0, 0.8318, 0.6104},  {1.955, 0.0, 0.8337, 0.6117}, {1.96, 0.0, 0.8378, 0.6156},
           {1.965, 0.0, 0.8449, 0.6204}, {1.97, 0.0, 0.8521, 0.6244},  {1.975, 0.0, 0.8578, 0.6303},
           {1.98, 0.0, 0.8626, 0.6344},  {1.985, 0.0, 0.8658, 0.6364}, {1.99, 0.0, 0.8678, 0.6382},
           {1.995, 0.0, 0.8685, 0.641},  {2.0, 0.0, 0.8711, 0.6428},   {2.005, 0.0, 0.8756, 0.6449},
           {2.01, 0.0, 0.8776, 0.6476},  {2.015, 0.0, 0.8782, 0.6488}, {2.02, 0.0, 0.8804, 0.6512},
           {2.025, 0.0, 0.8803, 0.6515}, {2.03, 0.0, 0.8791, 0.6516},  {2.035, 0.0, 0.8803, 0.6544},
           {2.04, 0.0, 0.8797, 0.6561},  {2.045, 0.0, 0.8795, 0.657},  {2.05, 0.0, 0.8809, 0.6583},
           {2.055, 0.0, 0.8818, 0.6587}, {2.06, 0.0, 0.8806, 0.6578},  {2.065, 0.0, 0.8777, 0.6592},
           {2.07, 0.0, 0.8681, 0.6553},  {2.075, 0.0, 0.8616, 0.6476}, {2.08, 0.0, 0.8556, 0.6435},
           {2.085, 0.0, 0.8503, 0.646},  {2.09, 0.0, 0.8474, 0.6464},  {2.095, 0.0, 0.851, 0.6485},
           {2.1, 0.0, 0.8471, 0.6444},   {2.105, 0.0, 0.8367, 0.6371}, {2.11, 0.0, 0.8184, 0.6242},
           {2.115, 0.0, 0.7964, 0.6108}, {2.12, 0.0, 0.7694, 0.596},   {2.125, 0.0, 0.7409, 0.572},
           {2.13, 0.0, 0.7156, 0.5579},  {2.135, 0.0, 0.7054, 0.5534}, {2.14, 0.0, 0.7276, 0.5668},
           {2.145, 0.0, 0.7438, 0.5783}, {2.15, 0.0, 0.751, 0.5758},   {2.155, 0.0, 0.745, 0.5845},
           {2.16, 0.0, 0.7412, 0.5841},  {2.165, 0.0, 0.7529, 0.5928}, {2.17, 0.0, 0.7705, 0.6003},
           {2.175, 0.0, 0.7749, 0.5999}, {2.18, 0.0, 0.7671, 0.5989},  {2.185, 0.0, 0.7612, 0.5951},
           {2.19, 0.0, 0.7661, 0.6027},  {2.195, 0.0, 0.7782, 0.6117}, {2.2, 0.0, 0.7861, 0.6156},
           {2.205, 0.0, 0.7872, 0.6168}, {2.21, 0.0, 0.7804, 0.6127},  {2.215, 0.0, 0.7654, 0.5993},
           {2.22, 0.0, 0.742, 0.5777},   {2.225, 0.0, 0.7072, 0.554},  {2.23, 0.0, 0.6783, 0.5291},
           {2.235, 0.0, 0.6453, 0.5009}, {2.24, 0.0, 0.6083, 0.4756},  {2.245, 0.0, 0.5777, 0.4527},
           {2.25, 0.0, 0.5574, 0.4348},  {2.255, 0.0, 0.5448, 0.4219}, {2.26, 0.0, 0.5386, 0.4137},
           {2.265, 0.0, 0.554, 0.4192},  {2.27, 0.0, 0.5583, 0.415},   {2.275, 0.0, 0.5736, 0.4229},
           {2.28, 0.0, 0.5922, 0.4156},  {2.285, 0.0, 0.5948, 0.4011}, {2.29, 0.0, 0.5951, 0.3768},
           {2.295, 0.0, 0.595, 0.3619},  {2.3, 0.0, 0.5922, 0.343},    {2.305, 0.0, 0.5831, 0.3256},
           {2.31, 0.0, 0.5791, 0.3252},  {2.315, 0.0, 0.5816, 0.3253}, {2.32, 0.0, 0.5773, 0.3314},
           {2.325, 0.0, 0.5645, 0.34},   {2.33, 0.0, 0.5516, 0.3379},  {2.335, 0.0, 0.5437, 0.33},
           {2.34, 0.0, 0.5544, 0.3443},  {2.345, 0.0, 0.5608, 0.354},  {2.35, 0.0, 0.5702, 0.3631},
           {2.355, 0.0, 0.5737, 0.3532}, {2.36, 0.0, 0.5822, 0.3593},  {2.365, 0.0, 0.5715, 0.3483},
           {2.37, 0.0, 0.5866, 0.3617},  {2.375, 0.0, 0.593, 0.3495},  {2.38, 0.0, 0.595, 0.3652},
           {2.385, 0.0, 0.5984, 0.3664}, {2.39, 0.0, 0.588, 0.3649},   {2.395, 0.0, 0.591, 0.3742},
           {2.4, 0.0, 0.5903, 0.3775},   {2.405, 0.0, 0.6053, 0.3805}, {2.41, 0.0, 0.6094, 0.3892},
           {2.415, 0.0, 0.5974, 0.3893}, {2.42, 0.0, 0.5962, 0.4049},  {2.425, 0.0, 0.5776, 0.3955},
           {2.43, 0.0, 0.5771, 0.3996},  {2.435, 0.0, 0.5455, 0.3828}, {2.44, 0.0, 0.5008, 0.3389},
           {2.445, 0.0, 0.48, 0.3346},   {2.45, 0.0, 0.4861, 0.3431},  {2.455, 0.0, 0.5135, 0.3627},
           {2.46, 0.0, 0.5599, 0.3991},  {2.465, 0.0, 0.5829, 0.4012}, {2.47, 0.0, 0.5987, 0.429},
           {2.475, 0.0, 0.6275, 0.4274}, {2.48, 0.0, 0.6452, 0.4492},  {2.485, 0.0, 0.6492, 0.4508},
           {2.49, 0.0, 0.6502, 0.4533},  {2.495, 0.0, 0.6618, 0.461},  {2.5, 0.0, 0.6568, 0.48}});

        return aMeasurements;
    }

    std::shared_ptr<CSpectralSampleData> loadSampleData_NFRC_102()
    {
        auto aMeasurements_102 = CSpectralSampleData::create(
          {{0.300, 0.0020, 0.0470, 0.0480}, {0.305, 0.0030, 0.0470, 0.0480},
           {0.310, 0.0090, 0.0470, 0.0480}, {0.315, 0.0350, 0.0470, 0.0480},
           {0.320, 0.1000, 0.0470, 0.0480}, {0.325, 0.2180, 0.0490, 0.0500},
           {0.330, 0.3560, 0.0530, 0.0540}, {0.335, 0.4980, 0.0600, 0.0610},
           {0.340, 0.6160, 0.0670, 0.0670}, {0.345, 0.7090, 0.0730, 0.0740},
           {0.350, 0.7740, 0.0780, 0.0790}, {0.355, 0.8180, 0.0820, 0.0820},
           {0.360, 0.8470, 0.0840, 0.0840}, {0.365, 0.8630, 0.0850, 0.0850},
           {0.370, 0.8690, 0.0850, 0.0860}, {0.375, 0.8610, 0.0850, 0.0850},
           {0.380, 0.8560, 0.0840, 0.0840}, {0.385, 0.8660, 0.0850, 0.0850},
           {0.390, 0.8810, 0.0860, 0.0860}, {0.395, 0.8890, 0.0860, 0.0860},
           {0.400, 0.8930, 0.0860, 0.0860}, {0.410, 0.8930, 0.0860, 0.0860},
           {0.420, 0.8920, 0.0860, 0.0860}, {0.430, 0.8920, 0.0850, 0.0850},
           {0.440, 0.8920, 0.0850, 0.0850}, {0.450, 0.8960, 0.0850, 0.0850},
           {0.460, 0.9000, 0.0850, 0.0850}, {0.470, 0.9020, 0.0840, 0.0840},
           {0.480, 0.9030, 0.0840, 0.0840}, {0.490, 0.9040, 0.0850, 0.0850},
           {0.500, 0.9050, 0.0840, 0.0840}, {0.510, 0.9050, 0.0840, 0.0840},
           {0.520, 0.9050, 0.0840, 0.0840}, {0.530, 0.9040, 0.0840, 0.0840},
           {0.540, 0.9040, 0.0830, 0.0830}, {0.550, 0.9030, 0.0830, 0.0830},
           {0.560, 0.9020, 0.0830, 0.0830}, {0.570, 0.9000, 0.0820, 0.0820},
           {0.580, 0.8980, 0.0820, 0.0820}, {0.590, 0.8960, 0.0810, 0.0810},
           {0.600, 0.8930, 0.0810, 0.0810}, {0.610, 0.8900, 0.0810, 0.0810},
           {0.620, 0.8860, 0.0800, 0.0800}, {0.630, 0.8830, 0.0800, 0.0800},
           {0.640, 0.8790, 0.0790, 0.0790}, {0.650, 0.8750, 0.0790, 0.0790},
           {0.660, 0.8720, 0.0790, 0.0790}, {0.670, 0.8680, 0.0780, 0.0780},
           {0.680, 0.8630, 0.0780, 0.0780}, {0.690, 0.8590, 0.0770, 0.0770},
           {0.700, 0.8540, 0.0760, 0.0770}, {0.710, 0.8500, 0.0760, 0.0760},
           {0.720, 0.8450, 0.0750, 0.0760}, {0.730, 0.8400, 0.0750, 0.0750},
           {0.740, 0.8350, 0.0750, 0.0750}, {0.750, 0.8310, 0.0740, 0.0740},
           {0.760, 0.8260, 0.0740, 0.0740}, {0.770, 0.8210, 0.0740, 0.0740},
           {0.780, 0.8160, 0.0730, 0.0730}, {0.790, 0.8120, 0.0730, 0.0730},
           {0.800, 0.8080, 0.0720, 0.0720}, {0.810, 0.8030, 0.0720, 0.0720},
           {0.820, 0.8000, 0.0720, 0.0720}, {0.830, 0.7960, 0.0710, 0.0710},
           {0.840, 0.7930, 0.0700, 0.0710}, {0.850, 0.7880, 0.0700, 0.0710},
           {0.860, 0.7860, 0.0700, 0.0700}, {0.870, 0.7820, 0.0740, 0.0740},
           {0.880, 0.7800, 0.0720, 0.0720}, {0.890, 0.7770, 0.0730, 0.0740},
           {0.900, 0.7760, 0.0720, 0.0720}, {0.910, 0.7730, 0.0720, 0.0720},
           {0.920, 0.7710, 0.0710, 0.0710}, {0.930, 0.7700, 0.0700, 0.0700},
           {0.940, 0.7680, 0.0690, 0.0690}, {0.950, 0.7660, 0.0680, 0.0680},
           {0.960, 0.7660, 0.0670, 0.0680}, {0.970, 0.7640, 0.0680, 0.0680},
           {0.980, 0.7630, 0.0680, 0.0680}, {0.990, 0.7620, 0.0670, 0.0670},
           {1.000, 0.7620, 0.0660, 0.0670}, {1.050, 0.7600, 0.0660, 0.0660},
           {1.100, 0.7590, 0.0660, 0.0660}, {1.150, 0.7610, 0.0660, 0.0660},
           {1.200, 0.7650, 0.0660, 0.0660}, {1.250, 0.7700, 0.0650, 0.0650},
           {1.300, 0.7770, 0.0670, 0.0670}, {1.350, 0.7860, 0.0660, 0.0670},
           {1.400, 0.7950, 0.0670, 0.0680}, {1.450, 0.8080, 0.0670, 0.0670},
           {1.500, 0.8190, 0.0690, 0.0690}, {1.550, 0.8290, 0.0690, 0.0690},
           {1.600, 0.8360, 0.0700, 0.0700}, {1.650, 0.8400, 0.0700, 0.0700},
           {1.700, 0.8420, 0.0690, 0.0700}, {1.750, 0.8420, 0.0690, 0.0700},
           {1.800, 0.8410, 0.0700, 0.0700}, {1.850, 0.8400, 0.0690, 0.0690},
           {1.900, 0.8390, 0.0680, 0.0680}, {1.950, 0.8390, 0.0710, 0.0710},
           {2.000, 0.8390, 0.0690, 0.0690}, {2.050, 0.8400, 0.0680, 0.0680},
           {2.100, 0.8410, 0.0680, 0.0680}, {2.150, 0.8390, 0.0690, 0.0690},
           {2.200, 0.8300, 0.0700, 0.0700}, {2.250, 0.8300, 0.0700, 0.0700},
           {2.300, 0.8320, 0.0690, 0.0690}, {2.350, 0.8320, 0.0690, 0.0700},
           {2.400, 0.8320, 0.0700, 0.0700}, {2.450, 0.8260, 0.0690, 0.0690},
           {2.500, 0.8220, 0.0680, 0.0680}});

        return aMeasurements_102;
    }

protected:
    virtual void SetUp()
    {
        const auto solarRadiation{loadSolarRadiationFile()};
        const auto wl{solarRadiation.getXArray()};

        const auto aBSDF = BSDFHemisphere::create(BSDFBasis::Quarter);

        auto thickness_102 = 3.048e-3;   // [m]

        auto aMaterial_102 = SingleLayerOptics::Material::nBandMaterial(loadSampleData_NFRC_102(),
                                                                        thickness_102,
                                                                        MaterialType::Monolithic,
                                                                        WavelengthRange::Solar);
        aMaterial_102->setBandWavelengths(wl);

        auto Layer_102 = CBSDFLayerMaker::getSpecularLayer(aMaterial_102, aBSDF);

        const auto thickness_31111{0.00023};
        auto aMaterial_31111 =
          SingleLayerOptics::Material::nBandMaterial(loadSampleData_NFRC_31111(),
                                                     thickness_31111,
                                                     MaterialType::Monolithic,
                                                     WavelengthRange::Solar);
        aMaterial_31111->setBandWavelengths(wl);

        // make cell geometry
        const auto x = 0.00169;        // m
        const auto y = 0.00169;        // m
        const auto radius = 0.00058;   // m

        // Perforated layer is created here
        const auto perforated = CBSDFLayerMaker::getCircularPerforatedLayer(
          aMaterial_31111, aBSDF, x, y, thickness_31111, radius);

        m_Layer = CMultiPaneBSDF::create({perforated, Layer_102});

        const CalculationProperties input{loadSolarRadiationFile(),
                                          loadSolarRadiationFile().getXArray()};
        m_Layer->setCalculationProperties(input);
    }

public:
    CMultiPaneBSDF & getLayer() const
    {
        return *m_Layer;
    };
};

TEST_F(MultiPaneBSDF_102_PerforatedCircular_NFRC18000, Test102PerofratedCircular)
{
    SCOPED_TRACE("Begin Test: Specular layer - BSDF.");

    const double minLambda = 0.3;
    const double maxLambda = 2.5;

    CMultiPaneBSDF & aLayer = getLayer();

    double tauDiff = aLayer.DiffDiff(minLambda, maxLambda, Side::Front, PropertySimple::T);
    EXPECT_NEAR(0.21249128345672524, tauDiff, 1e-6);

    double rhoDiff = aLayer.DiffDiff(minLambda, maxLambda, Side::Front, PropertySimple::R);
    EXPECT_NEAR(0.63539424553486545, rhoDiff, 1e-6);

    double absDiff1 = aLayer.AbsDiff(minLambda, maxLambda, Side::Front, 1);
    EXPECT_NEAR(0.1254558414581273, absDiff1, 1e-6);

    double absDiff2 = aLayer.AbsDiff(minLambda, maxLambda, Side::Front, 2);
    EXPECT_NEAR(0.02665862955027334, absDiff2, 1e-6);

    double theta = 0;
    double phi = 0;

    double tauHem = aLayer.DirHem(minLambda, maxLambda, Side::Front, PropertySimple::T, theta, phi);
    EXPECT_NEAR(0.31523300551099365, tauHem, 1e-6);

    double tauDir = aLayer.DirDir(minLambda, maxLambda, Side::Front, PropertySimple::T, theta, phi);
    EXPECT_NEAR(0.30872668074190435, tauDir, 1e-6);

    double rhoHem = aLayer.DirHem(minLambda, maxLambda, Side::Front, PropertySimple::R, theta, phi);
    EXPECT_NEAR(0.54336046822785034, rhoHem, 1e-6);

    double rhoDir = aLayer.DirDir(minLambda, maxLambda, Side::Front, PropertySimple::R, theta, phi);
    EXPECT_NEAR(0.023282859797643073, rhoDir, 1e-6);

    double abs1 = aLayer.Abs(minLambda, maxLambda, Side::Front, 1, theta, phi);
    EXPECT_NEAR(0.10672566908222603, abs1, 1e-6);

    double abs2 = aLayer.Abs(minLambda, maxLambda, Side::Front, 2, theta, phi);
    EXPECT_NEAR(0.034680857178930609, abs2, 1e-6);
}
