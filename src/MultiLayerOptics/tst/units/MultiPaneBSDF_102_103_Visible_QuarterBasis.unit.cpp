#include <memory>
#include <gtest/gtest.h>

#include "WCESpectralAveraging.hpp"
#include "WCEMultiLayerOptics.hpp"
#include "WCECommon.hpp"

using MultiLayerOptics::CMultiPaneBSDF;
using FenestrationCommon::CSeries;
using SpectralAveraging::CSpectralSampleData;

// Example on how to create multilayer BSDF from specular layers only

class MultiPaneBSDF_102_103_Visible_QuarterBasis : public testing::Test
{
private:
    std::unique_ptr<CMultiPaneBSDF> m_Layer;

    static CSeries loadSolarRadiationFile()
    {
        // Full CIE Illuminant D651 nm ssp table (used for PHOTOPIC properties)
        CSeries solarRadiation(
          {{0.300, 0.0341000}, {0.301, 0.3601400}, {0.302, 0.6861800}, {0.303, 1.0122200},
           {0.304, 1.3382600}, {0.305, 1.6643000}, {0.306, 1.9903400}, {0.307, 2.3163800},
           {0.308, 2.6424200}, {0.309, 2.9684600}, {0.310, 3.2945000}, {0.311, 4.9886500},
           {0.312, 6.6828000}, {0.313, 8.3769500}, {0.314, 10.071100}, {0.315, 11.765200},
           {0.316, 13.459400}, {0.317, 15.153500}, {0.318, 16.847700}, {0.319, 18.541800},
           {0.320, 20.236000}, {0.321, 21.917700}, {0.322, 23.599500}, {0.323, 25.281200},
           {0.324, 26.963000}, {0.325, 28.644700}, {0.326, 30.326500}, {0.327, 32.008200},
           {0.328, 33.690000}, {0.329, 35.371700}, {0.330, 37.053500}, {0.331, 37.343000},
           {0.332, 37.632600}, {0.333, 37.922100}, {0.334, 38.211600}, {0.335, 38.501100},
           {0.336, 38.790700}, {0.337, 39.080200}, {0.338, 39.369700}, {0.339, 39.659300},
           {0.340, 39.948800}, {0.341, 40.445100}, {0.342, 40.941400}, {0.343, 41.437700},
           {0.344, 41.934000}, {0.345, 42.430200}, {0.346, 42.926500}, {0.347, 43.422800},
           {0.348, 43.919100}, {0.349, 44.415400}, {0.350, 44.911700}, {0.351, 45.084400},
           {0.352, 45.257000}, {0.353, 45.429700}, {0.354, 45.602300}, {0.355, 45.775000},
           {0.356, 45.947700}, {0.357, 46.120300}, {0.358, 46.293000}, {0.359, 46.465600},
           {0.360, 46.638300}, {0.361, 47.183400}, {0.362, 47.728500}, {0.363, 48.273500},
           {0.364, 48.818600}, {0.365, 49.363700}, {0.366, 49.908800}, {0.367, 50.453900},
           {0.368, 50.998900}, {0.369, 51.544000}, {0.370, 52.089100}, {0.371, 51.877700},
           {0.372, 51.666400}, {0.373, 51.455000}, {0.374, 51.243700}, {0.375, 51.032300},
           {0.376, 50.820900}, {0.377, 50.609600}, {0.378, 50.398200}, {0.379, 50.186900},
           {0.380, 49.975500}, {0.381, 50.442800}, {0.382, 50.910000}, {0.383, 51.377300},
           {0.384, 51.844600}, {0.385, 52.311800}, {0.386, 52.779100}, {0.387, 53.246400},
           {0.388, 53.713700}, {0.389, 54.180900}, {0.390, 54.648200}, {0.391, 57.458900},
           {0.392, 60.269500}, {0.393, 63.080200}, {0.394, 65.890900}, {0.395, 68.701500},
           {0.396, 71.512200}, {0.397, 74.322900}, {0.398, 77.133600}, {0.399, 79.944200},
           {0.400, 82.754900}, {0.401, 83.628000}, {0.402, 84.501100}, {0.403, 85.374200},
           {0.404, 86.247300}, {0.405, 87.120400}, {0.406, 87.993600}, {0.407, 88.866700},
           {0.408, 89.739800}, {0.409, 90.612900}, {0.410, 91.486000}, {0.411, 91.680600},
           {0.412, 91.875200}, {0.413, 92.069700}, {0.414, 92.264300}, {0.415, 92.458900},
           {0.416, 92.653500}, {0.417, 92.848100}, {0.418, 93.042600}, {0.419, 93.237200},
           {0.420, 93.431800}, {0.421, 92.756800}, {0.422, 92.081900}, {0.423, 91.406900},
           {0.424, 90.732000}, {0.425, 90.057000}, {0.426, 89.382100}, {0.427, 88.707100},
           {0.428, 88.032200}, {0.429, 87.357200}, {0.430, 86.682300}, {0.431, 88.500600},
           {0.432, 90.318800}, {0.433, 92.137100}, {0.434, 93.955400}, {0.435, 95.773600},
           {0.436, 97.591900}, {0.437, 99.410200}, {0.438, 101.22800}, {0.439, 103.04700},
           {0.440, 104.86500}, {0.441, 106.07900}, {0.442, 107.29400}, {0.443, 108.50800},
           {0.444, 109.72200}, {0.445, 110.93600}, {0.446, 112.15100}, {0.447, 113.36500},
           {0.448, 114.57900}, {0.449, 115.79400}, {0.450, 117.00800}, {0.451, 117.08800},
           {0.452, 117.16900}, {0.453, 117.24900}, {0.454, 117.33000}, {0.455, 117.41000},
           {0.456, 117.49000}, {0.457, 117.57100}, {0.458, 117.65100}, {0.459, 117.73200},
           {0.460, 117.81200}, {0.461, 117.51700}, {0.462, 117.22200}, {0.463, 116.92700},
           {0.464, 116.63200}, {0.465, 116.33600}, {0.466, 116.04100}, {0.467, 115.74600},
           {0.468, 115.45100}, {0.469, 115.15600}, {0.470, 114.86100}, {0.471, 114.96700},
           {0.472, 115.07300}, {0.473, 115.18000}, {0.474, 115.28600}, {0.475, 115.39200},
           {0.476, 115.49800}, {0.477, 115.60400}, {0.478, 115.71100}, {0.479, 115.81700},
           {0.480, 115.92300}, {0.481, 115.21200}, {0.482, 114.50100}, {0.483, 113.78900},
           {0.484, 113.07800}, {0.485, 112.36700}, {0.486, 111.65600}, {0.487, 110.94500},
           {0.488, 110.23300}, {0.489, 109.52200}, {0.490, 108.81100}, {0.491, 108.86500},
           {0.492, 108.92000}, {0.493, 108.97400}, {0.494, 109.02800}, {0.495, 109.08200},
           {0.496, 109.13700}, {0.497, 109.19100}, {0.498, 109.24500}, {0.499, 109.30000},
           {0.500, 109.35400}, {0.501, 109.19900}, {0.502, 109.04400}, {0.503, 108.88800},
           {0.504, 108.73300}, {0.505, 108.57800}, {0.506, 108.42300}, {0.507, 108.26800},
           {0.508, 108.11200}, {0.509, 107.95700}, {0.510, 107.80200}, {0.511, 107.50100},
           {0.512, 107.20000}, {0.513, 106.89800}, {0.514, 106.59700}, {0.515, 106.29600},
           {0.516, 105.99500}, {0.517, 105.69400}, {0.518, 105.39200}, {0.519, 105.09100},
           {0.520, 104.79000}, {0.521, 105.08000}, {0.522, 105.37000}, {0.523, 105.66000},
           {0.524, 105.95000}, {0.525, 106.23900}, {0.526, 106.52900}, {0.527, 106.81900},
           {0.528, 107.10900}, {0.529, 107.39900}, {0.530, 107.68900}, {0.531, 107.36100},
           {0.532, 107.03200}, {0.533, 106.70400}, {0.534, 106.37500}, {0.535, 106.04700},
           {0.536, 105.71900}, {0.537, 105.39000}, {0.538, 105.06200}, {0.539, 104.73300},
           {0.540, 104.40500}, {0.541, 104.36900}, {0.542, 104.33300}, {0.543, 104.29700},
           {0.544, 104.26100}, {0.545, 104.22500}, {0.546, 104.19000}, {0.547, 104.15400},
           {0.548, 104.11800}, {0.549, 104.08200}, {0.550, 104.04600}, {0.551, 103.64100},
           {0.552, 103.23700}, {0.553, 102.83200}, {0.554, 102.42800}, {0.555, 102.02300},
           {0.556, 101.61800}, {0.557, 101.21400}, {0.558, 100.80900}, {0.559, 100.40500},
           {0.560, 100.00000}, {0.561, 99.633400}, {0.562, 99.266800}, {0.563, 98.900300},
           {0.564, 98.533700}, {0.565, 98.167100}, {0.566, 97.800500}, {0.567, 97.433900},
           {0.568, 97.067400}, {0.569, 96.700800}, {0.570, 96.334200}, {0.571, 96.279600},
           {0.572, 96.225000}, {0.573, 96.170300}, {0.574, 96.115700}, {0.575, 96.061100},
           {0.576, 96.006500}, {0.577, 95.951900}, {0.578, 95.897200}, {0.579, 95.842600},
           {0.580, 95.788000}, {0.581, 95.077800}, {0.582, 94.367500}, {0.583, 93.657300},
           {0.584, 92.947000}, {0.585, 92.236800}, {0.586, 91.526600}, {0.587, 90.816300},
           {0.588, 90.106100}, {0.589, 89.395800}, {0.590, 88.685600}, {0.591, 88.817700},
           {0.592, 88.949700}, {0.593, 89.081800}, {0.594, 89.213800}, {0.595, 89.345900},
           {0.596, 89.478000}, {0.597, 89.610000}, {0.598, 89.742100}, {0.599, 89.874100},
           {0.600, 90.006200}, {0.601, 89.965500}, {0.602, 89.924800}, {0.603, 89.884100},
           {0.604, 89.843400}, {0.605, 89.802600}, {0.606, 89.761900}, {0.607, 89.721200},
           {0.608, 89.680500}, {0.609, 89.639800}, {0.610, 89.599100}, {0.611, 89.409100},
           {0.612, 89.219000}, {0.613, 89.029000}, {0.614, 88.838900}, {0.615, 88.648900},
           {0.616, 88.458900}, {0.617, 88.268800}, {0.618, 88.078800}, {0.619, 87.888700},
           {0.620, 87.698700}, {0.621, 87.257700}, {0.622, 86.816700}, {0.623, 86.375700},
           {0.624, 85.934700}, {0.625, 85.493600}, {0.626, 85.052600}, {0.627, 84.611600},
           {0.628, 84.170600}, {0.629, 83.729600}, {0.630, 83.288600}, {0.631, 83.329700},
           {0.632, 83.370700}, {0.633, 83.411800}, {0.634, 83.452800}, {0.635, 83.493900},
           {0.636, 83.535000}, {0.637, 83.576000}, {0.638, 83.617100}, {0.639, 83.658100},
           {0.640, 83.699200}, {0.641, 83.332000}, {0.642, 82.964700}, {0.643, 82.597500},
           {0.644, 82.230200}, {0.645, 81.863000}, {0.646, 81.495800}, {0.647, 81.128500},
           {0.648, 80.761300}, {0.649, 80.394000}, {0.650, 80.026800}, {0.651, 80.045600},
           {0.652, 80.064400}, {0.653, 80.083100}, {0.654, 80.101900}, {0.655, 80.120700},
           {0.656, 80.139500}, {0.657, 80.158300}, {0.658, 80.177000}, {0.659, 80.195800},
           {0.660, 80.214600}, {0.661, 80.420900}, {0.662, 80.627200}, {0.663, 80.833600},
           {0.664, 81.039900}, {0.665, 81.246200}, {0.666, 81.452500}, {0.667, 81.658800},
           {0.668, 81.865200}, {0.669, 82.071500}, {0.670, 82.277800}, {0.671, 81.878400},
           {0.672, 81.479100}, {0.673, 81.079700}, {0.674, 80.680400}, {0.675, 80.281000},
           {0.676, 79.881600}, {0.677, 79.482300}, {0.678, 79.082900}, {0.679, 78.683600},
           {0.680, 78.284200}, {0.681, 77.427900}, {0.682, 76.571600}, {0.683, 75.715300},
           {0.684, 74.859000}, {0.685, 74.002700}, {0.686, 73.146500}, {0.687, 72.290200},
           {0.688, 71.433900}, {0.689, 70.577600}, {0.690, 69.721300}, {0.691, 69.910100},
           {0.692, 70.098900}, {0.693, 70.287600}, {0.694, 70.476400}, {0.695, 70.665200},
           {0.696, 70.854000}, {0.697, 71.042800}, {0.698, 71.231500}, {0.699, 71.420300},
           {0.700, 71.609100}, {0.701, 71.883100}, {0.702, 72.157100}, {0.703, 72.431100},
           {0.704, 72.705100}, {0.705, 72.979000}, {0.706, 73.253000}, {0.707, 73.527000},
           {0.708, 73.801000}, {0.709, 74.075000}, {0.710, 74.349000}, {0.711, 73.074500},
           {0.712, 71.800000}, {0.713, 70.525500}, {0.714, 69.251000}, {0.715, 67.976500},
           {0.716, 66.702000}, {0.717, 65.427500}, {0.718, 64.153000}, {0.719, 62.878500},
           {0.720, 61.604000}, {0.721, 62.432200}, {0.722, 63.260300}, {0.723, 64.088500},
           {0.724, 64.916600}, {0.725, 65.744800}, {0.726, 66.573000}, {0.727, 67.401100},
           {0.728, 68.229300}, {0.729, 69.057400}, {0.730, 69.885600}, {0.731, 70.405700},
           {0.732, 70.925900}, {0.733, 71.446000}, {0.734, 71.966200}, {0.735, 72.486300},
           {0.736, 73.006400}, {0.737, 73.526600}, {0.738, 74.046700}, {0.739, 74.566900},
           {0.740, 75.087000}, {0.741, 73.937600}, {0.742, 72.788100}, {0.743, 71.638700},
           {0.744, 70.489300}, {0.745, 69.339800}, {0.746, 68.190400}, {0.747, 67.041000},
           {0.748, 65.891600}, {0.749, 64.742100}, {0.750, 63.592700}, {0.751, 61.875200},
           {0.752, 60.157800}, {0.753, 58.440300}, {0.754, 56.722900}, {0.755, 55.005400},
           {0.756, 53.288000}, {0.757, 51.570500}, {0.758, 49.853100}, {0.759, 48.135600},
           {0.760, 46.418200}, {0.761, 48.456900}, {0.762, 50.495600}, {0.763, 52.534400},
           {0.764, 54.573100}, {0.765, 56.611800}, {0.766, 58.650500}, {0.767, 60.689200},
           {0.768, 62.728000}, {0.769, 64.766700}, {0.770, 66.805400}, {0.771, 66.463100},
           {0.772, 66.120900}, {0.773, 65.778600}, {0.774, 65.436400}, {0.775, 65.094100},
           {0.776, 64.751800}, {0.777, 64.409600}, {0.778, 64.067300}, {0.779, 63.725100},
           {0.780, 63.382800}, {0.781, 63.474900}, {0.782, 63.567000}, {0.783, 63.659200},
           {0.784, 63.751300}, {0.785, 63.843400}, {0.786, 63.935500}, {0.787, 64.027600},
           {0.788, 64.119800}, {0.789, 64.211900}, {0.790, 64.304000}, {0.791, 63.818800},
           {0.792, 63.333600}, {0.793, 62.848400}, {0.794, 62.363200}, {0.795, 61.877900},
           {0.796, 61.392700}, {0.797, 60.907500}, {0.798, 60.422300}, {0.799, 59.937100},
           {0.800, 59.451900}, {0.801, 58.702600}, {0.802, 57.953300}, {0.803, 57.204000},
           {0.804, 56.454700}, {0.805, 55.705400}, {0.806, 54.956200}, {0.807, 54.206900},
           {0.808, 53.457600}, {0.809, 52.708300}, {0.810, 51.959000}, {0.811, 52.507200},
           {0.812, 53.055300}, {0.813, 53.603500}, {0.814, 54.151600}, {0.815, 54.699800},
           {0.816, 55.248000}, {0.817, 55.796100}, {0.818, 56.344300}, {0.819, 56.892400},
           {0.820, 57.440600}, {0.821, 57.727800}, {0.822, 58.015000}, {0.823, 58.302200},
           {0.824, 58.589400}, {0.825, 58.876500}, {0.826, 59.163700}, {0.827, 59.450900},
           {0.828, 59.738100}, {0.829, 60.025300}, {0.830, 60.312500}});

        return solarRadiation;
    }

    static CSeries getDetectorData()
    {
        CSeries detectorData(std::initializer_list<std::pair<double, double>>({

          {0.380, 0.0000}, {0.385, 0.0001}, {0.390, 0.0001}, {0.395, 0.0002}, {0.400, 0.0004},
          {0.405, 0.0006}, {0.410, 0.0012}, {0.415, 0.0022}, {0.420, 0.0040}, {0.425, 0.0073},
          {0.430, 0.0116}, {0.435, 0.0168}, {0.440, 0.0230}, {0.445, 0.0298}, {0.450, 0.0380},
          {0.455, 0.0480}, {0.460, 0.0600}, {0.465, 0.0739}, {0.470, 0.0910}, {0.475, 0.1126},
          {0.480, 0.1390}, {0.485, 0.1693}, {0.490, 0.2080}, {0.495, 0.2586}, {0.500, 0.3230},
          {0.505, 0.4073}, {0.510, 0.5030}, {0.515, 0.6082}, {0.520, 0.7100}, {0.525, 0.7932},
          {0.530, 0.8620}, {0.535, 0.9149}, {0.540, 0.9540}, {0.545, 0.9803}, {0.550, 0.9950},
          {0.555, 1.0000}, {0.560, 0.9950}, {0.565, 0.9786}, {0.570, 0.9520}, {0.575, 0.9154},
          {0.580, 0.8700}, {0.585, 0.8163}, {0.590, 0.7570}, {0.595, 0.6949}, {0.600, 0.6310},
          {0.605, 0.5668}, {0.610, 0.5030}, {0.615, 0.4412}, {0.620, 0.3810}, {0.625, 0.3210},
          {0.630, 0.2650}, {0.635, 0.2170}, {0.640, 0.1750}, {0.645, 0.1382}, {0.650, 0.1070},
          {0.655, 0.0816}, {0.660, 0.0610}, {0.665, 0.0446}, {0.670, 0.0320}, {0.675, 0.0232},
          {0.680, 0.0170}, {0.685, 0.0119}, {0.690, 0.0082}, {0.695, 0.0057}, {0.700, 0.0041},
          {0.705, 0.0029}, {0.710, 0.0021}, {0.715, 0.0015}, {0.720, 0.0010}, {0.725, 0.0007},
          {0.730, 0.0005}, {0.735, 0.0004}, {0.740, 0.0002}, {0.745, 0.0002}, {0.750, 0.0001},
          {0.755, 0.0001}, {0.760, 0.0001}, {0.765, 0.0000}, {0.770, 0.0000}, {0.775, 0.0000},
          {0.780, 0.0000}}));

        return detectorData;
    }

    [[nodiscard]] static std::vector<double> color5nm()
    {
        return {0.38, 0.385, 0.39, 0.395, 0.4,  0.405, 0.41, 0.415, 0.42, 0.425, 0.43, 0.435,
                0.44, 0.445, 0.45, 0.455, 0.46, 0.465, 0.47, 0.475, 0.48, 0.485, 0.49, 0.495,
                0.5,  0.505, 0.51, 0.515, 0.52, 0.525, 0.53, 0.535, 0.54, 0.545, 0.55, 0.555,
                0.56, 0.565, 0.57, 0.575, 0.58, 0.585, 0.59, 0.595, 0.6,  0.605, 0.61, 0.615,
                0.62, 0.625, 0.63, 0.635, 0.64, 0.645, 0.65, 0.655, 0.66, 0.665, 0.67, 0.675,
                0.68, 0.685, 0.69, 0.695, 0.7,  0.705, 0.71, 0.715, 0.72, 0.725, 0.73, 0.735,
                0.74, 0.745, 0.75, 0.755, 0.76, 0.765, 0.77, 0.775, 0.78};
    };

    static std::shared_ptr<CSpectralSampleData> loadSampleData_NFRC_102()
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

    static std::shared_ptr<CSpectralSampleData> loadSampleData_NFRC_103()
    {
        auto aMeasurements_103 = CSpectralSampleData::create(
          {{0.300, 0.0000, 0.0470, 0.0490}, {0.305, 0.0050, 0.0470, 0.0490},
           {0.310, 0.0000, 0.0470, 0.0480}, {0.315, 0.0030, 0.0460, 0.0480},
           {0.320, 0.0190, 0.0460, 0.0480}, {0.325, 0.0660, 0.0450, 0.0460},
           {0.330, 0.1600, 0.0450, 0.0470}, {0.335, 0.2940, 0.0490, 0.0500},
           {0.340, 0.4370, 0.0550, 0.0560}, {0.345, 0.5660, 0.0620, 0.0620},
           {0.350, 0.6710, 0.0690, 0.0690}, {0.355, 0.7440, 0.0740, 0.0740},
           {0.360, 0.7930, 0.0780, 0.0780}, {0.365, 0.8220, 0.0800, 0.0800},
           {0.370, 0.8320, 0.0810, 0.0810}, {0.375, 0.8190, 0.0800, 0.0800},
           {0.380, 0.8090, 0.0790, 0.0790}, {0.385, 0.8290, 0.0800, 0.0800},
           {0.390, 0.8530, 0.0820, 0.0820}, {0.395, 0.8680, 0.0830, 0.0830},
           {0.400, 0.8750, 0.0830, 0.0830}, {0.410, 0.8750, 0.0830, 0.0830},
           {0.420, 0.8730, 0.0830, 0.0830}, {0.430, 0.8730, 0.0820, 0.0820},
           {0.440, 0.8730, 0.0820, 0.0820}, {0.450, 0.8800, 0.0820, 0.0820},
           {0.460, 0.8870, 0.0820, 0.0820}, {0.470, 0.8900, 0.0820, 0.0820},
           {0.480, 0.8920, 0.0830, 0.0830}, {0.490, 0.8930, 0.0820, 0.0820},
           {0.500, 0.8940, 0.0820, 0.0820}, {0.510, 0.8950, 0.0820, 0.0820},
           {0.520, 0.8950, 0.0820, 0.0820}, {0.530, 0.8940, 0.0820, 0.0820},
           {0.540, 0.8930, 0.0810, 0.0810}, {0.550, 0.8910, 0.0810, 0.0810},
           {0.560, 0.8880, 0.0810, 0.0810}, {0.570, 0.8840, 0.0800, 0.0800},
           {0.580, 0.8810, 0.0800, 0.0800}, {0.590, 0.8760, 0.0790, 0.0790},
           {0.600, 0.8710, 0.0790, 0.0790}, {0.610, 0.8650, 0.0780, 0.0780},
           {0.620, 0.8590, 0.0770, 0.0770}, {0.630, 0.8530, 0.0770, 0.0770},
           {0.640, 0.8470, 0.0760, 0.0760}, {0.650, 0.8400, 0.0750, 0.0750},
           {0.660, 0.8330, 0.0750, 0.0750}, {0.670, 0.8260, 0.0740, 0.0740},
           {0.680, 0.8180, 0.0730, 0.0730}, {0.690, 0.8100, 0.0730, 0.0730},
           {0.700, 0.8020, 0.0720, 0.0720}, {0.710, 0.7940, 0.0710, 0.0720},
           {0.720, 0.7860, 0.0710, 0.0710}, {0.730, 0.7770, 0.0700, 0.0700},
           {0.740, 0.7690, 0.0690, 0.0700}, {0.750, 0.7610, 0.0690, 0.0690},
           {0.760, 0.7520, 0.0680, 0.0680}, {0.770, 0.7440, 0.0670, 0.0680},
           {0.780, 0.7360, 0.0670, 0.0670}, {0.790, 0.7290, 0.0660, 0.0660},
           {0.800, 0.7220, 0.0660, 0.0660}, {0.810, 0.7150, 0.0650, 0.0660},
           {0.820, 0.7100, 0.0650, 0.0650}, {0.830, 0.7020, 0.0640, 0.0650},
           {0.840, 0.6980, 0.0640, 0.0640}, {0.850, 0.6900, 0.0630, 0.0640},
           {0.860, 0.6870, 0.0650, 0.0650}, {0.870, 0.6810, 0.0670, 0.0670},
           {0.880, 0.6770, 0.0650, 0.0660}, {0.890, 0.6730, 0.0660, 0.0660},
           {0.900, 0.6700, 0.0650, 0.0660}, {0.910, 0.6670, 0.0650, 0.0650},
           {0.920, 0.6640, 0.0640, 0.0640}, {0.930, 0.6600, 0.0630, 0.0630},
           {0.940, 0.6580, 0.0640, 0.0640}, {0.950, 0.6560, 0.0630, 0.0630},
           {0.960, 0.6540, 0.0610, 0.0610}, {0.970, 0.6530, 0.0620, 0.0620},
           {0.980, 0.6510, 0.0610, 0.0620}, {0.990, 0.6490, 0.0610, 0.0620},
           {1.000, 0.6480, 0.0590, 0.0600}, {1.050, 0.6450, 0.0590, 0.0600},
           {1.100, 0.6450, 0.0580, 0.0590}, {1.150, 0.6470, 0.0590, 0.0590},
           {1.200, 0.6530, 0.0590, 0.0590}, {1.250, 0.6610, 0.0580, 0.0590},
           {1.300, 0.6730, 0.0600, 0.0600}, {1.350, 0.6870, 0.0600, 0.0600},
           {1.400, 0.7020, 0.0610, 0.0610}, {1.450, 0.7220, 0.0610, 0.0620},
           {1.500, 0.7410, 0.0630, 0.0640}, {1.550, 0.7570, 0.0630, 0.0640},
           {1.600, 0.7690, 0.0650, 0.0650}, {1.650, 0.7750, 0.0650, 0.0640},
           {1.700, 0.7790, 0.0640, 0.0650}, {1.750, 0.7790, 0.0650, 0.0650},
           {1.800, 0.7770, 0.0650, 0.0650}, {1.850, 0.7760, 0.0650, 0.0630},
           {1.900, 0.7730, 0.0620, 0.0620}, {1.950, 0.7730, 0.0650, 0.0650},
           {2.000, 0.7720, 0.0650, 0.0650}, {2.050, 0.7740, 0.0640, 0.0640},
           {2.100, 0.7750, 0.0640, 0.0650}, {2.150, 0.7730, 0.0650, 0.0650},
           {2.200, 0.7580, 0.0640, 0.0650}, {2.250, 0.7590, 0.0640, 0.0640},
           {2.300, 0.7660, 0.0650, 0.0650}, {2.350, 0.7670, 0.0640, 0.0650},
           {2.400, 0.7660, 0.0640, 0.0640}, {2.450, 0.7570, 0.0640, 0.0640},
           {2.500, 0.7500, 0.0630, 0.0630}});

        return aMeasurements_103;
    }

protected:
    virtual void SetUp()
    {
        using FenestrationCommon::MaterialType;
        using FenestrationCommon::WavelengthRange;
        using SingleLayerOptics::BSDFHemisphere;
        using SingleLayerOptics::BSDFBasis;
        using SingleLayerOptics::CBSDFLayerMaker;
        using MultiLayerOptics::CalculationProperties;

        // Create material from samples
        auto thickness = 3.048e-3;   // [m]
        auto aMaterial_102 = SingleLayerOptics::Material::nBandMaterial(
          loadSampleData_NFRC_102(), thickness, MaterialType::Monolithic, WavelengthRange::Solar);
        thickness = 5.715e-3;   // [m]
        auto aMaterial_103 = SingleLayerOptics::Material::nBandMaterial(
          loadSampleData_NFRC_103(), thickness, MaterialType::Monolithic, WavelengthRange::Solar);

        // BSDF definition is needed as well as its material representation
        const auto aBSDF = BSDFHemisphere::create(BSDFBasis::Quarter);
        auto Layer_102 = CBSDFLayerMaker::getSpecularLayer(aMaterial_102, aBSDF);
        auto Layer_103 = CBSDFLayerMaker::getSpecularLayer(aMaterial_103, aBSDF);

        m_Layer = CMultiPaneBSDF::create({Layer_102, Layer_103}, color5nm());

        const CalculationProperties input{
          loadSolarRadiationFile(), loadSolarRadiationFile().getXArray(), getDetectorData()};
        m_Layer->setCalculationProperties(input);
    }

public:
    CMultiPaneBSDF & getLayer()
    {
        return *m_Layer;
    };
};

TEST_F(MultiPaneBSDF_102_103_Visible_QuarterBasis, TestSpecular1)
{
    using FenestrationCommon::Side;
    using FenestrationCommon::PropertySimple;

    const double minLambda = 0.38;
    const double maxLambda = 0.78;

    CMultiPaneBSDF & aLayer = getLayer();

    double tauDiff = aLayer.DiffDiff(minLambda, maxLambda, Side::Front, PropertySimple::T);
    EXPECT_NEAR(0.685078, tauDiff, 1e-6);

    double rhoDiff = aLayer.DiffDiff(minLambda, maxLambda, Side::Front, PropertySimple::R);
    EXPECT_NEAR(0.257564, rhoDiff, 1e-6);

    double theta = 0;
    double phi = 0;

    double tauHem = aLayer.DirHem(minLambda, maxLambda, Side::Front, PropertySimple::T, theta, phi);
    EXPECT_NEAR(0.800011, tauHem, 1e-6);

    double tauDir = aLayer.DirDir(minLambda, maxLambda, Side::Front, PropertySimple::T, theta, phi);
    EXPECT_NEAR(0.800011, tauDir, 1e-6);

    double rhoHem = aLayer.DirHem(minLambda, maxLambda, Side::Front, PropertySimple::R, theta, phi);
    EXPECT_NEAR(0.148032, rhoHem, 1e-6);

    double rhoDir = aLayer.DirDir(minLambda, maxLambda, Side::Front, PropertySimple::R, theta, phi);
    EXPECT_NEAR(0.148032, rhoDir, 1e-6);
}
