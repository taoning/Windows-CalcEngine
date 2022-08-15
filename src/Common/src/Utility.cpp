#include "Utility.hpp"

namespace FenestrationCommon
{
    IndexRange::IndexRange(size_t startIndex, size_t endIndex) : start(startIndex), end(endIndex)
    {}

    std::vector<IndexRange> chunkIt(size_t start, size_t end, size_t numberOfSplits)
    {
        size_t stepSize{numberOfSplits < (end - start) ? static_cast<size_t>((end - start) / numberOfSplits): 0u};

        std::vector<IndexRange> result;

        size_t currentStart{start};
        size_t currentEnd{currentStart};

        do {
            currentEnd = currentStart + stepSize < end ? currentStart + stepSize + 1u : end + 1u;
            result.emplace_back(currentStart, currentEnd);
            currentStart = currentEnd;
        } while(currentEnd < end + 1u);

        return result;
    }
}   // namespace FenestrationCommon