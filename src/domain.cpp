#include "mdsim/domain.hpp"

namespace mdsim {

std::vector<BlockRange> partitionRange(int total, int nblocks) {
    std::vector<BlockRange> out;
    if (nblocks < 1) nblocks = 1;
    if (nblocks > total) nblocks = total > 0 ? total : 1;
    out.resize(nblocks);
    for (int i = 0; i < nblocks; ++i) {
        out[i].begin = static_cast<int>((static_cast<long long>(i) * total) / nblocks);
        out[i].end   = static_cast<int>((static_cast<long long>(i + 1) * total) / nblocks);
    }
    return out;
}

} // namespace mdsim
