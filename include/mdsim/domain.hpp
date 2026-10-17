// Static domain decomposition.
//
// Space (and the cell list that indexes it) is split into contiguous blocks
// assigned to workers with static scheduling. Because cells are never split
// and every particle's force is summed by a single owning block in a fixed
// order, results are independent of the number of threads.
#pragma once

#include <vector>

namespace mdsim {

struct BlockRange {
    int begin;
    int end;
};

// Split [0, total) into nblocks contiguous ranges.
std::vector<BlockRange> partitionRange(int total, int nblocks);

} // namespace mdsim
