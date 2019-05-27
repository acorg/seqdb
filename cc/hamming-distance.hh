#pragma once

namespace acmacs::seqdb
{
    template <typename S1, typename S2> inline size_t hamming_distance(S1&& first, S2&& second)
    {
        const auto size = std::min(first.size(), second.size());
        size_t dist = 0;
        for (size_t index = 0; index < size; ++index) {
            if (first[index] != second[index])
                ++dist;
        }
        return dist;
    }

} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
