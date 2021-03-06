#pragma once

#include <string>
#include <vector>

#include "seqdb/sequence-shift.hh"

// ----------------------------------------------------------------------

namespace seqdb
{
    std::vector<std::string> clades_b_yamagata(std::string_view aSequence, Shift aShift, std::string_view aName);
    std::vector<std::string> clades_b_victoria(std::string_view aSequence, Shift aShift, std::string_view aName);
    std::vector<std::string> clades_h1pdm(std::string_view aSequence, Shift aShift, std::string_view aName);
    std::vector<std::string> clades_h3n2(std::string_view aSequence, Shift aShift, std::string_view aName);

      // Note aPos in a human notation, i.e. starts with 1
    inline char aa_at(int aPos, std::string_view aSequence, Shift aShift)
    {
        const size_t offset = static_cast<size_t>(aPos - 1 - aShift);
        return aSequence.size() > offset ? aSequence[offset] : ' ';
    }
}

// ----------------------------------------------------------------------
