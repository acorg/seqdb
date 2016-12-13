#include <iostream>

#include "clades.hh"

using namespace seqdb;

// ----------------------------------------------------------------------

std::vector<std::string> seqdb::clades_b_yamagata(std::string aSequence, Shift aShift)
{
      // 165N -> Y2, 165Y -> Y3
    auto r = std::vector<std::string>();
    auto const pos = 164 - aShift;
    if (aSequence.size() > static_cast<size_t>(pos) && pos > 0) {
        switch (aSequence[static_cast<size_t>(pos)]) {
          case 'N':
              r.push_back("Y2");
              break;
          case 'Y':
              r.push_back("Y3");
              break;
          default:
              break;
        }
    }
    return r;

} // clades_b_yamagata

// ----------------------------------------------------------------------

std::vector<std::string> seqdb::clades_h1pdm(std::string aSequence, Shift aShift)
{
      // 84N+162N+216T - 6B.1, 152T+173I+501E - 6B.2
    auto r = std::vector<std::string>();
    auto const pos84i = 83 - aShift;
    if (pos84i > 0) {
        const size_t pos84 = static_cast<size_t>(pos84i);
        if (aSequence.size() > pos84 && aSequence[pos84] == 'N') {
            const size_t pos162 = static_cast<size_t>(161 - aShift);
            const size_t pos216 = static_cast<size_t>(215 - aShift);
            if (aSequence.size() > pos216 && aSequence[pos162] == 'N' && aSequence[pos216] == 'T')
                r.push_back("6B1");
        }

        const size_t pos152 = static_cast<size_t>(151 - aShift);
        const size_t pos173 = static_cast<size_t>(172 - aShift);
        const size_t pos501 = static_cast<size_t>(500 - aShift);
        if (aSequence.size() > pos501 && aSequence[pos152] == 'T' && aSequence[pos173] == 'I' && aSequence[pos501] == 'E')
            r.push_back("6B2");
    }
    return r;

} // clades_h1pdm

// ----------------------------------------------------------------------

std::vector<std::string> seqdb::clades_h3n2(std::string aSequence, Shift aShift)
{
      // 158N, 159F -> 3C3, 159Y -> 3c2a, 159S -> 3c3a, 62K+83R+261Q -> 3C3b.
    auto r = std::vector<std::string>();
    auto const pos158i = 157 - aShift;
    if (pos158i > 0) {
        const size_t pos158 = static_cast<size_t>(pos158i);
        const size_t pos159 = static_cast<size_t>(158 - aShift);
        if (aSequence.size() > pos159 && aSequence[pos158] == 'N') {
            switch (aSequence[pos159]) {
              case 'F':
                  r.push_back("3C3");
                  break;
              case 'Y':
                  r.push_back("3C2a");
                  break;
              case 'S':
                  r.push_back("3C3a");
                  break;
              default:
                    // std::cerr << "@159: " << aSequence[pos159] << std::endl;
                  break;
            }
        }

          // auto const pos62 = 61 - aShift;
        const size_t pos83 =  static_cast<size_t>(82 - aShift);
        const size_t pos261 = static_cast<size_t>(260 - aShift);
        if (aSequence.size() > pos261 && aSequence[pos159] == 'F' && aSequence[pos83] == 'R' && aSequence[pos261] == 'Q') // && aSequence[pos62] == 'K')
            r.push_back("3C3b");
          // if (aSequence.size()  > pos261 && aSequence[pos62] == 'K' && aSequence[pos83] == 'R' && aSequence[pos261] == 'Q')
          //     r.push_back("3C3b?");

          // 160S -> gly, 160T -> gly, 160x -> no gly
        const size_t pos160 = static_cast<size_t>(159 - aShift);
        if (aSequence.size() > pos160 && pos160 > 0 && (aSequence[pos160] == 'S' || aSequence[pos160] == 'T'))
            r.push_back("gly");
        else
            r.push_back("no-gly");
    }
    return r;

} // clades_h3n2

// ----------------------------------------------------------------------
