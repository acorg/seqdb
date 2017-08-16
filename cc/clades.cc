#include <iostream>

#include "clades.hh"

using namespace seqdb;

// ----------------------------------------------------------------------

std::vector<std::string> seqdb::clades_b_yamagata(std::string aSequence, Shift aShift)
{
      // 165N -> Y2, 165Y -> Y3 (yamagata numeration, 163 is not -)
      // 166N -> Y2, 166Y -> Y3 (victoria numeration, 163 is -)
    auto r = std::vector<std::string>();
    auto const pos = 165 - aShift; // use victoria numeration
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

std::vector<std::string> seqdb::clades_b_victoria(std::string aSequence, Shift aShift)
{
      // mark clade 1A (all strains if empty clade removed, otherwise clade 1B defined by L58P, clade 1 defined by N75K, N165K, S172P)

      // 75K, 165K, 172P -> 1, 58P -> 1B, ? -> 1A
    auto r = std::vector<std::string>();
    const size_t pos58 = static_cast<size_t>(57 - aShift),
            pos75 = static_cast<size_t>(74 - aShift),
            pos165 = static_cast<size_t>(164 - aShift),
            pos172 = static_cast<size_t>(171 - aShift);
    if (aSequence.size() > pos172 && aSequence[pos75] == 'K' && aSequence[pos165] == 'K' && aSequence[pos172] == 'P' && aSequence[pos58] != 'P')
        r.push_back("1A");
    else if (aSequence.size() > pos58 && aSequence[pos58] == 'P')
        r.push_back("1B");
    else
        r.push_back("1");
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
    if (aa_at(158, aSequence, aShift) == 'N') {
        switch (aa_at(159, aSequence, aShift)) {
          case 'F':
              r.push_back("3C3");
              break;
          case 'Y':
              r.push_back("3C2a");
              if (aa_at(171, aSequence, aShift) == 'K' && aa_at(406, aSequence, aShift) == 'V' && aa_at(484, aSequence, aShift) == 'E') {
                    // Derek's message of 2016-12-23 10:32 "clade 3c.2a1"
                  r.push_back("3C2a1");
              }
              break;
          case 'S':
              r.push_back("3C3a");
              break;
          default:
                // std::cerr << "@159: " << aSequence[pos159] << std::endl;
              break;
        }
    }

    if (aa_at(159, aSequence, aShift) == 'F' && aa_at(83, aSequence, aShift) == 'R' && aa_at(261, aSequence, aShift) == 'Q') // && aa_at(62, aSequence, aShift) == 'K')
        r.push_back("3C3b");
      // if (aSequence.size()  > pos261 && aSequence[pos62] == 'K' && aSequence[pos83] == 'R' && aSequence[pos261] == 'Q')
      //     r.push_back("3C3b?");

      // 160S -> gly, 160T -> gly, 160x -> no gly
    if (aa_at(160, aSequence, aShift) == 'S' || aa_at(160, aSequence, aShift) == 'T')
        r.push_back("gly");
    else
        r.push_back("no-gly");

    return r;

} // clades_h3n2

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
