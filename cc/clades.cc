#include <iostream>
#include <algorithm>

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
            pos162 = static_cast<size_t>(161 - aShift),
            pos163 = static_cast<size_t>(162 - aShift),
            pos165 = static_cast<size_t>(164 - aShift),
            pos166 = static_cast<size_t>(165 - aShift),
            pos167 = static_cast<size_t>(166 - aShift),
            pos172 = static_cast<size_t>(171 - aShift);
    if (aSequence.size() > pos172 && aSequence[pos75] == 'K' && aSequence[pos165] == 'K' && aSequence[pos172] == 'P' && aSequence[pos58] != 'P')
        r.push_back("1A");
    else if (aSequence.size() > pos58 && aSequence[pos58] == 'P')
        r.push_back("1B");
    else
        r.push_back("1");
      // B/Vic always has '-' at pos163
    if (aSequence.size() > pos163 && aSequence[pos162] == '-' && aSequence[pos163] == '-')
        r.push_back("DEL2017"); // B/Vic deletion mutant 2017
    else if (aSequence.size() > pos167 && aSequence[pos163] == '-' && aSequence[pos166] == '-' && aSequence[pos167] == '-')
        r.push_back("TRIPLEDEL2017"); // B/Vic triple deletion mutant 2017
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

// gly: 160S or 160T
// no-gly: not gly
// 3C3: 158N 159F
// 3C3a: 158N 159S
// 3C3b: 62K 83R 158N 261Q
// 3C2a: 158N 159Y
// 3C2a1: 158N 159Y 171K 406V 484E
// defintions below provided by Sarah on 2018-02-17 21:53 based on the CDC TC2 report
// 3C2a1a: 158N 159Y 479E
// 3C2a1b: 92R 158N 159Y 311Q
// 3C2a2: 131K 142K 158N 159Y 261Q
// 3C2a3: 135K 150K 158N 159Y 261Q
// 3C2a4: 31S 53N 142G 144R 158N 159Y 171K 192T 197H

struct PosAA
{
    int pos;
    char aa;
};

struct CladeDesc
{
    const char* clade;
    std::vector<PosAA> pos_aa;
};

#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wexit-time-destructors"
#pragma GCC diagnostic ignored "-Wglobal-constructors"
#endif

static const std::vector<CladeDesc> sClades =
{
    {"3C3", {{158, 'N'}, {159, 'F'}}},
    {"3C3A", {{158, 'N'}, {159, 'S'}}},
    {"3C3B", {/* {62, 'K'}, */ {83, 'R'}, {158, 'N'}, {261, 'Q'}}},
    {"3C2A", {{158, 'N'}, {159, 'Y'}}},
    {"3C2A1", {{158, 'N'}, {159, 'Y'}, {171, 'K'}, {406, 'V'}, {484, 'E'}}},
    {"3C2A1A", {{158, 'N'}, {159, 'Y'}, {479, 'E'}}},
    {"3C2A1B", {{92, 'R'}, {158, 'N'}, {159, 'Y'}, {311, 'Q'}}},
    {"3C2A2", {{131, 'K'}, {142, 'K'}, {158, 'N'}, {159, 'Y'}, {261, 'Q'}}},
    {"3C2A3", {{135, 'K'}, {150, 'K'}, {158, 'N'}, {159, 'Y'}, {261, 'Q'}}},
    {"3C2A4", {{31, 'S'}, {53, 'N'}, {142, 'G'}, {144, 'R'}, {158, 'N'}, {159, 'Y'}, {171, 'K'}, {192, 'T'}, {197, 'H'}}},
// {"GLY: 160S OR 160T
// {"NO-GLY: NOT GLY
};

#pragma GCC diagnostic pop

std::vector<std::string> seqdb::clades_h3n2(std::string aSequence, Shift aShift)
{
    auto has_aa = [aShift,&aSequence](const auto& pos_aa) -> bool { return aa_at(pos_aa.pos, aSequence, aShift) == pos_aa.aa; };

    std::vector<std::string> r;
    for (const auto& clade_desc : sClades) {
        if (std::all_of(clade_desc.pos_aa.begin(), clade_desc.pos_aa.end(), has_aa))
            r.emplace_back(clade_desc.clade);
    }

      // 160S -> gly, 160T -> gly, 160x -> no gly
    if (aa_at(160, aSequence, aShift) == 'S' || aa_at(160, aSequence, aShift) == 'T')
        r.emplace_back("GLY");
    else
        r.emplace_back("NO-GLY");

    return r;

} // clades_h3n2

// std::vector<std::string> seqdb::clades_h3n2(std::string aSequence, Shift aShift)
// {
//     auto r = std::vector<std::string>();
//     if (aa_at(158, aSequence, aShift) == 'N') {
//         switch (aa_at(159, aSequence, aShift)) {
//           case 'F':
//               r.push_back("3C3");
//               break;
//           case 'Y':
//               r.push_back("3C2A");
//               if (aa_at(171, aSequence, aShift) == 'K' && aa_at(406, aSequence, aShift) == 'V' && aa_at(484, aSequence, aShift) == 'E') {
//                     // Derek's message of 2016-12-23 10:32 "clade 3c.2a1"
//                     // 3c.2a1 is a sub-clade of 3c2a
//                   r.push_back("3C2A1");
//               }
//               break;
//           case 'S':
//               r.push_back("3C3A");
//               break;
//           default:
//                 // std::cerr << "@159: " << aSequence[pos159] << std::endl;
//               break;
//         }
//     }

//     if (aa_at(159, aSequence, aShift) == 'F' && aa_at(83, aSequence, aShift) == 'R' && aa_at(261, aSequence, aShift) == 'Q') // && aa_at(62, aSequence, aShift) == 'K')
//         r.push_back("3C3B");
//       // if (aSequence.size()  > pos261 && aSequence[pos62] == 'K' && aSequence[pos83] == 'R' && aSequence[pos261] == 'Q')
//       //     r.push_back("3C3b?");

//       // 160S -> gly, 160T -> gly, 160x -> no gly
//     if (aa_at(160, aSequence, aShift) == 'S' || aa_at(160, aSequence, aShift) == 'T')
//         r.push_back("GLY");
//     else
//         r.push_back("NO-GLY");

//     return r;

// } // clades_h3n2

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
