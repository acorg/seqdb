#include <iostream>
#include <algorithm>

#include "clades.hh"

using namespace seqdb;

// ----------------------------------------------------------------------

std::vector<std::string> seqdb::clades_b_yamagata(std::string aSequence, Shift aShift, std::string /*aName*/)
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

std::vector<std::string> seqdb::clades_b_victoria(std::string aSequence, Shift aShift, std::string aName)
{
    auto r = std::vector<std::string>();
    const size_t
            pos58 = static_cast<size_t>(57 - aShift),
            pos75 = static_cast<size_t>(74 - aShift),
            pos162 = static_cast<size_t>(161 - aShift),
            pos163 = static_cast<size_t>(162 - aShift),
            pos164 = static_cast<size_t>(163 - aShift),
              // pos165 = static_cast<size_t>(164 - aShift),
              // pos166 = static_cast<size_t>(165 - aShift),
              // pos167 = static_cast<size_t>(166 - aShift),
            pos172 = static_cast<size_t>(171 - aShift);

      // before 2018-09-03
      // if (aSequence.size() > pos172 && aSequence[pos75] == 'K' && aSequence[pos165] == 'K' && aSequence[pos172] == 'P' && aSequence[pos58] != 'P')
      //     r.push_back("1A");

      // 2018-09-03, Sarah: clades should (technically) be defined by a phylogenetic tree rather than a set of amino acids
    if (aSequence.size() > pos172 && aSequence[pos75] == 'K' && aSequence[pos172] == 'P' && aSequence[pos58] != 'P')
        r.push_back("1A");
    else if (aSequence.size() > pos58 && aSequence[pos58] == 'P')
        r.push_back("1B");
    else
        r.push_back("1");

    if (aSequence.size() > pos164 && aSequence[pos162] == '-' && aSequence[pos163] == '-' && aSequence[pos164] == '-')
        r.push_back("TRIPLEDEL2017"); // B/Vic triple deletion mutant 2017 162,163,164 by convention
    else if (aSequence.size() > pos163 && aSequence[pos162] == '-' && aSequence[pos163] == '-')
        r.push_back("DEL2017"); // B/Vic deletion mutant 2017
    else if (aSequence.size() > pos164 && (aSequence[pos162] == '-' || aSequence[pos163] == '-' || aSequence[pos164] == '-'))
        std::cerr << "WARNING: [" << aName << "]: strange B/Vic deletion mutant: " << aSequence << '\n';
    return r;

      // earlier definitions:
      // mark clade 1A (all strains if empty clade removed, otherwise clade 1B defined by L58P, clade 1 defined by N75K, N165K, S172P)

} // clades_b_yamagata

// ----------------------------------------------------------------------

std::vector<std::string> seqdb::clades_h1pdm(std::string aSequence, Shift aShift, std::string /*aName*/)
{
    // ----------------------------------------------------------------------
    // 2018-09-19 clade definitions changed by Sarah before SSM
    // ----------------------------------------------------------------------
    // 6B: 163Q
    // 6B1: 162N, 163Q
    // 6B2: 152T, 163Q
    auto r = std::vector<std::string>();
    const auto pos152 = static_cast<size_t>(151 - aShift),
            pos162 = static_cast<size_t>(161 - aShift),
            pos163 = static_cast<size_t>(162 - aShift);
    if (pos163 > 0 && aSequence.size() > pos163) {
        if (aSequence[pos163] == 'Q') {
            r.push_back("6B");
            if (aSequence[pos162] == 'N')
                r.push_back("6B1");
            if (aSequence[pos152] == 'T')
                r.push_back("6B2");
        }
    }
    return r;

    // ----------------------------------------------------------------------
    // Before 2018-09-19
    // ----------------------------------------------------------------------
    //   // 84N+162N+216T - 6B.1, 152T+173I+501E - 6B.2
    //   // ? 156 (see A/PUERTO RICO/15/2018 of CDC:20180511)
    // auto r = std::vector<std::string>();
    // auto const pos84i = 83 - aShift;
    // if (pos84i > 0) {
    //     const size_t pos84 = static_cast<size_t>(pos84i);
    //     if (aSequence.size() > pos84 && aSequence[pos84] == 'N') {
    //         const size_t pos162 = static_cast<size_t>(161 - aShift);
    //         const size_t pos216 = static_cast<size_t>(215 - aShift);
    //         if (aSequence.size() > pos216 && aSequence[pos162] == 'N' && aSequence[pos216] == 'T')
    //             r.push_back("6B1");
    //     }

    //     const size_t pos152 = static_cast<size_t>(151 - aShift);
    //     const size_t pos173 = static_cast<size_t>(172 - aShift);
    //     const size_t pos501 = static_cast<size_t>(500 - aShift);
    //     if (aSequence.size() > pos501 && aSequence[pos152] == 'T' && aSequence[pos173] == 'I' && aSequence[pos501] == 'E')
    //         r.push_back("6B2");
    // }
    // return r;

} // clades_h1pdm

// ----------------------------------------------------------------------

// HK68 TY155
// EN72 Transition QK189
// VI75 Transition GE158+DN193
// TX77 Transition KE156
// BK79 Transition SY159 (within cluster)+YH155+KR189
// SI87 Transition NK145 to BE89
// SI87 Transition EK156 + SD133(for direction) to BE92
// BE92 Transition NK145 to WU95
// WU95 Transition KQ156+EK158 to SY97
// SY97 Transition HT155+QH156 to FU02

// Fu02-> Cal04 most likely KN145
// CA04 - WI05: K145N
// WI05 - PE09: S193F + K158N
// PE09 - SW13: F159S
// PE09 - HK14: F159Y

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

// clades per nextflu (2018-08-07, by Sarah):
// "3b":    [('HA2',158,'N'), ('HA1',198,'S'), ('HA1',312,'S'), ('HA1',223,'I'),     ('HA1',145,'S')],
// "3c":    [('HA1',45,'N'), ('HA1',48,'I'), ('nuc',456,'T'), ('HA1',198,'S'),     ('HA1',312,'S'), ('HA1',223,'I')],
// "3c2":   [('HA2',160,'N'), ('HA1',145,'S'), ('nuc',693,'A'), ('nuc',1518,'G')],
// "3c3":   [('nuc',285,'T'), ('nuc',430,'G'), ('nuc',472,'G'), ('nuc',1296,'A')],
// "3c3.A": 128A 142G 159S
// "3c2.A": [('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',225,'D'),     ('HA1',311,'H'), ('nuc',1491,'A'), ('nuc', 234, 'A')],
// "A1": 144S 159Y 171K 225D 311H   ('nuc',1491,'A'), ('HA2',77,'V'), ('HA2',155,'E'), ('HA2',160,'N')], # formerly clade 3c2.a1
// "3c2":   [('HA1',144,'N'), ('HA1',159,'F'), ('HA1',225,'N'), ('HA1',160,'T'), ('HA1',142,'R')],
// "3c3.B": [('HA1',83,'R'), ('HA1',261,'Q'), ('HA1',62,'K'),    ('HA1',122,'D')],
// "A2": 131K 142K 144S 159Y 225D 261Q 311H  ('nuc',1491,'A'), ('nuc', 234, 'A'),     ], # formerly clade 3
// "A3": 121K 144K 159Y 225D 311H     ('nuc',1491,'A'), ('nuc',234,'A'), ], # formerly clade 2
// "A4": 53N 144R 159Y 171K 192T 197H 225D 311H    ('nuc',1491,'A'), ('nuc',234,'A'), ] # formerly clade 1

// Clades defined by Sarah on 2018-08-09 for TC1 based on the tree made for TC1 and based on the older definitions above
// 3C.3    158N 159F
// 3C.3a   158N 159S
// 3C.3b   62K 83R 158N 261Q
// 3C.2a   158N 159Y
// 2a1     3I 158N 159Y 171K 406V 484E   (3C2a1)
// 2a1a    92K 158N 159Y 479E (3C2a1a)  (92K added by Sarah later: to eliminate the extra 2a1a label in the tree for TC1 2018-08-09)
// 2a1b    92R 158N 159Y 311Q (3C2a1b)
// 2a2     131K 142K 158N 159Y 261Q (3C2a2)
// 2a3     3I 121K 144K (former 3C2a3)
// 2a4     31S 53N 142G 144R 158N 159Y 171K 192T 197H

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
    {"3C.3", {{158, 'N'}, {159, 'F'}}},
    {"3C.3A", {{158, 'N'}, {159, 'S'}}},
    {"3C.3B", {{62, 'K'}, {83, 'R'}, {158, 'N'}, {261, 'Q'}}},
    {"3C.2A", {{158, 'N'}, {159, 'Y'}}},
    {"2A1", {{3, 'I'}, {158, 'N'}, {159, 'Y'}, {171, 'K'}, {406, 'V'}, {484, 'E'}}},
    {"2A1A", {{92, 'K'}, {158, 'N'}, {159, 'Y'}, {479, 'E'}}},
    {"2A1B", {{92, 'R'}, {158, 'N'}, {159, 'Y'}, {311, 'Q'}}},
    {"2A2", {{131, 'K'}, {142, 'K'}, {158, 'N'}, {159, 'Y'}, {261, 'Q'}}},
    {"2A3", {{3, 'I'}, {121, 'K'}, {144, 'K'}}},
    {"2A4", {{31, 'S'}, {53, 'N'}, {142, 'G'}, {144, 'R'}, {158, 'N'}, {159, 'Y'}, {171, 'K'}, {192, 'T'}, {197, 'H'}}},
};

#pragma GCC diagnostic pop

std::vector<std::string> seqdb::clades_h3n2(std::string aSequence, Shift aShift, std::string /*aName*/)
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
