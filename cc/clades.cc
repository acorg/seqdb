#include <iostream>
#include <algorithm>
#include <array>
#include <vector>

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

// https://notebooks.antigenic-cartography.org/eu/results/eu/2019-0118-clades/clades.org
// /scp:albertine:/syn/eu/ac/results/eu/2019-0118-clades/clades.org

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

static const std::array sClades
{
    CladeDesc{"3C.3",  {{158, 'N'}, {159, 'F'}}},
    CladeDesc{"3A",    {{138, 'S'}, {159, 'S'}, {225, 'D'}, {326, 'R'}}},
    CladeDesc{"3B",    {{ 62, 'K'}, { 83, 'R'}, {261, 'Q'}}},
    CladeDesc{"2A",    {{158, 'N'}, {159, 'Y'}}},
    CladeDesc{"2A1",   {{158, 'N'}, {159, 'Y'}, {171, 'K'}, {406, 'V'}, {484, 'E'}}},
    CladeDesc{"2A1A",  {{121, 'K'}, {135, 'K'}, {158, 'N'}, {159, 'Y'}, {171, 'K'}, {406, 'V'}, {479, 'E'}, {484, 'E'}}},
    CladeDesc{"2A1B",  {{ 92, 'R'}, {121, 'K'}, {158, 'N'}, {159, 'Y'}, {171, 'K'}, {311, 'Q'}, {406, 'V'}, {484, 'E'}}},
    CladeDesc{"2A2",   {{131, 'K'}, {142, 'K'}, {158, 'N'}, {159, 'Y'}, {261, 'Q'}}},
    CladeDesc{"2A3",   {{121, 'K'}, {135, 'K'}, {144, 'K'}, {150, 'K'}, {158, 'N'}, {159, 'Y'}, {261, 'Q'}}},
    CladeDesc{"2A4",   {{ 31, 'S'}, { 53, 'N'}, {142, 'G'}, {144, 'R'}, {158, 'N'}, {159, 'Y'}, {171, 'K'}, {192, 'T'}, {197, 'H'}}},
    CladeDesc{"GLY",   {{160, 'S'}}},
    CladeDesc{"GLY",   {{160, 'T'}}},
    CladeDesc{"159S",  {{159, 'S'}}}, // explicit Derek's request on 2019-04-18
    CladeDesc{"159F",  {{159, 'F'}}}, // explicit Derek's request on 2019-04-18
    CladeDesc{"159Y",  {{159, 'Y'}}}, // explicit Derek's request on 2019-04-18
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

    return r;

} // clades_h3n2

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
