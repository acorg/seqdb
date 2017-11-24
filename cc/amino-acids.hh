#pragma once

#include <string>
#include <set>

#include "seqdb/messages.hh"
#include "seqdb/sequence-shift.hh"

// ----------------------------------------------------------------------

namespace seqdb
{

    constexpr size_t MINIMUM_SEQUENCE_AA_LENGTH = 200; // actually H3 3C3b clade requires 261Q

// ----------------------------------------------------------------------

    struct AlignData
    {
        inline AlignData() : shift(Shift::AlignmentFailed) {}
          // inline AlignData(const AlignData&) = default;
        inline AlignData(std::string aSubtype, std::string aLineage, std::string aGene, Shift aShift)
            : subtype(aSubtype), lineage(aLineage), gene(aGene), shift(aShift) {}

        std::string subtype;
        std::string lineage;
        std::string gene;
        Shift shift;
    };

    struct AlignAminoAcidsData : public AlignData
    {
        inline AlignAminoAcidsData() = default;
        inline AlignAminoAcidsData(const AlignData& aAlignData, std::string aAminoAcids, int aOffset)
            : AlignData(aAlignData), amino_acids(aAminoAcids), offset(aOffset) {}
        inline AlignAminoAcidsData(AlignData&& aAlignData)
            : AlignData(aAlignData), offset(0) {}
          // inline AlignAminoAcidsData(const AlignAminoAcidsData&) = default;
          // inline AlignAminoAcidsData(AlignAminoAcidsData&&) = default;
          // inline AlignAminoAcidsData& operator=(const AlignAminoAcidsData&) = default;

        static inline AlignAminoAcidsData alignment_failed()
            {
                return AlignAminoAcidsData(AlignData(), std::string(), 0);
            }

        std::string amino_acids;
        int offset;
    };

// ----------------------------------------------------------------------

    AlignAminoAcidsData translate_and_align(std::string aNucleotides, Messages& aMessages);

    std::string translate_nucleotides_to_amino_acids(std::string aNucleotides, size_t aOffset, Messages& aMessages);
    AlignData align(std::string_view aAminoAcids, Messages& aMessages);

// ----------------------------------------------------------------------

    inline AlignAminoAcidsData align_amino_acids(std::string aAminoAcids, Messages& aMessages)
    {
        return AlignAminoAcidsData(align(aAminoAcids, aMessages));
    }

// ----------------------------------------------------------------------

    inline bool is_nucleotides(std::string aSequence)
    {
        constexpr char sNucleotideElements[] = "-ABCDGHKMNRSTUVWY"; // https://en.wikipedia.org/wiki/Nucleic_acid_notation
        const std::set<char> elements(aSequence.begin(), aSequence.end());
        return std::includes(std::begin(sNucleotideElements), std::end(sNucleotideElements), elements.begin(), elements.end());
    }

// ----------------------------------------------------------------------

} // namespace seqdb

inline std::ostream& operator << (std::ostream& out, const seqdb::AlignData& a)
{
    return out << "[AlignData " << a.subtype << " " << a.lineage << " " << a.gene << " shift:" << a.shift << "]";
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
