#pragma once

#include <iostream>
#include <string>

#include "seqdb/messages.hh"
#include "seqdb/sequence-shift.hh"

// ----------------------------------------------------------------------

namespace seqdb
{

    constexpr size_t NORMAL_SEQUENCE_AA_LENGTH_H1 = 549;
    constexpr size_t NORMAL_SEQUENCE_AA_LENGTH_H3 = 550;
    constexpr size_t NORMAL_SEQUENCE_AA_LENGTH_B = 570;

    constexpr size_t MINIMUM_SEQUENCE_AA_LENGTH = 400; // throw away everything shorter
    constexpr size_t MINIMUM_SEQUENCE_NUC_LENGTH = MINIMUM_SEQUENCE_AA_LENGTH * 3; // throw away everything shorter

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

    AlignAminoAcidsData translate_and_align(std::string aNucleotides, Messages& aMessages, std::string name);

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
        constexpr char sNucleotideElements[] = "-ABCDGHKMNRSTUVWXY"; // https://en.wikipedia.org/wiki/Nucleic_acid_notation + X (gisaid nuc seqs contain X)
        std::string sorted_seq = aSequence;
        std::sort(std::begin(sorted_seq), std::end(sorted_seq));
        const auto last = std::unique(std::begin(sorted_seq), std::end(sorted_seq));
        return std::includes(std::begin(sNucleotideElements), std::end(sNucleotideElements), begin(sorted_seq), last);
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
