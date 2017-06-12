#pragma once

#include "seqdb/seqdb.hh"

namespace seqdb
{
    class InsertionsDeletionsDetector
    {
     public:
        InsertionsDeletionsDetector(Seqdb& aSeqdb, std::string aVirusType);

        void detect();

        class Entry
        {
         public:
            inline Entry(SeqdbEntrySeq&& aEntrySeq) : entry_seq(aEntrySeq), amino_acids(entry_seq.seq().amino_acids(true)) {}
            inline void revert() { amino_acids = entry_seq.seq().amino_acids(true); }
            inline void insert_if(size_t pos, char aa, size_t num_insertions) { if (amino_acids[pos] == aa) amino_acids.insert(pos, num_insertions, '-'); }

            void align_to(std::string master, size_t min_common);
            size_t find_adjust_pos(std::string master) const;
            size_t number_of_common(std::string master) const;
            void apply_pos_number();
            static inline bool common(char a, char b) { return a == b && a != 'X' && a != '-'; }

            SeqdbEntrySeq entry_seq;
            std::string amino_acids;
            std::vector<std::pair<size_t, size_t>> pos_number; // to update amino acids in entry_seq.seq
        };
        using Entries = std::vector<Entry>;

        std::string mVirusType;
        Entries mEntries;
        std::string mLongest;

     private:
        void align_to_longest();

        inline void revert() { std::for_each(mEntries.begin(), mEntries.end(), [](auto& entry) { entry.revert(); }); }
        inline void insert_if(size_t pos, char aa, size_t num_insertions) { std::for_each(mEntries.begin(), mEntries.end(), [&](auto& entry) { entry.insert_if(pos, aa, num_insertions); }); }

    }; // class InsertionsDeletionsDetector

} // namespace seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
