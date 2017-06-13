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

              // void align_to(std::string master, size_t min_common);
            static inline bool common(char a, char b) { return a == b && a != 'X' && a != '-'; }
            static std::vector<std::pair<size_t, size_t>> align_to(std::string master, std::string& to_align, size_t min_common);
            void apply_pos_number();

            SeqdbEntrySeq entry_seq;
            std::string amino_acids;
            std::vector<std::pair<size_t, size_t>> pos_number; // to update amino acids in entry_seq.seq

         private:
            static std::pair<size_t, size_t> number_of_common(std::string a, std::string b, size_t start = 0);
            static size_t check_adjust_pos(std::string master, std::string& to_align, size_t adjust_pos);
        };
        using Entries = std::vector<Entry>;

        std::string mVirusType;
        Entries mEntries;
        std::string mMaster;

     private:
        void align_to_master();

        inline void revert() { std::for_each(mEntries.begin(), mEntries.end(), [](auto& entry) { entry.revert(); }); }
        // inline void insert_if(size_t pos, char aa, size_t num_insertions) { std::for_each(mEntries.begin(), mEntries.end(), [&](auto& entry) { entry.insert_if(pos, aa, num_insertions); }); }

    }; // class InsertionsDeletionsDetector

} // namespace seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
