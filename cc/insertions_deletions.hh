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

            static inline bool common(char a, char b) { return a == b && a != 'X' && a != '-'; }
            static std::vector<std::pair<size_t, size_t>> align_to(std::string master, std::string& to_align, size_t min_common);
            void apply_pos_number();

            SeqdbEntrySeq entry_seq;
            std::string amino_acids;
            std::vector<std::pair<size_t, size_t>> pos_number; // to update amino acids in entry_seq.seq
        };
        using Entries = std::vector<Entry>;

        std::string mVirusType;
        Entries mEntries;
        std::string mMaster;

     private:
        void align_to_master();

        inline void revert() { std::for_each(mEntries.begin(), mEntries.end(), [](auto& entry) { entry.revert(); }); }

    }; // class InsertionsDeletionsDetector

      // ----------------------------------------------------------------------

    class BLineageDetector
    {
     public:
        inline BLineageDetector(Seqdb& aSeqdb) : mSeqdb(aSeqdb) {}
        void detect();

     private:
        Seqdb& mSeqdb;

    }; // class BLineageDetector

} // namespace seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
