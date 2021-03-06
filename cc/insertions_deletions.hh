#pragma once

#include "seqdb/seqdb.hh"

namespace seqdb
{
    class InsertionsDeletionsDetector
    {
     public:
        InsertionsDeletionsDetector(Seqdb& aSeqdb, std::string_view aVirusType);

        void detect();

        class Entry
        {
         public:
            Entry(SeqdbEntrySeq&& aEntrySeq) : entry_seq(aEntrySeq), amino_acids(entry_seq.seq().amino_acids(true)) {}
            void revert() { amino_acids = entry_seq.seq().amino_acids(true); }
              // void insert_if(size_t pos, char aa, size_t num_insertions) { if (amino_acids[pos] == aa) amino_acids.insert(pos, num_insertions, '-'); }

            static inline bool common(char a, char b) { return a == b && a != 'X' && a != '-'; }
            static std::vector<std::pair<size_t, size_t>> align_to(const std::string_view master, std::string& to_align, const SeqdbEntrySeq& entry_seq);
            void apply_pos_number();

            SeqdbEntrySeq entry_seq;
            std::string amino_acids;
            std::vector<std::pair<size_t, size_t>> pos_number; // to update amino acids in entry_seq.seq
        };
        using Entries = std::vector<Entry>;

        std::string mVirusType;
        Entries mEntries;
        std::string mMaster;
        bool master_switching_allowed_ = true;

     private:
        void align_to_master();

        void revert() { std::for_each(mEntries.begin(), mEntries.end(), [](auto& entry) { entry.revert(); }); }
        void choose_master();

    }; // class InsertionsDeletionsDetector

      // ----------------------------------------------------------------------

    class BLineageDetector
    {
     public:
        BLineageDetector(Seqdb& aSeqdb) : mSeqdb(aSeqdb) {}
        void detect();

     private:
        Seqdb& mSeqdb;

    }; // class BLineageDetector

} // namespace seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
