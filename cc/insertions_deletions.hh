#pragma once

#include "seqdb/seqdb.hh"

namespace seqdb
{
    class InsertionsDeletionsDetector
    {
     public:
        InsertionsDeletionsDetector(Seqdb& aSeqdb, std::string aVirusType);

        void detect();

     private:
        class Entry
        {
         public:
            inline Entry(SeqdbEntrySeq&& aEntrySeq) : entry_seq(aEntrySeq), amino_acids(entry_seq.seq().amino_acids(true)) {}
            SeqdbEntrySeq entry_seq;
            std::string amino_acids;
        };

        std::vector<Entry> mEntries;

    }; // class InsertionsDeletionsDetector

} // namespace seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
