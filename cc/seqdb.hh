#pragma once

#include <iostream>
#include <string>
#include <functional>
#include <algorithm>
#include <regex>
#include <iterator>
#include <deque>
#include <map>
#include <vector>
#include <numeric>

#include "acmacs-base/stream.hh"
#include "acmacs-base/name-encode.hh"
#include "hidb/hidb.hh"
#include "seqdb/sequence-shift.hh"
#include "seqdb/amino-acids.hh"
#include "seqdb/messages.hh"

// ----------------------------------------------------------------------

class Antigens;                 // acmacs-chart

namespace seqdb
{
    class Seqdb;
    class SeqdbIterator;

// ----------------------------------------------------------------------

    class SequenceNotAligned : public std::runtime_error
    {
     public:
        inline SequenceNotAligned(std::string prefix) : std::runtime_error(prefix + ": sequence not aligned") {}
    };

// ----------------------------------------------------------------------

    class SeqdbSeq
    {
     public:
        using LabIds = std::map<std::string, std::vector<std::string>>;

        inline SeqdbSeq() : mGene("HA") {}

        inline SeqdbSeq(std::string aSequence, std::string aGene)
            : SeqdbSeq()
            {
                if (is_nucleotides(aSequence))
                    mNucleotides = aSequence;
                else
                    mAminoAcids = aSequence;
                if (!aGene.empty())
                    mGene = aGene;
            }

        // inline SeqdbSeq(bool aNucs, std::string aSequence, std::string aGene)
        //     : SeqdbSeq()
        //     {
        //         if (aNucs)
        //             mNucleotides = aSequence;
        //         else
        //             mAminoAcids = aSequence;
        //         if (!aGene.empty())
        //             mGene = aGene;
        //     }

        // inline SeqdbSeq(std::string aNucleotides, std::string aAminoAcids, std::string aGene)
        //     : SeqdbSeq()
        //     {
        //         mNucleotides = aNucleotides;
        //         mAminoAcids = aAminoAcids;
        //         if (!aGene.empty())
        //             mGene = aGene;
        //     }

        // inline SeqdbSeq(std::string aNucleotides, std::string aGene)
        //     : SeqdbSeq(aNucleotides, std::string(), aGene)
        //     {
        //     }

        AlignAminoAcidsData align(bool aForce, Messages& aMessages);

          // returns if aNucleotides matches mNucleotides or aAminoAcids matches mAminoAcids
        bool match_update(const SeqdbSeq& aNewSeq);
        bool match_update_nucleotides(const SeqdbSeq& aNewSeq);
        bool match_update_amino_acids(const SeqdbSeq& aNewSeq);
        void add_passage(std::string aPassage);
        void update_gene(std::string aGene, Messages& aMessages, bool replace_ha = false);
        void add_reassortant(std::string aReassortant);
        void add_lab_id(std::string aLab, std::string aLabId);
        void update_clades(std::string aVirusType, std::string aLineage);
        inline const std::vector<std::string>& clades() const { return mClades; }
        inline std::vector<std::string>& clades() { return mClades; }
        inline bool has_clade(std::string aClade) const { return std::find(std::begin(mClades), std::end(mClades), aClade) != std::end(mClades); }

        inline bool is_short() const { return mAminoAcids.empty() ? mNucleotides.size() < (MINIMUM_SEQUENCE_AA_LENGTH * 3) : mAminoAcids.size() < MINIMUM_SEQUENCE_AA_LENGTH; }
        inline bool translated() const { return !mAminoAcids.empty(); }
        inline bool aligned() const { return mAminoAcidsShift.aligned(); }
        inline bool matched() const { return !mHiNames.empty(); }

        inline bool has_lab(std::string aLab) const { return mLabIds.find(aLab) != mLabIds.end(); }
        inline std::string lab() const { return mLabIds.empty() ? std::string() : mLabIds.begin()->first; }
        inline std::string lab_id() const { return mLabIds.empty() ? std::string() : (mLabIds.begin()->second.empty() ? std::string() : mLabIds.begin()->second[0]); }
        inline const std::vector<std::string> cdcids() const { auto i = mLabIds.find("CDC"); return i == mLabIds.end() ? std::vector<std::string>{} : i->second; }
        inline const std::vector<std::string> lab_ids_for_lab(std::string lab) const { auto i = mLabIds.find(lab); return i == mLabIds.end() ? std::vector<std::string>{} : i->second; }
        inline const std::vector<std::string> lab_ids() const { std::vector<std::string> r; for (const auto& lid: mLabIds) { for (const auto& id: lid.second) { r.emplace_back(lid.first + "#" + id); } } return r; }
        inline const LabIds& lab_ids_raw() const { return mLabIds; }
        inline LabIds& lab_ids_raw() { return mLabIds; }
        inline bool match_labid(std::string lab, std::string id) const { auto i = mLabIds.find(lab); return i != mLabIds.end() && std::find(i->second.begin(), i->second.end(), id) != i->second.end(); }
        inline const auto& passages() const { return mPassages; }
        inline auto& passages() { return mPassages; }
        inline std::string passage() const { return mPassages.empty() ? std::string() : mPassages[0]; }
        inline bool passage_present(std::string aPassage) const { return mPassages.empty() ? aPassage.empty() : std::find(mPassages.begin(), mPassages.end(), aPassage) != mPassages.end(); }
        inline const auto& reassortant() const { return mReassortant; }
        inline auto& reassortant() { return mReassortant; }
        inline bool reassortant_match(std::string aReassortant) const { return mReassortant.empty() ? aReassortant.empty() : std::find(mReassortant.begin(), mReassortant.end(), aReassortant) != mReassortant.end(); }
        inline std::string gene() const { return mGene; }
        inline void gene(const char* str, size_t length) { mGene.assign(str, length); }

        inline const std::vector<std::string>& hi_names() const { return mHiNames; }
        inline std::vector<std::string>& hi_names() { return mHiNames; }
        inline void add_hi_name(std::string aHiName) { mHiNames.push_back(aHiName); }
        inline bool hi_name_present(std::string aHiName) const { return std::find(mHiNames.begin(), mHiNames.end(), aHiName) != mHiNames.end(); }

          // if aAligned && aLeftPartSize > 0 - include signal peptide and other stuff to the left from the beginning of the aligned sequence
        std::string amino_acids(bool aAligned, size_t aLeftPartSize = 0) const;
        std::string nucleotides(bool aAligned, size_t aLeftPartSize = 0) const;
        inline Shift amino_acids_shift() const { return mAminoAcidsShift; } // throws if sequence was not aligned
        inline Shift nucleotides_shift() const { return mNucleotidesShift; }  // throws if sequence was not aligned
        // inline int& amino_acids_shift_raw() { return mAminoAcidsShift.raw(); }
        // inline int& nucleotides_shift_raw() { return mNucleotidesShift.raw(); }
        inline void amino_acids_shift_raw(int shift) { mAminoAcidsShift.raw() = shift; }
        inline void nucleotides_shift_raw(int shift) { mNucleotidesShift.raw() = shift; }
        char amino_acid_at(size_t aPos) const; // aPos counts from 1!

        void amino_acids(const char* str, size_t length) { mAminoAcids.assign(str, length); }
        void nucleotides(const char* str, size_t length) { mNucleotides.assign(str, length); }

        std::vector<std::string> make_all_reassortant_passage_variants() const;

          // amino_acid_pos counts from 0!
        void add_deletions(size_t amino_acid_pos, size_t num_amino_acid_deletions);

          //   // Empty passages must not be removed! this is just for testing purposes
          // inline void remove_empty_passages()
          //     {
          //         mPassages.erase(std::remove(mPassages.begin(), mPassages.end(), std::string()), mPassages.end());
          //     }

     private:
        std::vector<std::string> mPassages;
        std::string mNucleotides;
        std::string mAminoAcids;
        Shift mNucleotidesShift;
        Shift mAminoAcidsShift;
        LabIds mLabIds;
        std::string mGene;
        std::vector<std::string> mHiNames;
        std::vector<std::string> mReassortant;
        std::vector<std::string> mClades;

        static inline std::string shift(std::string aSource, int aShift, char aFill)
            {
                std::string r = aSource;
                if (aShift < 0)
                    r.erase(0, static_cast<size_t>(-aShift));
                else if (aShift > 0)
                    r.insert(0, static_cast<size_t>(aShift), aFill);
                return r;
            }

        friend class Seqdb;
        friend class SeqdbIterator;
        friend class SeqdbIteratorBase;

    }; // class SeqdbSeq

// ----------------------------------------------------------------------

    class SeqdbEntry
    {
     public:
        inline SeqdbEntry() {}
        inline SeqdbEntry(std::string aName) : mName(aName) {}
        inline SeqdbEntry(std::string aName, std::string aVirusType, std::string aLineage) : mName(aName), mVirusType(aVirusType), mLineage(aLineage) {}

        inline const std::string name() const { return mName; }
        inline void name(const char* str, size_t length) { mName.assign(str, length); }
        inline const std::string country() const { return mCountry; }
        inline void country(std::string aCountry) { mCountry = aCountry; }
        inline void country(const char* str, size_t length) { mCountry.assign(str, length); }
        inline const std::string continent() const { return mContinent; }
        inline void continent(std::string aContinent) { mContinent = aContinent; }
        inline void continent(const char* str, size_t length) { mContinent.assign(str, length); }
        inline bool empty() const { return mSeq.empty(); }

        inline const std::string virus_type() const { return mVirusType; }
        inline void virus_type(std::string aVirusType) { mVirusType = aVirusType; }
        inline void virus_type(const char* str, size_t length) { mVirusType.assign(str, length); }
        void add_date(std::string aDate);
        inline const auto& dates() const { return mDates; }
        inline auto& dates() { return mDates; }
        inline const std::string date() const { return mDates.empty() ? std::string() : mDates.back(); }
        inline const std::string lineage() const { return mLineage; }
        inline void lineage(std::string aLineage) { mLineage = aLineage; }
        inline void lineage(const char* str, size_t length) { mLineage.assign(str, length); }
        void update_lineage(std::string aLineage, Messages& aMessages);
        void update_subtype_name(std::string aSubtype, Messages& aMessages);
          // returns warning message or an empty string
          // std::string add_or_update_sequence(std::string aSequence, std::string aPassage, std::string aReassortant, std::string aLab, std::string aLabId, std::string aGene);

        inline bool date_within_range(std::string aBegin, std::string aEnd) const
            {
                const std::string date = mDates.size() > 0 ? mDates.back() : "0000-00-00";
                return (aBegin.empty() || date >= aBegin) && (aEnd.empty() || date < aEnd);
            }

        inline void remove_short_sequences()
            {
                auto short_seq = [this](auto& seq) {
                    if (seq.is_short())
                        std::cerr << "Warning: removing too short sequence in " << mName << std::endl;
                    return seq.is_short();
                };
                mSeq.erase(std::remove_if(mSeq.begin(), mSeq.end(), short_seq), mSeq.end());
            }

        inline void remove_not_translated_sequences()
            {
                auto not_translated = [this](auto& seq) {
                    if (!seq.translated())
                        std::cerr << "Warning: removing not translated sequence in " << mName << std::endl;
                    return !seq.translated();
                };
                mSeq.erase(std::remove_if(mSeq.begin(), mSeq.end(), not_translated), mSeq.end());
            }

        inline std::vector<std::string> cdcids() const
            {
                std::vector<std::string> r;
                std::for_each(mSeq.begin(), mSeq.end(), [&r](auto const & seq) {auto seq_cdcids = seq.cdcids(); r.insert(r.end(), std::make_move_iterator(seq_cdcids.begin()), std::make_move_iterator(seq_cdcids.end())); });
                std::sort(r.begin(), r.end());
                r.erase(std::unique(r.begin(), r.end()), r.end());
                return r;
            }

        std::vector<std::string> make_all_names() const;
        std::vector<std::string> make_all_variants() const;

        inline const auto& seqs() const { return mSeq; }
        inline auto& seqs() { return mSeq; }
        inline auto begin_seq() { return mSeq.begin(); }
        inline auto end_seq() { return mSeq.end(); }
        inline auto begin_seq() const { return mSeq.begin(); }
        inline auto end_seq() const { return mSeq.end(); }
        inline size_t number_of_seqs() const { return seqs().size(); }

        const SeqdbSeq* find_by_hi_name(std::string aHiName) const;

          //   // Empty passages must not be removed! this is just for testing purposes
          // inline void remove_empty_passages()
          //     {
          //         std::for_each(mSeq.begin(), mSeq.end(), std::mem_fn(&SeqdbSeq::remove_empty_passages));
          //     }

     private:
        std::string mName;
        std::string mVirusType;
        std::string mLineage;
        std::string mCountry;
        std::string mContinent;
        std::vector<std::string> mDates;
        std::vector<SeqdbSeq> mSeq;

        friend class Seqdb;
        friend class SeqdbIteratorBase;
        friend class SeqdbIterator;
        friend class ConstSeqdbIterator;

    }; // class SeqdbEntry

// ----------------------------------------------------------------------

    class SeqdbEntrySeq
    {
     public:
        inline SeqdbEntrySeq() : mEntry(nullptr), mSeq(nullptr) {}
        inline SeqdbEntrySeq(SeqdbEntry& aEntry, SeqdbSeq& aSeq) : mEntry(&aEntry), mSeq(&aSeq) {}
        inline SeqdbEntrySeq(const SeqdbEntry& aEntry, const SeqdbSeq& aSeq) : mEntry(const_cast<SeqdbEntry*>(&aEntry)), mSeq(const_cast<SeqdbSeq*>(&aSeq)) {}

        inline void assign(SeqdbEntrySeq&& aEntrySeq) { mEntry = aEntrySeq.mEntry; mSeq = aEntrySeq.mSeq; }
        inline void assign(SeqdbEntry& aEntry, SeqdbSeq& aSeq) { mEntry = &aEntry; mSeq = &aSeq; }
        inline void assign(const SeqdbEntry& aEntry, const SeqdbSeq& aSeq) { mEntry = const_cast<SeqdbEntry*>(&aEntry); mSeq = const_cast<SeqdbSeq*>(&aSeq); }
        inline void assign(SeqdbEntry* aEntry, SeqdbSeq* aSeq) { mEntry = aEntry; mSeq = aSeq; }
        inline void assign(const SeqdbEntry* aEntry, const SeqdbSeq* aSeq) { mEntry = const_cast<SeqdbEntry*>(aEntry); mSeq = const_cast<SeqdbSeq*>(aSeq); }

        inline operator bool() const { return mEntry != nullptr && mSeq != nullptr; }

        inline SeqdbEntry& entry() { return *mEntry; }
        inline SeqdbSeq& seq() { return *mSeq; }
        inline const SeqdbEntry& entry() const { return *mEntry; }
        inline const SeqdbSeq& seq() const { return *mSeq; }

        inline std::string make_name(std::string aPassageSeparator = " ") const
            {
                return mEntry && mSeq ? (mSeq->hi_names().empty() ? string::strip(mEntry->name() + aPassageSeparator + mSeq->passage()) : mSeq->hi_names()[0]) : "*NOT-FOUND*";
            }

          // seq_id is concatenation of sequence name and passage separeted by __
        inline std::string seq_id(bool encoded) const
            {
                std::string r = (mEntry && mSeq) ? (string::strip(mEntry->name() + "__" + mSeq->passage())) : "*NOT-FOUND*";
                if (encoded)
                    r = name_encode(r);
                return r;
            }

     private:
        SeqdbEntry* mEntry;
        SeqdbSeq* mSeq;

    }; // class SeqdbEntrySeq

// ----------------------------------------------------------------------

    class SeqdbIteratorBase : public std::iterator<std::input_iterator_tag, SeqdbEntrySeq>
    {
     public:
        inline SeqdbIteratorBase(const SeqdbIteratorBase& a) = default;
        inline SeqdbIteratorBase(SeqdbIteratorBase&& a) = default;
        inline virtual ~SeqdbIteratorBase() {}

        inline virtual bool operator==(const SeqdbIteratorBase& aNother) const { return mEntryNo == aNother.mEntryNo && mSeqNo == aNother.mSeqNo; }
        inline virtual bool operator!=(const SeqdbIteratorBase& aNother) const { return ! operator==(aNother); }

        inline SeqdbIteratorBase& filter_lab(std::string aLab) { mLab = aLab; filter_added(); return *this; }
        inline SeqdbIteratorBase& filter_labid(std::string aLab, std::string aId) { mLabId = std::make_pair(aLab, aId); filter_added(); return *this; }
        inline SeqdbIteratorBase& filter_subtype(std::string aSubtype) { mSubtype = aSubtype; filter_added(); return *this; }
        inline SeqdbIteratorBase& filter_lineage(std::string aLineage) { mLineage = aLineage; filter_added(); return *this; }
        inline SeqdbIteratorBase& filter_aligned(bool aAligned) { mAligned = aAligned; filter_added(); return *this; }
        inline SeqdbIteratorBase& filter_gene(std::string aGene) { mGene = aGene; filter_added(); return *this; }
        inline SeqdbIteratorBase& filter_date_range(std::string aBegin, std::string aEnd) { mBegin = aBegin; mEnd = aEnd; filter_added(); return *this; }
        inline SeqdbIteratorBase& filter_hi_name(bool aHasHiName) { mHasHiName = aHasHiName; filter_added(); return *this; }
        inline SeqdbIteratorBase& filter_name_regex(std::string aNameRegex) { mNameMatcher.assign(aNameRegex, std::regex::icase); mNameMatcherSet = true; filter_added(); return *this; }

        virtual const Seqdb& seqdb() const = 0;
        virtual std::string make_name(std::string aPassageSeparator = " ") const = 0;

        inline SeqdbIteratorBase& operator ++ ();

        inline void validate() const;

     protected:
        inline SeqdbIteratorBase() : mNameMatcherSet(false) { end(); }
        inline SeqdbIteratorBase(size_t aEntryNo, size_t aSeqNo) : mEntryNo(aEntryNo), mSeqNo(aSeqNo), mAligned(false), mHasHiName(false), mNameMatcherSet(false) /*, mNameMatcher(".")*/ {}

        inline bool suitable_entry() const;
        inline bool suitable_seq() const;
        inline bool next_seq();
        inline void next_entry();

        inline size_t entry_no() const { return mEntryNo; }
        inline size_t seq_no() const { return mSeqNo; }

     private:
        size_t mEntryNo;
        size_t mSeqNo;

          // filter
        std::string mLab;
        std::string mSubtype;
        std::string mLineage;
        bool mAligned;
        std::string mGene;
        std::string mBegin;
        std::string mEnd;
        bool mHasHiName;
        bool mNameMatcherSet;
        std::regex mNameMatcher;
        std::pair<std::string, std::string> mLabId;

        inline void end() { mEntryNo = mSeqNo = std::numeric_limits<size_t>::max(); }
        inline void filter_added() { if (!suitable_entry() || !suitable_seq()) operator ++(); }

    }; // class SeqdbIteratorBase

// ----------------------------------------------------------------------

    class SeqdbIterator : public SeqdbIteratorBase
    {
     public:
        inline SeqdbIterator(const SeqdbIterator& a) = default;
        inline SeqdbIterator(SeqdbIterator&& a) = default;
        inline virtual ~SeqdbIterator() {}

        inline bool operator==(const SeqdbIterator& aNother) const { return &mSeqdb == &aNother.mSeqdb && SeqdbIteratorBase::operator==(aNother); }

        inline SeqdbEntrySeq operator*();

        inline virtual const Seqdb& seqdb() const { return mSeqdb; }
        inline virtual std::string make_name(std::string aPassageSeparator = " ") const { return const_cast<SeqdbIterator*>(this)->operator*().make_name(aPassageSeparator); }

     private:
        inline SeqdbIterator(Seqdb& aSeqdb) : SeqdbIteratorBase(), mSeqdb(aSeqdb) {}
        inline SeqdbIterator(Seqdb& aSeqdb, size_t aEntryNo, size_t aSeqNo) : SeqdbIteratorBase(aEntryNo, aSeqNo), mSeqdb(aSeqdb) {}

        Seqdb& mSeqdb;

        friend class Seqdb;

    }; // class SeqdbIterator

    class ConstSeqdbIterator : public SeqdbIteratorBase
    {
     public:
        inline ConstSeqdbIterator(const ConstSeqdbIterator& a) = default;
        inline ConstSeqdbIterator(ConstSeqdbIterator&& a) = default;
        inline virtual ~ConstSeqdbIterator() {}

        inline bool operator==(const ConstSeqdbIterator& aNother) const { return &mSeqdb == &aNother.mSeqdb && SeqdbIteratorBase::operator==(aNother); }

        inline const SeqdbEntrySeq operator*() const;

        inline virtual const Seqdb& seqdb() const { return mSeqdb; }
        inline virtual std::string make_name(std::string aPassageSeparator = " ") const { return operator*().make_name(aPassageSeparator); }

     private:
        inline ConstSeqdbIterator(const Seqdb& aSeqdb) : SeqdbIteratorBase(), mSeqdb(aSeqdb) {}
        inline ConstSeqdbIterator(const Seqdb& aSeqdb, size_t aEntryNo, size_t aSeqNo) : SeqdbIteratorBase(aEntryNo, aSeqNo), mSeqdb(aSeqdb) {}

        const Seqdb& mSeqdb;

        friend class Seqdb;

    }; // class ConstSeqdbIterator

// ----------------------------------------------------------------------

    class Seqdb
    {
     public:
        // inline Seqdb() = default;

        void load(std::string filename);
        void save(std::string filename, size_t indent = 0) const;

        inline size_t number_of_entries() const { return mEntries.size(); }
        inline size_t number_of_seqs() const { return std::accumulate(mEntries.begin(), mEntries.end(), 0U, [](size_t acc, const auto& e) { return acc + e.seqs().size(); }); }

        inline SeqdbEntry* find_by_name(std::string aName)
            {
                auto const first = find_insertion_place(aName);
                // if (first != mEntries.end() && aName != first->name())
                //     std::cerr << "Warining: looking for: \"" << aName << "\" found: \"" << first->name() << "\"" << std::endl;
                return (first != mEntries.end() && aName == first->name()) ? &(*first) : nullptr;
            }

        inline const SeqdbEntry* find_by_name(std::string aName) const
            {
                auto const first = find_insertion_place(aName);
                return (first != mEntries.end() && aName == first->name()) ? &(*first) : nullptr;
            }

        SeqdbEntrySeq find_by_seq_id(std::string aSeqId) const;

        // SeqdbEntry* new_entry(std::string aName);
        std::string add_sequence(std::string aName, std::string aVirusType, std::string aLineage, std::string aLab, std::string aDate, std::string aLabId, std::string aPassage, std::string aReassortant, std::string aSequence, std::string aGene);

          // fills by_virus_type that maps virus type to the list of indices of mEntries
        std::set<std::string> virus_types() const;
        void detect_insertions_deletions();
        void detect_b_lineage();

          // removes short sequences, removes entries having no sequences. returns messages
        std::string cleanup(bool remove_short_sequences);

          // returns db stat
        std::string report() const;
        std::string report_identical() const;
        std::string report_not_aligned(size_t prefix_size) const;
        std::vector<std::string> all_hi_names() const;
        std::vector<std::string> all_passages() const;
        void remove_hi_names();
        std::vector<std::string> match_hidb(bool aVerbose, bool aGreedy = true); // seqdb-hidb.cc  returns list of not found location names

          // iterating over sequences with filtering
        inline SeqdbIterator begin() { return SeqdbIterator(*this, 0, 0); }
        inline SeqdbIterator end() { return SeqdbIterator(*this); }
        inline ConstSeqdbIterator begin() const { return ConstSeqdbIterator(*this, 0, 0); }
        inline ConstSeqdbIterator end() const { return ConstSeqdbIterator(*this); }
        inline const auto& entries() const { return mEntries; }
        inline auto& entries() { return mEntries; }
        inline auto begin_entry() { return mEntries.begin(); }
        inline auto end_entry() { return mEntries.end(); }

        template <typename Value> std::deque<std::vector<SeqdbEntrySeq>> find_identical_sequences(Value value) const;

        void build_hi_name_index();
        const SeqdbEntrySeq* find_hi_name(std::string aHiName) const;

          // Matches antigens of a chart against seqdb, returns number of antigens matched.
          // Fills aPerAntigen with EntrySeq for each antigen.
        size_t match(const Antigens& aAntigens, std::vector<SeqdbEntrySeq>& aPerAntigen, bool aVerbose = true) const;
        std::vector<SeqdbEntrySeq> match(const Antigens& aAntigens, bool aVerbose = true) const;
          // Matches antigens of a chart against seqdb, for each matched antigen extract AA at the passed positions
        void aa_at_positions_for_antigens(const Antigens& aAntigens, const std::vector<size_t>& aPositions, std::map<std::string, std::vector<size_t>>& aa_indices, bool aVerbose) const;

     private:
        using HiNameIndex = std::map<std::string, SeqdbEntrySeq>;

        std::vector<SeqdbEntry> mEntries;
        const std::regex sReYearSpace = std::regex("/[12][0-9][0-9][0-9] ");
        HiNameIndex mHiNameIndex;

        inline std::vector<SeqdbEntry>::iterator find_insertion_place(std::string aName)
            {
                return std::lower_bound(mEntries.begin(), mEntries.end(), aName, [](const SeqdbEntry& entry, std::string name) -> bool { return entry.name() < name; });
            }

        inline std::vector<SeqdbEntry>::const_iterator find_insertion_place(std::string aName) const
            {
                return std::lower_bound(mEntries.begin(), mEntries.end(), aName, [](const SeqdbEntry& entry, std::string name) -> bool { return entry.name() < name; });
            }

        friend class SeqdbIteratorBase;
        friend class SeqdbIterator;
        friend class ConstSeqdbIterator;

        class LocationNotFound : public std::exception
        {
         public:
            inline LocationNotFound(std::string aLocation) : mLocation{aLocation} {}
            inline operator std::string() const { return mLocation; }
         private:
            std::string mLocation;
        };

          // throws LocationNotFound
        void find_in_hidb_update_country_lineage_date(std::vector<const hidb::AntigenData*>& found, SeqdbEntry& entry) const;
        // void split_by_virus_type(std::map<std::string, std::vector<size_t>>& by_virus_type) const;

    }; // class Seqdb

// ----------------------------------------------------------------------

    inline void SeqdbIteratorBase::validate() const
    {
        if (mEntryNo >= seqdb().mEntries.size() || mSeqNo >= seqdb().mEntries[mEntryNo].mSeq.size())
            throw std::out_of_range("SeqdbIterator is out of range");

    } // SeqdbIteratorBase::valid

// ----------------------------------------------------------------------

    inline SeqdbEntrySeq SeqdbIterator::operator*()
    {
        validate();
        return SeqdbEntrySeq(mSeqdb.mEntries[entry_no()], mSeqdb.mEntries[entry_no()].mSeq[seq_no()]);

    } // SeqdbIterator::operator*

// ----------------------------------------------------------------------

    inline const SeqdbEntrySeq ConstSeqdbIterator::operator*() const
    {
        validate();
        return SeqdbEntrySeq(mSeqdb.mEntries[entry_no()], mSeqdb.mEntries[entry_no()].mSeq[seq_no()]);

    } // SeqdbIterator::operator*

// ----------------------------------------------------------------------

// inline SeqdbIterator::ConstDereferenceType SeqdbIterator::operator*() const
// {
//     if (mEntryNo >= mSeqdb.mEntries.size() || mSeqNo >= mSeqdb.mEntries[mEntryNo].mSeq.size())
//         throw std::out_of_range("SeqdbIterator is out of range");
//     return std::make_pair(&mSeqdb.mEntries[mEntryNo], &mSeqdb.mEntries[mEntryNo].mSeq[mSeqNo]);

// } // SeqdbIterator::operator*

// ----------------------------------------------------------------------

    inline bool SeqdbIteratorBase::suitable_entry() const
    {
        auto const & entry = seqdb().mEntries[mEntryNo];
        return (mSubtype.empty() || entry.mVirusType == mSubtype)
                && (mLineage.empty() || entry.mLineage == mLineage)
                && entry.date_within_range(mBegin, mEnd)
                ;

    } // SeqdbIterator::suitable_entry

// ----------------------------------------------------------------------

    inline bool SeqdbIteratorBase::suitable_seq() const
    {
        auto const & seq = seqdb().mEntries[mEntryNo].mSeq[mSeqNo];
        return (!mAligned || seq.aligned())
                && (mGene.empty() || seq.mGene == mGene)
                && (!mHasHiName || !seq.mHiNames.empty())
                && (mLab.empty() || seq.has_lab(mLab))
                && (mLabId.first.empty() || seq.match_labid(mLabId.first, mLabId.second))
                && (!mNameMatcherSet || std::regex_search(make_name(), mNameMatcher))
                ;

    } // SeqdbIterator::suitable_seq

// ----------------------------------------------------------------------

    inline bool SeqdbIteratorBase::next_seq()
    {
        auto const & entry = seqdb().mEntries[mEntryNo];
        ++mSeqNo;
        while (mSeqNo < entry.mSeq.size() && !suitable_seq())
            ++mSeqNo;
        return mSeqNo < entry.mSeq.size();

    } // SeqdbIterator::next_seq

// ----------------------------------------------------------------------

    inline void SeqdbIteratorBase::next_entry()
    {
        while (true) {
            ++mEntryNo;
            while (mEntryNo < seqdb().mEntries.size() && !suitable_entry())
                ++mEntryNo;
            if (mEntryNo >= seqdb().mEntries.size()) {
                end();
                break;
            }
            else {
                mSeqNo = static_cast<size_t>(-1);
                if (next_seq())
                    break;
            }
        }

    } // SeqdbIterator::next_entry

// ----------------------------------------------------------------------

    inline SeqdbIteratorBase& SeqdbIteratorBase::operator ++ ()
    {
        if (!next_seq()) {
            next_entry();
        }
        return *this;

    } // SeqdbIterator::operator ++

// ----------------------------------------------------------------------

    template <typename Value> std::deque<std::vector<SeqdbEntrySeq>> Seqdb::find_identical_sequences(Value value) const
    {
        std::vector<SeqdbEntrySeq> refs(begin(), end());
        sort(refs.begin(), refs.end(), [&value](const auto& a, const auto b) { return value(a) < value(b); });
        std::deque<std::vector<SeqdbEntrySeq>> identical = {{}};
        for (auto previous = refs.begin(), current = previous + 1; current != refs.end(); ++current) {
            if (value(*previous) == value(*current)) {
                if (!value(*previous).empty()) { // empty means not aligned, ignore them
                    if (identical.back().empty())
                        identical.back().push_back(*previous);
                    identical.back().push_back(*current);
                }
            }
            else {
                previous = current;
                if (!identical.back().empty())
                    identical.push_back(std::vector<SeqdbEntrySeq>());
            }
        }
        if (identical.back().empty())
            identical.pop_back();
        return identical;

    } // Seqdb::find_identical_sequences

// ----------------------------------------------------------------------

    void setup(std::string aFilename);
    void setup_dbs(std::string aDbDir);
    const Seqdb& get(report_time aTimeit = report_time::No);

// ----------------------------------------------------------------------

} // namespace seqdb

// ----------------------------------------------------------------------

// template<> inline void swap(seqdb::SeqdbEntrySeq& a, seqdb::SeqdbEntrySeq& b) noexcept(::std::is_nothrow_move_constructible<seqdb::SeqdbEntrySeq>::value && ::std::is_nothrow_move_assignable<seqdb::SeqdbEntrySeq>::value)
// {
//     auto z = ::std::move(a);
//     a = ::std::move(b);
//     b = ::std::move(z);
// }

// ----------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& out, const seqdb::SeqdbEntry& entry)
{
    out << entry.virus_type() << " " << entry.name();
    for (auto& seq: entry.seqs()) {
        out << " {";
        if (!seq.reassortant().empty())
            out << 'R' << seq.reassortant().size() << seq.reassortant() << " ";
        out << seq.passages().size() << seq.passages() << " #" << seq.cdcids().size() << seq.cdcids() << '}';
    }
    return out;
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
