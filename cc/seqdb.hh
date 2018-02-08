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
#include "hidb-5/hidb.hh"
#include "seqdb/sequence-shift.hh"
#include "seqdb/amino-acids.hh"
#include "seqdb/messages.hh"

// ----------------------------------------------------------------------

namespace acmacs::chart { class Antigens; }

namespace seqdb
{
    class Seqdb;
    class SeqdbIterator;

    class import_error : public std::runtime_error { public: using std::runtime_error::runtime_error; };

// ----------------------------------------------------------------------

    class SequenceNotAligned : public std::runtime_error
    {
     public:
        SequenceNotAligned(std::string prefix) : std::runtime_error(prefix + ": sequence not aligned") {}
    };

// ----------------------------------------------------------------------

    class SeqdbSeq
    {
     public:
        using LabIds = std::map<std::string, std::vector<std::string>>;

        SeqdbSeq() : mGene("HA") {}

        SeqdbSeq(const std::string& aSequence, const std::string& aGene)
            : SeqdbSeq()
            {
                if (is_nucleotides(aSequence))
                    mNucleotides = aSequence;
                else
                    mAminoAcids = aSequence;
                if (!aGene.empty())
                    mGene = aGene;
            }

        // SeqdbSeq(bool aNucs, const std::string& aSequence, const std::string& aGene)
        //     : SeqdbSeq()
        //     {
        //         if (aNucs)
        //             mNucleotides = aSequence;
        //         else
        //             mAminoAcids = aSequence;
        //         if (!aGene.empty())
        //             mGene = aGene;
        //     }

        // SeqdbSeq(const std::string& aNucleotides, const std::string& aAminoAcids, const std::string& aGene)
        //     : SeqdbSeq()
        //     {
        //         mNucleotides = aNucleotides;
        //         mAminoAcids = aAminoAcids;
        //         if (!aGene.empty())
        //             mGene = aGene;
        //     }

        // SeqdbSeq(const std::string& aNucleotides, const std::string& aGene)
        //     : SeqdbSeq(aNucleotides, std::string(), aGene)
        //     {
        //     }

        AlignAminoAcidsData align(bool aForce, Messages& aMessages);

          // returns if aNucleotides matches mNucleotides or aAminoAcids matches mAminoAcids
        bool match_update(const SeqdbSeq& aNewSeq);
        bool match_update_nucleotides(const SeqdbSeq& aNewSeq);
        bool match_update_amino_acids(const SeqdbSeq& aNewSeq);
        void add_passage(const std::string& aPassage);
        void update_gene(const std::string& aGene, Messages& aMessages, bool replace_ha = false);
        void add_reassortant(const std::string& aReassortant);
        void add_lab_id(const std::string& aLab, const std::string& aLabId);
        const std::vector<std::string>& update_clades(const std::string& aVirusType, const std::string& aLineage);
        const std::vector<std::string>& clades() const { return mClades; }
        std::vector<std::string>& clades() { return mClades; }
        bool has_clade(const std::string& aClade) const { return std::find(std::begin(mClades), std::end(mClades), aClade) != std::end(mClades); }

        bool is_short() const { return mAminoAcids.empty() ? mNucleotides.size() < (MINIMUM_SEQUENCE_AA_LENGTH * 3) : mAminoAcids.size() < MINIMUM_SEQUENCE_AA_LENGTH; }
        bool translated() const { return !mAminoAcids.empty(); }
        bool aligned() const { return mAminoAcidsShift.aligned(); }
        bool matched() const { return !mHiNames.empty(); }

        bool has_lab(const std::string& aLab) const { return mLabIds.find(aLab) != mLabIds.end(); }
        std::string lab() const { return mLabIds.empty() ? std::string{} : mLabIds.begin()->first; }
        std::string lab_id() const { return mLabIds.empty() ? std::string{} : (mLabIds.begin()->second.empty() ? std::string{} : mLabIds.begin()->second[0]); }
        const std::vector<std::string> cdcids() const { auto i = mLabIds.find("CDC"); return i == mLabIds.end() ? std::vector<std::string>{} : i->second; }
        const std::vector<std::string> lab_ids_for_lab(const std::string& lab) const { auto i = mLabIds.find(lab); return i == mLabIds.end() ? std::vector<std::string>{} : i->second; }
        const std::vector<std::string> lab_ids() const { std::vector<std::string> r; for (const auto& lid: mLabIds) { for (const auto& id: lid.second) { r.emplace_back(lid.first + "#" + id); } } return r; }
        const LabIds& lab_ids_raw() const { return mLabIds; }
        LabIds& lab_ids_raw() { return mLabIds; }
        bool match_labid(const std::string& lab, const std::string& id) const { auto i = mLabIds.find(lab); return i != mLabIds.end() && std::find(i->second.begin(), i->second.end(), id) != i->second.end(); }
        const auto& passages() const { return mPassages; }
        auto& passages() { return mPassages; }
        std::string passage() const { return mPassages.empty() ? std::string{} : mPassages[0]; }
        bool passage_present(const std::string& aPassage) const { return mPassages.empty() ? aPassage.empty() : std::find(mPassages.begin(), mPassages.end(), aPassage) != mPassages.end(); }
        const auto& reassortant() const { return mReassortant; }
        auto& reassortant() { return mReassortant; }
        bool reassortant_match(const std::string& aReassortant) const { return mReassortant.empty() ? aReassortant.empty() : std::find(mReassortant.begin(), mReassortant.end(), aReassortant) != mReassortant.end(); }
        const std::string& gene() const { return mGene; }
        void gene(const char* str, size_t length) { mGene.assign(str, length); }

        const std::vector<std::string>& hi_names() const { return mHiNames; }
        std::vector<std::string>& hi_names() { return mHiNames; }
        void add_hi_name(const std::string& aHiName) { mHiNames.push_back(aHiName); }
        bool hi_name_present(const std::string& aHiName) const { return std::find(mHiNames.begin(), mHiNames.end(), aHiName) != mHiNames.end(); }

          // if aAligned && aLeftPartSize > 0 - include signal peptide and other stuff to the left from the beginning of the aligned sequence
        std::string amino_acids(bool aAligned, size_t aLeftPartSize = 0) const;
        std::string nucleotides(bool aAligned, size_t aLeftPartSize = 0) const;
        Shift amino_acids_shift() const { return mAminoAcidsShift; } // throws if sequence was not aligned
        Shift nucleotides_shift() const { return mNucleotidesShift; }  // throws if sequence was not aligned
        // int& amino_acids_shift_raw() { return mAminoAcidsShift.raw(); }
        // int& nucleotides_shift_raw() { return mNucleotidesShift.raw(); }
        void amino_acids_shift_raw(int shift) { mAminoAcidsShift.raw() = shift; }
        void nucleotides_shift_raw(int shift) { mNucleotidesShift.raw() = shift; }
        char amino_acid_at(size_t aPos) const; // aPos counts from 1!

        void amino_acids(const char* str, size_t length) { mAminoAcids.assign(str, length); }
        void nucleotides(const char* str, size_t length) { mNucleotides.assign(str, length); }

        std::vector<std::string> make_all_reassortant_passage_variants() const;

          // amino_acid_pos counts from 0!
        void add_deletions(size_t amino_acid_pos, size_t num_amino_acid_deletions);

          //   // Empty passages must not be removed! this is just for testing purposes
          // void remove_empty_passages()
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

        static inline std::string shift(const std::string& aSource, int aShift, char aFill)
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
        SeqdbEntry() {}
        SeqdbEntry(const std::string& aName) : mName(aName) {}
        SeqdbEntry(const std::string& aName, const std::string& aVirusType, const std::string& aLineage) : mName(aName), mVirusType(aVirusType), mLineage(aLineage) {}

        const std::string& name() const { return mName; }
        void name(const char* str, size_t length) { mName.assign(str, length); }
        const std::string& country() const { return mCountry; }
        void country(const std::string& aCountry) { mCountry = aCountry; }
        void country(const char* str, size_t length) { mCountry.assign(str, length); }
        const std::string& continent() const { return mContinent; }
        void continent(const std::string& aContinent) { mContinent = aContinent; }
        void continent(const char* str, size_t length) { mContinent.assign(str, length); }
        bool empty() const { return mSeq.empty(); }

        const std::string& virus_type() const { return mVirusType; }
        void virus_type(const std::string& aVirusType) { mVirusType = aVirusType; }
        void virus_type(const char* str, size_t length) { mVirusType.assign(str, length); }
        void add_date(const std::string& aDate);
        const auto& dates() const { return mDates; }
        auto& dates() { return mDates; }
        std::string date() const { return mDates.empty() ? std::string() : mDates.back(); }
        const std::string& lineage() const { return mLineage; }
        void lineage(const std::string& aLineage) { mLineage = aLineage; }
        void lineage(const char* str, size_t length) { mLineage.assign(str, length); }
        void update_lineage(const std::string& aLineage, Messages& aMessages);
        void update_subtype_name(const std::string& aSubtype, Messages& aMessages);
          // returns warning message or an empty string
          // std::string add_or_update_sequence(const std::string& aSequence, const std::string& aPassage, const std::string& aReassortant, const std::string& aLab, const std::string& aLabId, const std::string& aGene);

        bool date_within_range(const std::string& aBegin, const std::string& aEnd) const
            {
                const std::string date = mDates.size() > 0 ? mDates.back() : "0000-00-00";
                return (aBegin.empty() || date >= aBegin) && (aEnd.empty() || date < aEnd);
            }

        void remove_short_sequences()
            {
                auto short_seq = [this](auto& seq) {
                    if (seq.is_short())
                        std::cerr << "Warning: removing too short sequence in " << mName << std::endl;
                    return seq.is_short();
                };
                mSeq.erase(std::remove_if(mSeq.begin(), mSeq.end(), short_seq), mSeq.end());
            }

        void remove_not_translated_sequences()
            {
                auto not_translated = [this](auto& seq) {
                    if (!seq.translated())
                        std::cerr << "Warning: removing not translated sequence in " << mName << std::endl;
                    return !seq.translated();
                };
                mSeq.erase(std::remove_if(mSeq.begin(), mSeq.end(), not_translated), mSeq.end());
            }

        std::vector<std::string> cdcids() const
            {
                std::vector<std::string> r;
                std::for_each(mSeq.begin(), mSeq.end(), [&r](auto const & seq) {auto seq_cdcids = seq.cdcids(); r.insert(r.end(), std::make_move_iterator(seq_cdcids.begin()), std::make_move_iterator(seq_cdcids.end())); });
                std::sort(r.begin(), r.end());
                r.erase(std::unique(r.begin(), r.end()), r.end());
                return r;
            }

        std::vector<std::string> lab_ids() const
            {
                std::vector<std::string> r;
                std::for_each(mSeq.begin(), mSeq.end(), [&r](auto const & seq) {auto seq_lab_ids = seq.lab_ids(); r.insert(r.end(), std::make_move_iterator(seq_lab_ids.begin()), std::make_move_iterator(seq_lab_ids.end())); });
                std::sort(r.begin(), r.end());
                r.erase(std::unique(r.begin(), r.end()), r.end());
                return r;
            }

        std::vector<std::string> make_all_names() const;
        std::vector<std::string> make_all_variants() const;

        const auto& seqs() const { return mSeq; }
        auto& seqs() { return mSeq; }
        auto begin_seq() { return mSeq.begin(); }
        auto end_seq() { return mSeq.end(); }
        auto begin_seq() const { return mSeq.begin(); }
        auto end_seq() const { return mSeq.end(); }
        size_t number_of_seqs() const { return seqs().size(); }

        const SeqdbSeq* find_by_hi_name(const std::string& aHiName) const;

          //   // Empty passages must not be removed! this is just for testing purposes
          // void remove_empty_passages()
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
        SeqdbEntrySeq() : mEntry(nullptr), mSeq(nullptr) {}
        SeqdbEntrySeq(SeqdbEntry& aEntry, SeqdbSeq& aSeq) : mEntry(&aEntry), mSeq(&aSeq) {}
        SeqdbEntrySeq(const SeqdbEntry& aEntry, const SeqdbSeq& aSeq) : mEntry(const_cast<SeqdbEntry*>(&aEntry)), mSeq(const_cast<SeqdbSeq*>(&aSeq)) {}

        void assign(SeqdbEntrySeq&& aEntrySeq) { mEntry = aEntrySeq.mEntry; mSeq = aEntrySeq.mSeq; }
        void assign(SeqdbEntry& aEntry, SeqdbSeq& aSeq) { mEntry = &aEntry; mSeq = &aSeq; }
        void assign(const SeqdbEntry& aEntry, const SeqdbSeq& aSeq) { mEntry = const_cast<SeqdbEntry*>(&aEntry); mSeq = const_cast<SeqdbSeq*>(&aSeq); }
        void assign(SeqdbEntry* aEntry, SeqdbSeq* aSeq) { mEntry = aEntry; mSeq = aSeq; }
        void assign(const SeqdbEntry* aEntry, const SeqdbSeq* aSeq) { mEntry = const_cast<SeqdbEntry*>(aEntry); mSeq = const_cast<SeqdbSeq*>(aSeq); }

        operator bool() const { return mEntry != nullptr && mSeq != nullptr; }

        SeqdbEntry& entry() { return *mEntry; }
        SeqdbSeq& seq() { return *mSeq; }
        const SeqdbEntry& entry() const { return *mEntry; }
        const SeqdbSeq& seq() const { return *mSeq; }

        std::string make_name(const std::string& aPassageSeparator = " ") const
            {
                return mEntry && mSeq ? (mSeq->hi_names().empty() ? string::strip(mEntry->name() + aPassageSeparator + mSeq->passage()) : mSeq->hi_names()[0]) : "*NOT-FOUND*";
            }

        enum class encoded_t { no, yes };
          // seq_id is concatenation of sequence name and passage separeted by __
        std::string seq_id(encoded_t encoded) const
            {
                std::string r = (mEntry && mSeq) ? (string::strip(mEntry->name() + "__" + mSeq->passage())) : "*NOT-FOUND*";
                if (encoded == encoded_t::yes)
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
        SeqdbIteratorBase(const SeqdbIteratorBase& a) = default;
        SeqdbIteratorBase(SeqdbIteratorBase&& a) = default;
        virtual ~SeqdbIteratorBase() {}

        virtual bool operator==(const SeqdbIteratorBase& aNother) const { return mEntryNo == aNother.mEntryNo && mSeqNo == aNother.mSeqNo; }
        virtual bool operator!=(const SeqdbIteratorBase& aNother) const { return ! operator==(aNother); }

        SeqdbIteratorBase& filter_lab(const std::string& aLab) { mLab = aLab; filter_added(); return *this; }
        SeqdbIteratorBase& filter_labid(const std::string& aLab, const std::string& aId) { mLabId = std::make_pair(aLab, aId); filter_added(); return *this; }
        SeqdbIteratorBase& filter_subtype(const std::string& aSubtype) { mSubtype = aSubtype; filter_added(); return *this; }
        SeqdbIteratorBase& filter_lineage(const std::string& aLineage) { mLineage = aLineage; filter_added(); return *this; }
        SeqdbIteratorBase& filter_aligned(bool aAligned) { mAligned = aAligned; filter_added(); return *this; }
        SeqdbIteratorBase& filter_gene(const std::string& aGene) { mGene = aGene; filter_added(); return *this; }
        SeqdbIteratorBase& filter_date_range(const std::string& aBegin, const std::string& aEnd) { mBegin = aBegin; mEnd = aEnd; filter_added(); return *this; }
        SeqdbIteratorBase& filter_hi_name(bool aHasHiName) { mHasHiName = aHasHiName; filter_added(); return *this; }
        SeqdbIteratorBase& filter_name_regex(const std::string& aNameRegex) { mNameMatcher.assign(aNameRegex, std::regex::icase); mNameMatcherSet = true; filter_added(); return *this; }

        virtual const Seqdb& seqdb() const = 0;
        virtual std::string make_name(const std::string& aPassageSeparator = " ") const = 0;

        SeqdbIteratorBase& operator ++ ();

        void validate() const;

     protected:
        SeqdbIteratorBase() : mNameMatcherSet(false) { end(); }
        SeqdbIteratorBase(size_t aEntryNo, size_t aSeqNo) : mEntryNo(aEntryNo), mSeqNo(aSeqNo), mAligned(false), mHasHiName(false), mNameMatcherSet(false) /*, mNameMatcher(".")*/ {}

        bool suitable_entry() const;
        bool suitable_seq() const;
        bool next_seq();
        void next_entry();

        size_t entry_no() const { return mEntryNo; }
        size_t seq_no() const { return mSeqNo; }

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

        void end() { mEntryNo = mSeqNo = std::numeric_limits<size_t>::max(); }
        void filter_added() { if (!suitable_entry() || !suitable_seq()) operator ++(); }

    }; // class SeqdbIteratorBase

// ----------------------------------------------------------------------

    class SeqdbIterator : public SeqdbIteratorBase
    {
     public:
        SeqdbIterator(const SeqdbIterator& a) = default;
        SeqdbIterator(SeqdbIterator&& a) = default;
        virtual ~SeqdbIterator() {}

        bool operator==(const SeqdbIterator& aNother) const { return &mSeqdb == &aNother.mSeqdb && SeqdbIteratorBase::operator==(aNother); }

        SeqdbEntrySeq operator*();

        virtual const Seqdb& seqdb() const { return mSeqdb; }
        virtual std::string make_name(const std::string& aPassageSeparator = " ") const { return const_cast<SeqdbIterator*>(this)->operator*().make_name(aPassageSeparator); }

     private:
        SeqdbIterator(Seqdb& aSeqdb) : SeqdbIteratorBase(), mSeqdb(aSeqdb) {}
        SeqdbIterator(Seqdb& aSeqdb, size_t aEntryNo, size_t aSeqNo) : SeqdbIteratorBase(aEntryNo, aSeqNo), mSeqdb(aSeqdb) {}

        Seqdb& mSeqdb;

        friend class Seqdb;

    }; // class SeqdbIterator

    class ConstSeqdbIterator : public SeqdbIteratorBase
    {
     public:
        ConstSeqdbIterator(const ConstSeqdbIterator& a) = default;
        ConstSeqdbIterator(ConstSeqdbIterator&& a) = default;
        virtual ~ConstSeqdbIterator() {}

        bool operator==(const ConstSeqdbIterator& aNother) const { return &mSeqdb == &aNother.mSeqdb && SeqdbIteratorBase::operator==(aNother); }

        const SeqdbEntrySeq operator*() const;

        virtual const Seqdb& seqdb() const { return mSeqdb; }
        virtual std::string make_name(const std::string& aPassageSeparator = " ") const { return operator*().make_name(aPassageSeparator); }

     private:
        ConstSeqdbIterator(const Seqdb& aSeqdb) : SeqdbIteratorBase(), mSeqdb(aSeqdb) {}
        ConstSeqdbIterator(const Seqdb& aSeqdb, size_t aEntryNo, size_t aSeqNo) : SeqdbIteratorBase(aEntryNo, aSeqNo), mSeqdb(aSeqdb) {}

        const Seqdb& mSeqdb;

        friend class Seqdb;

    }; // class ConstSeqdbIterator

// ----------------------------------------------------------------------

    class Seqdb
    {
     public:
        // Seqdb() = default;

        void load(const std::string& filename);
        void save(const std::string& filename = std::string{}, size_t indent = 0) const;

        size_t number_of_entries() const { return mEntries.size(); }
        size_t number_of_seqs() const { return std::accumulate(mEntries.begin(), mEntries.end(), 0U, [](size_t acc, const auto& e) { return acc + e.seqs().size(); }); }

        SeqdbEntry* find_by_name(const std::string& aName)
            {
                auto const first = find_insertion_place(aName);
                // if (first != mEntries.end() && aName != first->name())
                //     std::cerr << "Warining: looking for: \"" << aName << "\" found: \"" << first->name() << "\"" << std::endl;
                return (first != mEntries.end() && aName == first->name()) ? &(*first) : nullptr;
            }

        const SeqdbEntry* find_by_name(const std::string& aName) const
            {
                auto const first = find_insertion_place(aName);
                return (first != mEntries.end() && aName == first->name()) ? &(*first) : nullptr;
            }

        SeqdbEntrySeq find_by_seq_id(const std::string& aSeqId) const;

        // SeqdbEntry* new_entry(const std::string& aName);
        std::string add_sequence(const std::string& aName, const std::string& aVirusType, const std::string& aLineage, const std::string& aLab, const std::string& aDate, const std::string& aLabId, const std::string& aPassage, const std::string& aReassortant, const std::string& aSequence, const std::string& aGene);

          // fills by_virus_type that maps virus type to the list of indices of mEntries
        std::set<std::string> virus_types() const;
        void detect_insertions_deletions();
        void detect_b_lineage();
        void update_clades(bool aVerbose);

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
        SeqdbIterator begin() { return SeqdbIterator(*this, 0, 0); }
        SeqdbIterator end() { return SeqdbIterator(*this); }
        ConstSeqdbIterator begin() const { return ConstSeqdbIterator(*this, 0, 0); }
        ConstSeqdbIterator end() const { return ConstSeqdbIterator(*this); }
        const auto& entries() const { return mEntries; }
        auto& entries() { return mEntries; }
        auto begin_entry() { return mEntries.begin(); }
        auto end_entry() { return mEntries.end(); }

        template <typename Value> std::deque<std::vector<SeqdbEntrySeq>> find_identical_sequences(Value value) const;

        void build_hi_name_index();
        const SeqdbEntrySeq* find_hi_name(const std::string& aHiName) const;

          // Matches antigens of a chart against seqdb, returns number of antigens matched.
          // Fills aPerAntigen with EntrySeq for each antigen.
        size_t match(const acmacs::chart::Antigens& aAntigens, std::vector<SeqdbEntrySeq>& aPerAntigen, const std::string& aChartVirusType, bool aVerbose = true) const;
        std::vector<SeqdbEntrySeq> match(const acmacs::chart::Antigens& aAntigens, const std::string& aChartVirusType, bool aVerbose = true) const;

          // Matches antigens of a chart against seqdb, for each matched antigen extract AA at the passed positions
        void aa_at_positions_for_antigens(const acmacs::chart::Antigens& aAntigens, const std::vector<size_t>& aPositions, std::map<std::string, std::vector<size_t>>& aa_indices, bool aVerbose) const;
        std::map<std::string, std::vector<size_t>> aa_at_positions_for_antigens(const acmacs::chart::Antigens& aAntigens, const std::vector<size_t>& aPositions, bool aVerbose) const
            {
                std::map<std::string, std::vector<size_t>> result;
                aa_at_positions_for_antigens(aAntigens, aPositions, result, aVerbose);
                return result;
            }

     private:
        using HiNameIndex = std::map<std::string, SeqdbEntrySeq>;

        std::vector<SeqdbEntry> mEntries;
        const std::regex sReYearSpace = std::regex("/[12][0-9][0-9][0-9] ");
        HiNameIndex mHiNameIndex;
        std::string mLoadedFromFilename;

        std::vector<SeqdbEntry>::iterator find_insertion_place(const std::string& aName)
            {
                return std::lower_bound(mEntries.begin(), mEntries.end(), aName, [](const SeqdbEntry& entry, const std::string& name) -> bool { return entry.name() < name; });
            }

        std::vector<SeqdbEntry>::const_iterator find_insertion_place(const std::string& aName) const
            {
                return std::lower_bound(mEntries.begin(), mEntries.end(), aName, [](const SeqdbEntry& entry, const std::string& name) -> bool { return entry.name() < name; });
            }

        friend class SeqdbIteratorBase;
        friend class SeqdbIterator;
        friend class ConstSeqdbIterator;

          // throws LocationNotFound
        void find_in_hidb_update_country_lineage_date(hidb::AntigenPList& found, SeqdbEntry& entry) const;
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

    void setup(const std::string& aFilename, bool aVerbose);
    void setup_dbs(const std::string& aDbDir, bool aVerbose);
    const Seqdb& get(report_time aTimeit = report_time::No);
    Seqdb& get_for_updating(report_time aTimeit = report_time::No);

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
