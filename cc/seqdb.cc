#include <typeinfo>
#include <tuple>
#include <set>

#include "acmacs-base/acmacsd.hh"
#include "acmacs-base/read-file.hh"
#include "acmacs-base/timeit.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/virus-name.hh"
#include "locationdb/locdb.hh"
#include "acmacs-chart-2/chart.hh"
#include "seqdb/seqdb.hh"
#include "clades.hh"
#include "seqdb-export.hh"
#include "seqdb-import.hh"
#include "insertions_deletions.hh"

using namespace seqdb;

// ----------------------------------------------------------------------

#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wglobal-constructors"
#pragma GCC diagnostic ignored "-Wexit-time-destructors"
#endif

static std::unique_ptr<Seqdb> sSeqdb;
static std::string sSeqdbFilename = acmacs::acmacsd_root() + "/data/seqdb.json.xz";
static bool sVerbose = false;

#pragma GCC diagnostic pop

void seqdb::setup(const std::string& aFilename, bool aVerbose)
{
    sVerbose = aVerbose;
    if (!aFilename.empty())
        sSeqdbFilename = aFilename;
}

const Seqdb& seqdb::get(report_time aTimeit)
{
    using namespace std::string_literals;
    if (!sSeqdb) {
        try {
            Timeit ti_seqdb{"DEBUG: SeqDb loading from " + sSeqdbFilename + ": ", sVerbose ? report_time::Yes : aTimeit};
            sSeqdb = std::make_unique<Seqdb>();
            sSeqdb->load(sSeqdbFilename);
            sSeqdb->build_hi_name_index();
        }
        catch (std::exception& err) {
            throw import_error("seqdb import: "s + err.what());
        }
    }
    return *sSeqdb;

} // seqdb::get

Seqdb& seqdb::get_for_updating(report_time aTimeit)
{
    return const_cast<Seqdb&>(get(aTimeit));

} // seqdb::get_for_updating

void seqdb::setup_dbs(const std::string& aDbDir, bool aVerbose)
{
    if (!aDbDir.empty()) {
        setup(aDbDir + "/seqdb.json.xz", aVerbose);
        locdb_setup(aDbDir + "/locationdb.json.xz", aVerbose);
    }
    else {
        setup(std::string{}, aVerbose);
        locdb_setup(std::string{}, aVerbose);
    }
    hidb::setup(aDbDir, {}, aVerbose);
}

// ----------------------------------------------------------------------

bool SeqdbSeq::match_update(const SeqdbSeq& aNewSeq)
{
    bool result = false;
    if (!aNewSeq.mNucleotides.empty())
        result = match_update_nucleotides(aNewSeq);
    else
        result = match_update_amino_acids(aNewSeq);
    return result;

} // SeqdbSeq::match_update

// ----------------------------------------------------------------------

bool SeqdbSeq::match_update_nucleotides(const SeqdbSeq& aNewSeq)
{
    bool matches = false;
    if (mNucleotides == aNewSeq.mNucleotides) {
        matches = true;
    }
    else if (mNucleotides.find(aNewSeq.mNucleotides) != std::string::npos) { // sub
        matches = true;
    }
    else if (aNewSeq.mNucleotides.find(mNucleotides) != std::string::npos) { // super
        matches = true;
        mNucleotides = aNewSeq.mNucleotides;
        mNucleotidesShift = aNewSeq.mNucleotidesShift;
        mAminoAcids = aNewSeq.mAminoAcids;
        mAminoAcidsShift = aNewSeq.mAminoAcidsShift;
    }
    return matches;

} // SeqdbSeq::match_update_nucleotides

// ----------------------------------------------------------------------

bool SeqdbSeq::match_update_amino_acids(const SeqdbSeq& aNewSeq)
{
    bool matches = false;
    if (mAminoAcids == aNewSeq.mAminoAcids) {
        matches = true;
    }
    else if (mAminoAcids.find(aNewSeq.mAminoAcids) != std::string::npos) { // sub
        matches = true;
    }
    else if (aNewSeq.mAminoAcids.find(mAminoAcids) != std::string::npos) { // super
        matches = true;
        mNucleotides = aNewSeq.mNucleotides;
        mNucleotidesShift = aNewSeq.mNucleotidesShift;
        mAminoAcids = aNewSeq.mAminoAcids;
        mAminoAcidsShift = aNewSeq.mAminoAcidsShift;
    }
    return matches;

} // SeqdbSeq::match_update_amino_acids

// ----------------------------------------------------------------------

void SeqdbSeq::add_passage(const std::string& aPassage)
{
    if (!aPassage.empty() && std::find(mPassages.begin(), mPassages.end(), aPassage) == mPassages.end())
        mPassages.push_back(aPassage);

} // SeqdbSeq::add_passage

// ----------------------------------------------------------------------

void SeqdbSeq::update_gene(const std::string& aGene, Messages& aMessages, bool replace_ha)
{
    if (!aGene.empty()) {
        if (mGene.empty())
            mGene = aGene;
        else if (aGene != mGene) {
            if (replace_ha && mGene == "HA")
                mGene = aGene;
            else
                aMessages.warning() << "[SAMESEQ] different genes " << mGene << " vs. " << aGene << '\n';
        }
    }

} // SeqdbSeq::update_gene

// ----------------------------------------------------------------------

void SeqdbSeq::add_reassortant(const std::string& aReassortant)
{
    if (!aReassortant.empty() && std::find(mReassortant.begin(), mReassortant.end(), aReassortant) == mReassortant.end())
        mReassortant.push_back(aReassortant);

} // SeqdbSeq::add_reassortant

// ----------------------------------------------------------------------

void SeqdbSeq::add_lab_id(const std::string& aLab, const std::string& aLabId)
{
    if (!aLab.empty()) {
        auto& lab_ids = mLabIds[aLab];
        if (!aLabId.empty() && std::find(lab_ids.begin(), lab_ids.end(), aLabId) == lab_ids.end()) {
            lab_ids.push_back(aLabId);
        }
    }

} // SeqdbSeq::add_lab_id

// ----------------------------------------------------------------------

AlignAminoAcidsData SeqdbSeq::align(bool aForce, Messages& aMessages)
{
    AlignAminoAcidsData align_data;

    enum WhatAlign {no_align, align_nucleotides, aling_amino_acids};
    WhatAlign what_align = no_align;
    if (!mNucleotides.empty() && (mAminoAcids.empty() || aForce))
        what_align = align_nucleotides;
    else if (!mAminoAcids.empty() && mNucleotides.empty() && (!mAminoAcidsShift.aligned() || aForce))
        what_align = aling_amino_acids;

    switch (what_align) {
      case no_align:
          break;
      case align_nucleotides:
          mAminoAcidsShift.reset();
          align_data = translate_and_align(mNucleotides, aMessages);
          if (!align_data.amino_acids.empty())
              mAminoAcids = align_data.amino_acids;
          if (!align_data.shift.alignment_failed()) {
              update_gene(align_data.gene, aMessages, true);
              if (align_data.shift.aligned()) {
                  mAminoAcidsShift = align_data.shift;
                  mNucleotidesShift = - align_data.offset + align_data.shift * 3;
              }
          }
          else {
              aMessages.warning() << "Nucs not translated/aligned" /* << mNucleotides */ << '\n';
          }
          break;
      case aling_amino_acids:
          mAminoAcidsShift.reset();
          align_data = align_amino_acids(mAminoAcids, aMessages);
          if (align_data.shift.aligned()) {
              mAminoAcidsShift = align_data.shift;
              update_gene(align_data.gene, aMessages, true);
          }
          else {
              aMessages.warning() << "AA not aligned" /* << mAminoAcids */ << '\n';
          }
          break;
    }

    return align_data;

} // SeqdbSeq::align

// ----------------------------------------------------------------------

const std::vector<std::string>& SeqdbSeq::update_clades(const std::string& aVirusType, const std::string& aLineage)
{
    if (aligned()) {
        if (aVirusType == "B") {
            if (aLineage == "YAMAGATA") {
                mClades = clades_b_yamagata(mAminoAcids, mAminoAcidsShift);
            }
            else if (aLineage == "VICTORIA") {
                mClades = clades_b_victoria(mAminoAcids, mAminoAcidsShift);
            }
        }
        else if (aVirusType == "A(H1N1)") {
            mClades = clades_h1pdm(mAminoAcids, mAminoAcidsShift);
        }
        else if (aVirusType == "A(H3N2)") {
            mClades = clades_h3n2(mAminoAcids, mAminoAcidsShift);
        }
        // else {
        //     std::cerr << "Cannot update clades for virus type " << aVirusType << '\n';
        // }
    }
    return mClades;

} // SeqdbSeq::update_clades

// ----------------------------------------------------------------------

std::string SeqdbSeq::amino_acids(bool aAligned, size_t aLeftPartSize) const
{
    std::string r = mAminoAcids;
    if (aAligned) {
        if (!aligned())
            throw SequenceNotAligned("amino_acids()");
        r = shift(r, mAminoAcidsShift + static_cast<int>(aLeftPartSize), 'X');

          // find the longest part not having *, replace parts before longest with X, truncate traling parts
        size_t longest_part_start = 0;
        size_t longest_part_len = 0;
        for (size_t start = 0; start != std::string::npos; ) {
            const size_t found = r.find_first_of('*', start);
            const size_t len = (found == std::string::npos ? r.size() : found) - start;
            if (longest_part_len < len) {
                longest_part_len = len;
                longest_part_start = start;
            }
            start = found == std::string::npos ? found : found + 1;
        }
        r.replace(0, longest_part_start, longest_part_start, 'X');
        r.resize(longest_part_start + longest_part_len, 'X');
    }
    return r;

} // SeqdbSeq::amino_acids

// ----------------------------------------------------------------------

// aPos counts from 1!
char SeqdbSeq::amino_acid_at(size_t aPos) const
{
    if (!aligned())
        throw SequenceNotAligned("amino_acid_at()");
    const int offset = static_cast<int>(aPos) - 1 - mAminoAcidsShift;
    if (offset < 0 || offset >= static_cast<int>(mAminoAcids.size()))
        throw std::runtime_error("amino_acid_at(): Invalid pos");
    return mAminoAcids[static_cast<size_t>(offset)];

} // SeqdbSeq::amino_acid_at

// ----------------------------------------------------------------------

std::string SeqdbSeq::nucleotides(bool aAligned, size_t aLeftPartSize) const
{
    std::string r = mNucleotides;
    if (aAligned) {
        if (!aligned())
            throw SequenceNotAligned("nucleotides()");
        r = shift(r, mNucleotidesShift + static_cast<int>(aLeftPartSize), '-');
    }
    return r;

} // SeqdbSeq::nucleotides

// ----------------------------------------------------------------------

std::vector<std::string> SeqdbSeq::make_all_reassortant_passage_variants() const
{
    std::vector<std::string> result;
    if (mReassortant.empty()) {
        result = mPassages;
    }
    else {
        for (const auto& reassortant: mReassortant) {
            std::transform(mPassages.begin(), mPassages.end(), std::back_inserter(result), [&reassortant](const auto& passage) -> std::string { return reassortant + " " + passage; });
        }
    }
    return result;

} // SeqdbSeq::make_all_reassortant_passage_variants

// ----------------------------------------------------------------------

  // amino_acid_pos counts from 0!
void SeqdbSeq::add_deletions(size_t amino_acid_pos, size_t num_amino_acid_deletions)
{
    const size_t aa_offset = static_cast<size_t>(static_cast<int>(amino_acid_pos) - amino_acids_shift());
    mAminoAcids.insert(aa_offset, num_amino_acid_deletions, '-');
    const size_t nuc_offset = static_cast<size_t>(static_cast<int>(amino_acid_pos * 3) - nucleotides_shift());
    mNucleotides.insert(nuc_offset, num_amino_acid_deletions * 3, '-');

} // SeqdbSeq::add_deletions

// ----------------------------------------------------------------------

void SeqdbEntry::add_date(const std::string& aDate)
{
    if (!aDate.empty()) {
        auto insertion_pos = std::lower_bound(mDates.begin(), mDates.end(), aDate);
        if (insertion_pos == mDates.end() || aDate != *insertion_pos) {
            mDates.insert(insertion_pos, aDate);
        }
    }

} // SeqdbEntry::add_date

// ----------------------------------------------------------------------

void SeqdbEntry::update_lineage(const std::string& aLineage, Messages& aMessages)
{
      // std::cerr << "Lineage " << mName << " " << (aLineage.empty() ? std::string("?") : aLineage) << '\n';
    if (!aLineage.empty()) {
        if (mLineage.empty())
            mLineage = aLineage;
        else if (aLineage != mLineage)
            aMessages.warning() << "Different lineages " << mLineage << " (stored) vs. " << aLineage << " (ignored)" << '\n';
    }

} // SeqdbEntry::update_lineage

// ----------------------------------------------------------------------

void SeqdbEntry::update_subtype_name(const std::string& aSubtype, Messages& aMessages)
{
    if (!aSubtype.empty()) {
        if (mVirusType.empty()) {
            mVirusType = aSubtype;
        }
        else if (aSubtype != mVirusType) {
            if (mVirusType == "A(H3N0)" && aSubtype == "A(H3N2)") {
                  // NIMR sent few sequences to gisaid having H3N0 while they are really H3N2 (detected by our aligner)
                mVirusType = aSubtype;
                  // replace virus type in the name too
                if (mName.find("A(H3N0)") == 0)
                    mName[5] = '2';
            }
            else {
                aMessages.warning() << "Different subtypes " << mVirusType << " (stored) vs. " << aSubtype << " (ignored)" << '\n';
            }
        }
          // fix subtype in the name too
        if (mName.find("A/") == 0)
            mName.replace(0, 1, mVirusType);
    }

} // SeqdbEntry::update_subtype_name

// ----------------------------------------------------------------------

const SeqdbSeq* SeqdbEntry::find_by_hi_name(const std::string& aHiName) const
{
    auto found = std::find_if(begin_seq(), end_seq(), [aHiName](const auto& seq) -> bool { return seq.hi_name_present(aHiName); });
    return found == end_seq() ? nullptr : &*found;

} // SeqdbEntry::find_by_hi_name

// ----------------------------------------------------------------------

// std::string SeqdbEntry::add_or_update_sequence(const std::string& aSequence, const std::string& aPassage, const std::string& aReassortant, const std::string& aLab, const std::string& aLabId, const std::string& aGene)
// {
//     Messages messages;
//     const bool nucs = is_nucleotides(aSequence);
//     decltype(mSeq.begin()) found;
//     if (nucs)
//         found = std::find_if(mSeq.begin(), mSeq.end(), [&aSequence](SeqdbSeq& seq) { return seq.match_update_nucleotides(aSequence); });
//     else
//         found = std::find_if(mSeq.begin(), mSeq.end(), [&aSequence](SeqdbSeq& seq) { return seq.match_update_amino_acids(aSequence); });
//     if (found != mSeq.end()) {  // update
//         found->update_gene(aGene, messages);
//     }
//     else { // add
//         if (nucs)
//             mSeq.push_back(SeqdbSeq(aSequence, aGene));
//         else
//             mSeq.push_back(SeqdbSeq(std::string(), aSequence, aGene));
//         found = mSeq.end() - 1;
//     }
//     if (found != mSeq.end()) {
//         found->add_passage(aPassage);
//         found->add_reassortant(aReassortant);
//         found->add_lab_id(aLab, aLabId);
//         auto align_data = found->align(false, messages);
//         update_subtype_name(align_data.subtype, messages);
//         update_lineage(align_data.lineage, messages);
//     }
//     return messages;

// } // SeqdbEntry::add_or_update_sequence

// ----------------------------------------------------------------------

std::vector<std::string> SeqdbEntry::make_all_names() const
{
    std::vector<std::string> result;
    for (const auto& seq: seqs()) {
        const auto variants = seq.make_all_reassortant_passage_variants();
        std::transform(variants.begin(), variants.end(), std::back_inserter(result), [this](const auto& var) -> std::string { return mName + " " + var; });
    }
    return result;

} // SeqdbEntry::make_all_names

// ----------------------------------------------------------------------

std::vector<std::string> SeqdbEntry::make_all_variants() const
{
    std::vector<std::string> result;
    for (const auto& seq: seqs()) {
        const auto variants = seq.make_all_reassortant_passage_variants();
        result.insert(result.end(), variants.begin(), variants.end());
    }
    return result;

} // SeqdbEntry::make_variants

// ----------------------------------------------------------------------

std::string Seqdb::add_sequence(const std::string& aName, const std::string& aVirusType, const std::string& aLineage, const std::string& aLab, const std::string& aDate, const std::string& aLabId, const std::string& aPassage, const std::string& aReassortant, const std::string& aSequence, const std::string& aGene)
{
    Messages messages;
    std::string name = virus_name::normalize(aName);
    SeqdbEntry entry(name, aVirusType, aLineage);

    SeqdbSeq new_seq(aSequence, aGene);
    auto align_data = new_seq.align(false, messages);
    entry.update_subtype_name(align_data.subtype, messages); // updates entry.mName!
    // std::cerr << "add " << align_data.subtype << ' ' << entry.name() << '\n';

    auto inserted_entry = find_insertion_place(entry.name());
    SeqdbSeq* inserted_seq = nullptr;
    if (inserted_entry == mEntries.end() || entry.name() != inserted_entry->name()) {
        inserted_entry = mEntries.insert(inserted_entry, std::move(entry));
        inserted_entry->seqs().push_back(std::move(new_seq));
        inserted_seq = &inserted_entry->seqs().back();
    }
    else {
        inserted_entry->update_subtype_name(align_data.subtype, messages);
        auto& seqs = inserted_entry->seqs();
        auto found_seq = std::find_if(seqs.begin(), seqs.end(), [&new_seq](SeqdbSeq& seq) { return seq.match_update(new_seq); });
        if (found_seq == seqs.end()) {
            seqs.push_back(std::move(new_seq));
            inserted_seq = &seqs.back();
        }
        else {
            inserted_seq = &*found_seq;
        }
    }
    inserted_entry->update_lineage(align_data.lineage, messages);
    if (!aDate.empty())
        inserted_entry->add_date(aDate);

    inserted_seq->add_reassortant(aReassortant);
    inserted_seq->add_passage(aPassage);
    inserted_seq->add_lab_id(aLab, aLabId);

    return messages;

} // Seqdb::add_sequence

// ----------------------------------------------------------------------

// SeqdbEntry* Seqdb::new_entry(const std::string& aName)
// {
//     auto const first = find_insertion_place(aName);
//     if (first != mEntries.end() && aName == first->name())
//         throw std::runtime_error(std::string("Entry for \"") + aName + "\" already exists");
//     auto inserted = mEntries.insert(first, SeqdbEntry(aName));
//     return &*inserted;

// } // Seqdb::new_entry

// ----------------------------------------------------------------------

std::string Seqdb::cleanup(bool remove_short_sequences)
{
    Messages messages;
    if (remove_short_sequences) {
        std::for_each(mEntries.begin(), mEntries.end(), std::mem_fn(&SeqdbEntry::remove_short_sequences));
    }

    std::for_each(mEntries.begin(), mEntries.end(), std::mem_fn(&SeqdbEntry::remove_not_translated_sequences));

      // Note empty passages must not be removed, this is just for testing purposes
      // std::for_each(mEntries.begin(), mEntries.end(), std::mem_fn(&SeqdbEntry::remove_empty_passages));

      // remove empty entries
    auto const num_entries_before = mEntries.size();
    mEntries.erase(std::remove_if(mEntries.begin(), mEntries.end(), std::mem_fn(&SeqdbEntry::empty)), mEntries.end());
    if (mEntries.size() != num_entries_before)
        messages.warning() << (num_entries_before - mEntries.size()) << " entries removed during cleanup" << '\n';
    return messages;

} // Seqdb::cleanup

// ----------------------------------------------------------------------

std::string Seqdb::report() const
{
    std::ostringstream os;
    os << "Entries: " << mEntries.size() << '\n';

    std::map<std::string, size_t> virus_types;
    for_each(mEntries.begin(), mEntries.end(), [&virus_types](auto& e) { ++virus_types[e.virus_type()]; });
    os << "Virus types: " << virus_types << '\n';

    std::map<std::string, size_t> lineages;
    for_each(mEntries.begin(), mEntries.end(), [&lineages](auto& e) { ++lineages[e.lineage()]; });
    os << "Lineages: " << lineages << '\n';

    std::map<std::string, size_t> aligned;
    for_each(mEntries.begin(), mEntries.end(), [&aligned](auto& e) { if (any_of(e.mSeq.begin(), e.mSeq.end(), std::mem_fn(&SeqdbSeq::aligned))) ++aligned[e.virus_type()]; });
    os << "Aligned: " << aligned << '\n';

    std::map<std::string, size_t> matched;
    for_each(mEntries.begin(), mEntries.end(), [&matched](auto& e) { if (any_of(e.mSeq.begin(), e.mSeq.end(), std::mem_fn(&SeqdbSeq::matched))) ++matched[e.virus_type()]; });
    os << "Matched: " << matched << '\n';

    std::map<std::string, size_t> have_dates;
    for_each(mEntries.begin(), mEntries.end(), [&have_dates](auto& e) { if (!e.mDates.empty()) ++have_dates[e.virus_type()]; });
    os << "Have dates: " << have_dates << '\n';

    std::map<std::string, size_t> have_clades;
    for_each(mEntries.begin(), mEntries.end(), [&have_clades](auto& e) { if (any_of(e.mSeq.begin(), e.mSeq.end(), [](auto& seq) { return !seq.mClades.empty(); })) ++have_clades[e.virus_type()]; });
    os << "Have clades: " << have_clades << '\n';

    return os.str();

} // Seqdb::report

// ----------------------------------------------------------------------

std::string Seqdb::report_identical() const
{
    std::ostringstream os;

    auto report = [&os](const std::string& prefix, const auto& groups) {
        if (!groups.empty()) {
            os << prefix << '\n';
            for (auto const& group: groups) {
                for (auto const& entry: group) {
                    os << entry.make_name() << ' ';
                }
                os << '\n';
            }
        }
    };

    try {
        report("Identical nucleotides:", find_identical_sequences([](const SeqdbEntrySeq& e) -> std::string { try { return e.seq().nucleotides(true); } catch (SequenceNotAligned&) { return std::string(); } }));
        os << '\n';
    }
    catch (std::exception& err) {
        os << "Cannot report identical nucleotides: " << typeid(err).name() << ": " << err.what() << '\n';
    }

    try {
        report("Identical amino-acids:", find_identical_sequences([](const SeqdbEntrySeq& e) -> std::string { try { return e.seq().amino_acids(true); } catch (SequenceNotAligned&) { return std::string(); } }));
        os << '\n';
    }
    catch (std::exception& err) {
        os << "Cannot report identical amino-acids: " << typeid(err).name() << ": " << err.what() << '\n';
    }

    return os.str();

} // Seqdb::report_identical

// ----------------------------------------------------------------------

std::string Seqdb::report_not_aligned(size_t prefix_size) const
{
    std::vector<SeqdbEntrySeq> not_aligned;
    std::copy_if(begin(), end(), std::back_inserter(not_aligned), [](const auto& e) -> bool { return !e.seq().aligned(); });
    std::vector<std::string> prefixes;
    std::transform(not_aligned.begin(), not_aligned.end(), std::back_inserter(prefixes), [&prefix_size](const auto& e) -> std::string { return std::string(e.seq().amino_acids(false), 0, prefix_size); });
    std::sort(prefixes.begin(), prefixes.end());
    const auto p_end = std::unique(prefixes.begin(), prefixes.end());

    std::ostringstream os;
    os << "Prefixes of not aligned sequences of length " << prefix_size << ": " << p_end - prefixes.begin() << '\n';
    std::copy(prefixes.begin(), p_end, std::ostream_iterator<std::string>(os, "\n"));

    return os.str();

} // Seqdb::report_not_aligned

// ----------------------------------------------------------------------

std::vector<std::string> Seqdb::all_hi_names() const
{
    std::vector<std::string> r;
    for (auto const& entry: mEntries) {
        for (auto const& seq: entry.mSeq) {
            std::copy(seq.hi_names().begin(), seq.hi_names().end(), std::back_inserter(r));
        }
    }
    std::sort(r.begin(), r.end());
    auto last = std::unique(r.begin(), r.end());
    r.erase(last, r.end());
    return r;

} // Seqdb::all_hi_names

// ----------------------------------------------------------------------

void Seqdb::remove_hi_names()
{
    for (auto& entry: mEntries) {
        for (auto& seq: entry.mSeq) {
            seq.hi_names().clear();
        }
    }

} // Seqdb::remove_hi_names

// ----------------------------------------------------------------------

std::vector<std::string> Seqdb::all_passages() const
{
    std::vector<std::string> passages;
    for (const auto& entry: mEntries) {
        for (const auto& seq: entry.seqs()) {
            std::copy(seq.passages().begin(), seq.passages().end(), std::back_inserter(passages));
        }
    }
    std::sort(passages.begin(), passages.end());
    passages.erase(std::unique(passages.begin(), passages.end()), passages.end());
    return passages;

} // Seqdb::all_passages

// ----------------------------------------------------------------------

void Seqdb::find_in_hidb_update_country_lineage_date(hidb::AntigenPList& found, SeqdbEntry& entry) const
{
    try {
        auto hidb_antigens = hidb::get(entry.virus_type()).antigens();
        if (const auto cdcids = entry.cdcids(); !cdcids.empty()) {
            for (const auto& cdcid: cdcids) {
                const auto f_cdcid = hidb_antigens->find_labid(cdcid);
                std::copy(f_cdcid.begin(), f_cdcid.end(), std::back_inserter(found));
            }
        }

        std::string name = entry.name();
        hidb::AntigenPIndexList antigen_index_list;
        try {
            antigen_index_list = hidb_antigens->find(name, hidb::FixLocation::Yes, hidb::FindFuzzy::No);
        }
        catch (LocationNotFound&) {
            if (std::count(name.begin(), name.end(), '/') == 4) {
                  // perhaps isolation split with / and there is no host
                  // in that case HI table usually has - instead of / (fixed manually by Eu on table parsing)
                const auto parts = acmacs::string::split(name, "/", acmacs::string::Split::KeepEmpty);
                name = string::join("/", {parts[0], parts[1], string::join("-", {parts[2], parts[3]}), parts[4]});
                antigen_index_list = hidb_antigens->find(name, hidb::FixLocation::Yes, hidb::FindFuzzy::No);
            }
            else {
                  // try without fixing location
                antigen_index_list = hidb_antigens->find(name, hidb::FixLocation::No, hidb::FindFuzzy::No);
            }
        }
        std::transform(antigen_index_list.begin(), antigen_index_list.end(), std::back_inserter(found), [](const hidb::AntigenPIndex& antigen_index) -> hidb::AntigenP { return antigen_index.first; });
        std::sort(found.begin(), found.end());
        found.erase(std::unique(found.begin(), found.end()), found.end());
    }
    catch (hidb::get_error& err) {
        std::cerr << "WARNING: no HiDb for " << entry.name() << ": " << err.what() << '\n';
    }

    if (!found.empty()) {
          // update country and continent
        if (entry.country().empty()) {
            try {
                const std::string country = get_locdb().find(virus_name::location(entry.name())).country();
                entry.country(country);
                entry.continent(get_locdb().continent_of_country(country));
            }
            catch (LocationNotFound&) {
            }
            catch (virus_name::Unrecognized&) {
            }
        }

          // update lineage
        if (entry.virus_type() == "B") {
            Messages messages;
            entry.update_lineage(found.front()->lineage(), messages);
            if (messages)
                std::cerr << messages << '\n';
        }

          // update date
        for (const auto& e: found) {
            entry.add_date(e->date());
        }
    }

} // Seqdb::find_in_hidb_update_country_lineage_date

// ----------------------------------------------------------------------

SeqdbEntrySeq Seqdb::find_by_seq_id(const std::string& aSeqId) const
{
    SeqdbEntrySeq result;
    const std::string seq_id = name_decode(aSeqId);
    auto passage_separator = seq_id.find("__");
    if (passage_separator != std::string::npos) { // seq_id
        if (const auto entry = find_by_name(std::string(seq_id, 0, passage_separator)); entry != nullptr) {
            const auto passage_distinct = acmacs::string::split(string::string_view(seq_id, passage_separator + 2), "__", acmacs::string::Split::KeepEmpty);
            auto index = passage_distinct.size() == 1 ? 0 : std::stoi(std::string(passage_distinct[1]));
            for (const auto& seq: entry->seqs()) {
                if (seq.passage() == passage_distinct[0]) {
                    if (index == 0) {
                        result.assign(*entry, seq);
                        break;
                    }
                    else
                        --index;
                }
            }
        }
        else {
            std::cerr << "Error: no entry for \"" << std::string(seq_id, 0, passage_separator) << "\" in seqdb [" << __FILE__ << ":" << __LINE__ << ']' << '\n';
        }
    }
    else {
        std::smatch year_space;
        const auto year_space_present = std::regex_search(seq_id, year_space, sReYearSpace);
        const std::string look_for = year_space_present ? std::string(seq_id, 0, static_cast<std::string::size_type>(year_space.position(0) + year_space.length(0)) - 1) : seq_id;
        if (const auto entry = find_by_name(look_for); entry != nullptr) {
            auto found = std::find_if(entry->begin_seq(), entry->end_seq(), [&seq_id](const auto& seq) -> bool { return seq.hi_name_present(seq_id); });
            if (found == entry->end_seq()) { // not found by hi_name, look by passage (or empty passage)
                const std::string passage = year_space_present ? std::string(seq_id, static_cast<std::string::size_type>(year_space.position(0) + year_space.length(0))) : std::string();
                found = std::find_if(entry->begin_seq(), entry->end_seq(), [&passage](const auto& seq) -> bool { return seq.passage_present(passage); });
            }
            if (found != entry->end_seq()) {
                result.assign(*entry, *found);
            }
        }
    }

    if (!result) {
        std::cerr << "Error: \"" << seq_id << "\" not in seqdb" << '\n';
    }
    return result;

} // Seqdb::find_by_seq_id

// ----------------------------------------------------------------------

void Seqdb::build_hi_name_index()
{
    mHiNameIndex.clear();
    for (auto entry_seq: *this) {
        for (auto hi_name: entry_seq.seq().hi_names()) {
              // const auto pos_inserted =
            mHiNameIndex.emplace(hi_name, entry_seq);
              // if (!pos_inserted.second)
              //     std::cerr << "warning:0: " << hi_name << " was already in mHiNameIndex [" << pos_inserted.first->second.entry().name() << "] [" << entry_seq.entry().name() << ']' << '\n';
        }
    }

} // Seqdb::build_hi_name_index

// ----------------------------------------------------------------------

const SeqdbEntrySeq* Seqdb::find_hi_name(const std::string& aHiName) const
{
    const auto it = mHiNameIndex.find(aHiName);
    return it == mHiNameIndex.end() ? nullptr : &it->second;

} // Seqdb::find_hi_name

// ----------------------------------------------------------------------

size_t Seqdb::match(const acmacs::chart::Antigens& aAntigens, std::vector<SeqdbEntrySeq>& aPerAntigen, const std::string& aChartVirusType, bool aVerbose) const
{
    size_t matched = 0;
    aPerAntigen.clear();
    for (auto antigen: aAntigens) {
        const SeqdbEntrySeq* entry = find_hi_name(antigen->full_name());
        if (!entry)
            entry = find_hi_name(antigen->full_name_for_seqdb_matching());
        if (entry) {
            if (!aChartVirusType.empty() && aChartVirusType != entry->entry().virus_type()) {
                if (aVerbose)
                    std::cerr << "WARNING: Seqdb::match: virus type mismatch: chart:" << aChartVirusType << " seq:" << entry->entry().virus_type() << " name: " << antigen->full_name() << '\n';
            }
            else if (antigen->lineage() != acmacs::chart::BLineage::Unknown && static_cast<std::string>(antigen->lineage()) != entry->entry().lineage()) {
                std::cerr << "WARNING: Seqdb::match: lineage mismatch: antigen:" << antigen->lineage() << " seq:" << entry->entry().lineage() << " name: " << antigen->full_name() << '\n';
            }
            else {
                aPerAntigen.push_back(*entry);
                ++matched;
            }
        }
        else {
            aPerAntigen.emplace_back();
            // if (aVerbose)
            //     std::cerr << "WARNING: seqdb::match failed for \"" << antigen->full_name() << "\"" << '\n';
        }
    }
    if (aVerbose)
        std::cerr << "INFO: " << matched << " antigens from chart have sequences in seqdb" << '\n';
    return matched;

} // Seqdb::match

// ----------------------------------------------------------------------

std::vector<SeqdbEntrySeq> Seqdb::match(const acmacs::chart::Antigens& aAntigens, const std::string& aChartVirusType, bool aVerbose) const
{
    std::vector<SeqdbEntrySeq> per_antigen;
    match(aAntigens, per_antigen, aChartVirusType, aVerbose);
    return per_antigen;

} // Seqdb::match

// ----------------------------------------------------------------------

void Seqdb::aa_at_positions_for_antigens(const acmacs::chart::Antigens& aAntigens, const std::vector<size_t>& aPositions, std::map<std::string, std::vector<size_t>>& aa_indices, bool aVerbose) const
{
    size_t matched = 0;
    for (auto ag = aAntigens.begin(); ag != aAntigens.end(); ++ag) {
        const SeqdbEntrySeq* entry = find_hi_name((*ag)->full_name());
        if (!entry)
            entry = find_hi_name((*ag)->full_name_for_seqdb_matching());
        if (entry) {
            std::string aa(aPositions.size(), 'X');
            std::transform(aPositions.begin(), aPositions.end(), aa.begin(), [&entry](size_t pos) { return entry->seq().amino_acid_at(pos); });
            aa_indices[aa].push_back(ag.index());
            ++matched;
        }
    }
    if (aVerbose)
        std::cerr << "INFO: " << matched << " antigens from chart have sequences in seqdb" << '\n';

} // Seqdb::aa_at_positions_for_antigens

// ----------------------------------------------------------------------

void Seqdb::load(const std::string& filename)
{
    seqdb_import(filename, *this);
    mLoadedFromFilename = filename;

} // Seqdb::from_json_file

// ----------------------------------------------------------------------

void Seqdb::save(const std::string& filename, size_t indent) const
{
    seqdb_export(filename.empty() ? mLoadedFromFilename : filename, *this, indent);

} // Seqdb::save

// ----------------------------------------------------------------------

std::set<std::string> Seqdb::virus_types() const
{
    std::set<std::string> result;
    std::transform(mEntries.begin(), mEntries.end(), std::inserter(result, result.end()), [](const auto& entry) { return entry.virus_type(); });
    return result;

} // Seqdb::virus_types

// ----------------------------------------------------------------------

// void Seqdb::split_by_virus_type(std::map<std::string, std::vector<size_t>>& by_virus_type) const
// {
//     for (size_t index = 0; index < mEntries.size(); ++index) {
//         auto ptr_inserted = by_virus_type.emplace(mEntries[index].virus_type(), std::vector<size_t>{});
//         ptr_inserted.first->second.push_back(index);
//     }

// } // Seqdb::split_by_virus_type

// ----------------------------------------------------------------------

void Seqdb::detect_insertions_deletions()
{
    for (std::string virus_type: virus_types()) {
        if (!virus_type.empty()) {
            std::cout << "Detect insertions/deletions for " << virus_type << '\n';
            InsertionsDeletionsDetector detector(*this, virus_type);
            detector.detect();
        }
    }

} // Seqdb::detect_insertions_deletions

// ----------------------------------------------------------------------

void Seqdb::update_clades(bool aVerbose)
{
    std::map<std::string, size_t> clade_count;
    for (auto entry_seq: *this) {
        const auto& clades = entry_seq.seq().update_clades(entry_seq.entry().virus_type(), entry_seq.entry().lineage());
        for (const auto& clade: clades)
            ++clade_count[clade];
    }
    if (aVerbose)
        std::cerr << "INFO: clades: " << clade_count << '\n';

} // Seqdb::update_clades

// ----------------------------------------------------------------------

void Seqdb::detect_b_lineage()
{
    BLineageDetector detector(*this);
    detector.detect();

} // Seqdb::detect_b_lineage

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
