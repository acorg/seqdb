#include <typeinfo>
#include <tuple>
#include <set>

#include "acmacs-base/acmacsd.hh"
#include "acmacs-base/read-file.hh"
#include "acmacs-base/timeit.hh"
#include "acmacs-base/string.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/to-json.hh"
#include "acmacs-base/enumerate.hh"
#include "locationdb/locdb.hh"
#include "acmacs-virus/virus-name.hh"
#include "acmacs-chart-2/chart-modify.hh"
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
static seqdb::report sReport = seqdb::report::no;

#pragma GCC diagnostic pop

void seqdb::setup(std::string_view aFilename, seqdb::report aReport)
{
    sReport = aReport;
    if (!aFilename.empty())
        sSeqdbFilename = aFilename;
}

const Seqdb& seqdb::get(ignore_errors ignore_err, report_time aTimeit)
{
    using namespace std::string_literals;
    if (!sSeqdb) {
        Timeit ti_seqdb{"DEBUG: SeqDb loading from " + sSeqdbFilename + ": ", sReport == report::yes ? report_time::yes : aTimeit};
            sSeqdb = std::make_unique<Seqdb>();
            try {
                sSeqdb->load(sSeqdbFilename);
                sSeqdb->build_hi_name_index();
            }
            catch (std::exception& err) {
                if (ignore_err == ignore_errors::no)
                    throw import_error("seqdb import: "s + err.what());
                else
                    std::cerr << "WARNING: seqdb import error (ignored): " << err.what() << '\n';
            }
    }
    return *sSeqdb;

} // seqdb::get

Seqdb& seqdb::get_for_updating(report_time aTimeit)
{
    return const_cast<Seqdb&>(get(ignore_errors::no, aTimeit));

} // seqdb::get_for_updating

void seqdb::setup_dbs(std::string_view aDbDir, seqdb::report aReport)
{
    if (!aDbDir.empty()) {
        setup(string::concat(aDbDir, "/seqdb.json.xz"), aReport);
        locdb_setup(string::concat(aDbDir, "/locationdb.json.xz"), aReport == report::yes ? true : false);
    }
    else {
        setup(std::string{}, aReport);
        locdb_setup(std::string{}, aReport == report::yes ? true : false);
    }
    hidb::setup(aDbDir, {}, aReport == report::yes ? true : false);
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

void SeqdbSeq::add_passage(std::string aPassage)
{
    if (!aPassage.empty()) {
        if (aPassage.substr(0, 8) == "PASSAGE-") // artefact in few (400+) names in gisaid
            aPassage = aPassage.substr(8);
        if (std::find(mPassages.begin(), mPassages.end(), aPassage) == mPassages.end())
            mPassages.push_back(aPassage);
    }

} // SeqdbSeq::add_passage

// ----------------------------------------------------------------------

void SeqdbSeq::update_gene(std::string aGene, Messages& aMessages, bool replace_ha)
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

void SeqdbSeq::add_reassortant(std::string aReassortant)
{
    // std::cerr << "DEBUG: add_reassortant " << aReassortant << '\n';
    if (!aReassortant.empty() && std::find(mReassortant.begin(), mReassortant.end(), aReassortant) == mReassortant.end())
        mReassortant.push_back(aReassortant);

} // SeqdbSeq::add_reassortant

// ----------------------------------------------------------------------

void SeqdbSeq::add_lab_id(std::string aLab, std::string aLabId)
{
    if (!aLab.empty()) {
        auto& lab_ids = mLabIds[aLab];
        if (!aLabId.empty() && std::find(lab_ids.begin(), lab_ids.end(), aLabId) == lab_ids.end()) {
            lab_ids.push_back(aLabId);
        }
    }

} // SeqdbSeq::add_lab_id

// ----------------------------------------------------------------------

AlignAminoAcidsData SeqdbSeq::align(bool aForce, Messages& aMessages, std::string name)
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
          align_data = translate_and_align(mNucleotides, aMessages, name);
          if (!align_data.amino_acids.empty())
              mAminoAcids = align_data.amino_acids;
          if (!align_data.shift.alignment_failed()) {
              update_gene(align_data.gene, aMessages, true);
              if (align_data.shift.aligned()) {
                  mAminoAcidsShift = align_data.shift;
                  mNucleotidesShift = - align_data.offset + align_data.shift * 3;
              }
          }
          // else {
          //     aMessages.warning() << "Nucs not translated/aligned\n"; //     AA:      " << mAminoAcids << '\n';
          // }
          break;
      case aling_amino_acids:
          mAminoAcidsShift.reset();
          align_data = align_amino_acids(mAminoAcids, aMessages);
          if (align_data.shift.aligned()) {
              mAminoAcidsShift = align_data.shift;
              update_gene(align_data.gene, aMessages, true);
          }
          else {
              aMessages.warning() << " AA not aligned\n"; // << mAminoAcids << '\n';
          }
          break;
    }

    return align_data;

} // SeqdbSeq::align

// ----------------------------------------------------------------------

const clades_t& SeqdbSeq::update_clades(std::string aVirusType, std::string aLineage, std::string aName)
{
    if (aligned()) {
        if (aVirusType == "B") {
            if (aLineage == "YAMAGATA") {
                mClades = clades_b_yamagata(mAminoAcids, mAminoAcidsShift, aName);
            }
            else if (aLineage == "VICTORIA") {
                mClades = clades_b_victoria(mAminoAcids, mAminoAcidsShift, aName);
            }
        }
        else if (aVirusType == "A(H1N1)") {
            mClades = clades_h1pdm(mAminoAcids, mAminoAcidsShift, aName);
        }
        else if (aVirusType == "A(H3N2)") {
            mClades = clades_h3n2(mAminoAcids, mAminoAcidsShift, aName);
        }
        // else {
        //     std::cerr << "Cannot update clades for virus type " << aVirusType << '\n';
        // }
    }
    return mClades;

} // SeqdbSeq::update_clades

// ----------------------------------------------------------------------

std::string SeqdbSeq::amino_acids(bool aAligned, size_t aLeftPartSize, size_t aResize) const
{
    std::string r = mAminoAcids;
    if (aAligned) {
        if (!aligned())
            throw SequenceNotAligned("SeqdbSeq::amino_acids()");
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

        if (aResize > 0)
            r.resize(aResize, 'X');
    }
    return r;

} // SeqdbSeq::amino_acids

// ----------------------------------------------------------------------

// aPos counts from 1!
char SeqdbSeq::amino_acid_at(size_t aPos, bool ignore_errors) const
{
    try {
    if (!aligned())
        throw SequenceNotAligned("SeqdbSeq::amino_acid_at()");
    const int offset = static_cast<int>(aPos) - 1 - mAminoAcidsShift;
    if (offset < 0 || offset >= static_cast<int>(mAminoAcids.size()))
        throw std::runtime_error("SeqdbSeq::amino_acid_at(): Invalid pos");
    return mAminoAcids[static_cast<size_t>(offset)];
    }
    catch (std::exception& err) {
        if (ignore_errors) {
            std::cerr << "WARNING: " << err.what() << '\n';
            return '?';
        }
        else {
            throw;
        }
    }

} // SeqdbSeq::amino_acid_at

// ----------------------------------------------------------------------

std::string SeqdbSeq::nucleotides(bool aAligned, size_t aLeftPartSize, size_t aResize) const
{
    std::string r = mNucleotides;
    if (aAligned) {
        if (!aligned())
            throw SequenceNotAligned("nucleotides()");
        r = shift(r, mNucleotidesShift + static_cast<int>(aLeftPartSize), '-');
        if (aResize > 0)
            r.resize(aResize, '-');
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

void SeqdbEntry::add_date(std::string aDate)
{
    if (!aDate.empty()) {
        auto insertion_pos = std::lower_bound(mDates.begin(), mDates.end(), aDate);
        if (insertion_pos == mDates.end() || aDate != *insertion_pos) {
            mDates.insert(insertion_pos, aDate);
        }
    }

} // SeqdbEntry::add_date

// ----------------------------------------------------------------------

void SeqdbEntry::update_lineage(std::string aLineage, Messages& aMessages)
{
      // std::cerr << "Lineage " << mName << " " << (aLineage.empty() ? std::string("?") : aLineage) << '\n';
    if (!aLineage.empty()) {
        if (mLineage.empty())
            mLineage = aLineage;
        else if (aLineage != mLineage)
            aMessages.warning() << mName << ": different lineages " << mLineage << " (stored) vs. " << aLineage << " (ignored)" << '\n';
    }

} // SeqdbEntry::update_lineage

// ----------------------------------------------------------------------

void SeqdbEntry::update_subtype_name(std::string aSubtype, Messages& aMessages)
{
    if (!aSubtype.empty() && aSubtype[0] != '*') { // do not update subtypes if it starts with *, it is a more general name (e.g. lacks N part)
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
                if (!(mVirusType == "A(H1N2)" && aSubtype == "A(H1N1)")) {
                    aMessages.warning() << mName << ": different subtypes " << mVirusType << " (stored) vs. " << aSubtype << " (ignored)" << '\n';
                }
            }
        }
          // fix subtype in the name too
        if (mName.find("A/") == 0)
            mName.replace(0, 1, mVirusType);
    }

} // SeqdbEntry::update_subtype_name

// ----------------------------------------------------------------------

const SeqdbSeq* SeqdbEntry::find_by_hi_name(std::string aHiName) const
{
    auto found = std::find_if(begin_seq(), end_seq(), [aHiName](const auto& seq) -> bool { return seq.hi_name_present(aHiName); });
    return found == end_seq() ? nullptr : &*found;

} // SeqdbEntry::find_by_hi_name

// ----------------------------------------------------------------------

// std::string SeqdbEntry::add_or_update_sequence(std::string aSequence, std::string aPassage, std::string aReassortant, std::string aLab, std::string aLabId, std::string aGene)
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

std::string Seqdb::add_sequence(std::string aName, std::string aVirusType, std::string aLineage, std::string aLab, std::string aDate, std::string aLabId, std::string aPassage,
                                std::string aReassortant, std::string aSequence, std::string aGene)
{
    Messages messages;
    virus_name::Name name_fields(aName);
    name_fields.fix_extra();
    try {
        ::virus_name::fix_location(name_fields);
    }
    catch (LocationNotFound&) {
        throw std::runtime_error("unrecognized location in " + aName);
    }
    const std::string name = name_fields.name_extra();
    SeqdbEntry entry(name, aVirusType, aLineage);

    SeqdbSeq new_seq(aSequence, aGene);
    const auto align_data = new_seq.align(false, messages, name);
    // std::cerr << "DEBUG: Seqdb::add_sequence: " << name << ' ' << aPassage << " nucs:" << new_seq.nucleotides_size() << " aa:" << new_seq.amino_acids_size() << '\n';
    if (!align_data.shift.aligned() && (new_seq.amino_acids_size() > MINIMUM_SEQUENCE_AA_LENGTH || new_seq.nucleotides_size() > (MINIMUM_SEQUENCE_AA_LENGTH * 3)))
        not_aligned_.emplace_back(aVirusType, name + ' ' + aPassage, new_seq.nucleotides_raw(), new_seq.amino_acids_raw());
    // std::cerr << "DEBUG: Seqdb::add_sequence: aligned: " << align_data.shift.aligned() << " nucs:" << new_seq.nucleotides_size() << '\n';
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

    if (aReassortant.empty() && !name_fields.reassortant.empty())
        inserted_seq->add_reassortant(name_fields.reassortant);
    else
        inserted_seq->add_reassortant(aReassortant);
    inserted_seq->add_passage(aPassage);
    inserted_seq->add_lab_id(aLab, aLabId);

    return messages;

} // Seqdb::add_sequence

// ----------------------------------------------------------------------

void Seqdb::report_not_aligned_after_adding() const
{
    if (!not_aligned_.empty()) {
        std::cerr << "WARNING: " << not_aligned_.size() << " sequences were not aligned\n";
        for (const auto& [virus_type, name, nucs, aa] : not_aligned_) {
            std::cerr << "    " << std::setw(9) << std::left << virus_type << ' ' << std::setw(60) << std::left << name << " aa:" << std::setw(3) << aa.size() << " nucs:" << nucs.size();
            if (!aa.empty())
                std::cerr << "   " << aa;
            std::cerr << '\n';
        }
    }

} // Seqdb::report_not_aligned_after_adding

// ----------------------------------------------------------------------

// SeqdbEntry* Seqdb::new_entry(std::string aName)
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
        size_t num_short_sequences = 0;
        std::for_each(mEntries.begin(), mEntries.end(), [&num_short_sequences](auto& entry) { if (entry.remove_short_sequences()) ++num_short_sequences; });
        if (num_short_sequences)
            std::cerr << "INFO: too short sequences removed: " << num_short_sequences << '\n';
    }

    {
        size_t num_not_translated_sequences = 0;
        std::for_each(mEntries.begin(), mEntries.end(), [&num_not_translated_sequences](auto& entry) { if (entry.remove_not_translated_sequences()) ++num_not_translated_sequences; });
        if (num_not_translated_sequences)
            std::cerr << "INFO: not translated sequences removed: " << num_not_translated_sequences << '\n';
    }

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

    auto report = [&os](std::string prefix, const auto& groups) {
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
    if (p_end - prefixes.begin()) {
        os << "WARNING: Prefixes of not aligned sequences of length " << prefix_size << ": " << p_end - prefixes.begin() << '\n';
        std::copy(prefixes.begin(), p_end, std::ostream_iterator<std::string>(os, "\n"));
    }
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
            antigen_index_list = hidb_antigens->find(name, hidb::fix_location::yes, hidb::find_fuzzy::no);
        }
        catch (LocationNotFound&) {
            if (std::count(name.begin(), name.end(), '/') == 4) {
                  // perhaps isolation split with / and there is no host
                  // in that case HI table usually has - instead of / (fixed manually by Eu on table parsing)
                const auto parts = acmacs::string::split(name, "/", acmacs::string::Split::KeepEmpty);
                name = string::join("/", {parts[0], parts[1], string::join("-", {parts[2], parts[3]}), parts[4]});
                antigen_index_list = hidb_antigens->find(name, hidb::fix_location::yes, hidb::find_fuzzy::no);
            }
            else {
                  // try without fixing location
                antigen_index_list = hidb_antigens->find(name, hidb::fix_location::no, hidb::find_fuzzy::no);
            }
        }
        std::transform(antigen_index_list.begin(), antigen_index_list.end(), std::back_inserter(found), [](const hidb::AntigenPIndex& antigen_index) -> hidb::AntigenP { return antigen_index.first; });
        std::sort(found.begin(), found.end());
        found.erase(std::unique(found.begin(), found.end()), found.end());
        // std::cerr << "DEBUG: find_in_hidb: " << found << '\n';
    }
    catch (hidb::get_error& /*err*/) {
        // std::cerr << "WARNING: no HiDb for " << entry.name() << ": " << err.what() << '\n';
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

SeqdbEntrySeq Seqdb::find_by_seq_id(std::string aSeqId, ignore_not_found ignore) const
{
    SeqdbEntrySeq result;
    const std::string seq_id = name_decode(aSeqId);
    auto passage_separator = seq_id.find("__");
    if (passage_separator != std::string::npos) { // seq_id
        if (const auto entry = find_by_name(std::string(seq_id, 0, passage_separator)); entry != nullptr) {
            const auto passage_distinct = acmacs::string::split(string::string_view(seq_id, passage_separator + 2), "__", acmacs::string::Split::KeepEmpty);
            auto index = passage_distinct.size() == 1 ? 0 : std::stoi(std::string(passage_distinct[1]));
            for (const auto& seq : entry->seqs()) {
                if (seq.passages().empty() ? passage_distinct[0].empty() : std::find(seq.passages().begin(), seq.passages().end(), passage_distinct[0]) != seq.passages().end()) {
                      // passage matched
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
            if (ignore == ignore_not_found::no)
                std::cerr << "ERROR: no entry for \"" << std::string(seq_id, 0, passage_separator) << "\" in seqdb [" << __FILE__ << ":" << __LINE__ << ']' << '\n';
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
        if (ignore == ignore_not_found::no)
            std::cerr << "ERROR: \"" << seq_id << "\" not in seqdb" << '\n';
    }
    return result;

} // Seqdb::find_by_seq_id

// ----------------------------------------------------------------------

void Seqdb::build_hi_name_index()
{
    mHiNameIndex.clear();
    for (auto entry_seq: *this) {
        for (const auto& hi_name: entry_seq.seq().hi_names()) {
              // const auto pos_inserted =
            mHiNameIndex.emplace(hi_name, entry_seq);
              // if (!pos_inserted.second)
              //     std::cerr << "warning:0: " << hi_name << " was already in mHiNameIndex [" << pos_inserted.first->second.entry().name() << "] [" << entry_seq.entry().name() << ']' << '\n';
        }
    }

} // Seqdb::build_hi_name_index

// ----------------------------------------------------------------------

size_t Seqdb::match(const acmacs::chart::Antigens& aAntigens, std::vector<SeqdbEntrySeq>& aPerAntigen, std::string aChartVirusType, seqdb::report aReport) const
{
    size_t matched = 0;
    aPerAntigen.clear();
    for (auto antigen: aAntigens) {
        bool destroy_entry = false;
        bool found = false;
        const SeqdbEntrySeq* entry = find_hi_name(antigen->full_name());
        if (!entry)
            entry = find_hi_name(antigen->full_name_for_seqdb_matching());
        if (!entry && antigen->passage().empty()) {
            if (const auto* s_entry = find_by_name(antigen->name()); s_entry) {
                for (const auto& seq : s_entry->seqs()) {
                    if (seq.reassortant_match(antigen->reassortant())) {
                        entry = new SeqdbEntrySeq(*s_entry, seq);
                        destroy_entry = true;
                    }
                }
            }
        }
        if (entry) {
            if (!aChartVirusType.empty() && aChartVirusType != entry->entry().virus_type()) {
                if (aReport == report::yes)
                    std::cerr << "WARNING: Seqdb::match: virus type mismatch: chart:" << aChartVirusType << " seq:" << entry->entry().virus_type() << " name: " << antigen->full_name() << '\n';
            }
            else if (antigen->lineage() != acmacs::chart::BLineage::Unknown && static_cast<std::string>(antigen->lineage()) != entry->entry().lineage()) {
                std::cerr << "WARNING: Seqdb::match: lineage mismatch: antigen:" << antigen->lineage() << " seq:" << entry->entry().lineage() << " name: " << antigen->full_name() << '\n';
            }
            else {
                aPerAntigen.push_back(*entry);
                found = true;
                ++matched;
            }
            if (destroy_entry)
                delete entry;
        }
        if (!found) {
            aPerAntigen.emplace_back();
            // if (aReport == report::yes)
            //     std::cerr << "WARNING: seqdb::match failed for \"" << antigen->full_name() << "\"" << '\n';
        }
    }
    if (aReport == report::yes)
        std::cerr << "INFO: " << matched << " antigens from chart have sequences in seqdb" << '\n';
    return matched;

} // Seqdb::match

// ----------------------------------------------------------------------

std::vector<SeqdbEntrySeq> Seqdb::match(const acmacs::chart::Antigens& aAntigens, std::string aChartVirusType, seqdb::report aReport) const
{
    std::vector<SeqdbEntrySeq> per_antigen;
    match(aAntigens, per_antigen, aChartVirusType, aReport);
    return per_antigen;

} // Seqdb::match

// ----------------------------------------------------------------------

clades_t Seqdb::clades_for_name(std::string name, Seqdb::clades_for_name_inclusive inclusive) const
{
    clades_t result;
    if (const auto* entry = find_by_name(name); entry) {
        bool clades_found = false;
        for (const auto& seq : entry->seqs()) {
            if (inclusive == clades_for_name_inclusive::yes || !clades_found) {
                for (const auto& clade : seq.clades())
                    result.push_back(clade);
            }
            else {
                result.erase(std::remove_if(std::begin(result), std::end(result), [&seq](const auto& clade) { return !seq.has_clade(clade); }), std::end(result));
            }
            clades_found |= !seq.clades().empty();
        }
    }
    return result;

} // Seqdb::clades_for_name

// ----------------------------------------------------------------------

void Seqdb::aa_at_positions_for_antigens(const acmacs::chart::Antigens& aAntigens, const std::vector<size_t>& aPositions, std::map<std::string, std::vector<size_t>>& aa_indices, seqdb::report aReport) const
{
    size_t matched = 0;
    for (auto ag = aAntigens.begin(); ag != aAntigens.end(); ++ag) {
        const SeqdbEntrySeq* entry = find_hi_name((*ag)->full_name());
        if (!entry)
            entry = find_hi_name((*ag)->full_name_for_seqdb_matching());
        if (entry) {
            std::string aa(aPositions.size(), 'X');
            std::transform(aPositions.begin(), aPositions.end(), aa.begin(), [&entry](size_t pos) { return entry->seq().amino_acid_at(pos, true); });
            aa_indices[aa].push_back(ag.index());
            ++matched;
        }
    }
    if (aReport == report::yes)
        std::cerr << "INFO: " << matched << " antigens from chart have sequences in seqdb" << '\n';

} // Seqdb::aa_at_positions_for_antigens

// ----------------------------------------------------------------------

void Seqdb::load(std::string filename)
{
    seqdb_import(filename, *this);
    mLoadedFromFilename = filename;

} // Seqdb::from_json_file

// ----------------------------------------------------------------------

void Seqdb::save(std::string filename, size_t indent) const
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
    std::cerr << "========== Deletions/insertions ==========\n";
    for (std::string virus_type: virus_types()) {
        if (!virus_type.empty()) {
            std::cout << "Detect insertions/deletions for " << virus_type << '\n';
            InsertionsDeletionsDetector detector(*this, virus_type);
            detector.detect();
        }
    }
    std::cerr << "========== Deletions/insertions done ==========\n";

} // Seqdb::detect_insertions_deletions

// ----------------------------------------------------------------------

void Seqdb::update_clades(seqdb::report aReport)
{
    std::cerr << "========== Clades ==========\n";
    std::map<std::string, size_t> clade_count;
    for (auto entry_seq: *this) {
        const auto& clades = entry_seq.seq().update_clades(entry_seq.entry().virus_type(), entry_seq.entry().lineage(), entry_seq.make_name());
        for (const auto& clade: clades)
            ++clade_count[clade];
    }
    if (aReport == report::yes)
        std::cerr << "INFO: clades: " << clade_count << '\n';
    std::cerr << "========== Clades done ==========\n";

} // Seqdb::update_clades

// ----------------------------------------------------------------------

void Seqdb::detect_b_lineage()
{
    BLineageDetector detector(*this);
    detector.detect();

} // Seqdb::detect_b_lineage

// ----------------------------------------------------------------------

void seqdb::add_clades(acmacs::chart::ChartModify& chart, ignore_errors ignore_err, report a_report)
{
    if (const auto& seqdb = get(ignore_err, report_time::no); seqdb) {
        auto antigens = chart.antigens_modify();
        const auto per_antigen = seqdb.match(*antigens, chart.info()->virus_type(acmacs::chart::Info::Compute::Yes), a_report);
        size_t sequenced = 0;
        for (auto ag_no : acmacs::range(per_antigen.size())) {
            if (const auto& entry_seq = per_antigen[ag_no]; entry_seq) {
                ++sequenced;
                auto& antigen = antigens->at(ag_no);
                if (const auto clades = entry_seq.seq().clades(); !clades.empty()) {
                    try {
                        for (const auto& clade : clades)
                            antigen.add_clade(clade);
                    }
                    catch (std::exception& err) {
                        std::cerr << "WARNING: cannot figure out clade for \"" << antigen.name() << "\": " << err.what() << '\n';
                    }
                    catch (...) {
                        std::cerr << "WARNING: cannot figure out clade for \"" << antigen.name() << "\": unknown exception\n";
                    }
                }
                else {
                    antigen.add_clade("SEQUENCED");
                }
            }
        }
        if (a_report == report::yes)
            std::cerr << "INFO: sequenced: " << sequenced << '\n';
    }

} // seqdb::add_clades

// ----------------------------------------------------------------------

struct StatPerPos
{
    ssize_t shannon_index = 0;  // https://en.wikipedia.org/wiki/Diversity_index
    std::map<char, size_t> aa_count;
};

std::string seqdb::sequences_of_chart_for_ace_view_1(acmacs::chart::Chart& chart)
{
    const auto matches = get().match(*chart.antigens(), chart.info()->virus_type());
    constexpr size_t max_num_pos = 1000;
    std::vector<StatPerPos> stat_per_pos(max_num_pos);
    auto json_antigens = to_json::object();
    for (auto [ag_no, entry_seq] : acmacs::enumerate(matches)) {
        if (entry_seq) {
            try {
                const auto sequence = entry_seq.seq().amino_acids(true);
                json_antigens = to_json::object_append(json_antigens, ag_no, sequence);
                for (auto [pos, aa] : acmacs::enumerate(sequence, 1))
                    ++stat_per_pos[pos].aa_count[aa];
            }
            catch (seqdb::SequenceNotAligned& err) {
                std::cerr << "WARNING: " << err.what() << ' ' << entry_seq.entry().name() << '\n';
            }
        }
    }
    for (auto& per_pos : stat_per_pos) {
        const auto sum = std::accumulate(per_pos.aa_count.begin(), per_pos.aa_count.end(), 0UL, [](auto accum, const auto& entry) { return accum + entry.second; });
        const auto shannon_index = -std::accumulate(per_pos.aa_count.begin(), per_pos.aa_count.end(), 0.0, [sum = double(sum)](auto accum, const auto& entry) {
            const double p = entry.second / sum;
            return accum + p * std::log(p);
        });
        per_pos.shannon_index = std::lround(shannon_index * 100);
    }
    auto json_per_pos = to_json::object();
    for (auto [pos, entry] : acmacs::enumerate(stat_per_pos)) {
        // if (entry.size() > 1) // && (entry.find('X') == entry.end() || entry.size() > 2))
        json_per_pos = to_json::object_append(json_per_pos, pos, to_json::raw(to_json::object("shannon", entry.shannon_index, "aa_count", to_json::raw(to_json::object(entry.aa_count)))));
    }
    return to_json::object("sequences", to_json::raw(to_json::object("antigens", to_json::raw(json_antigens), "per_pos", to_json::raw(json_per_pos))));

} // seqdb::sequences_of_chart_for_ace_view_1

// ----------------------------------------------------------------------

std::string seqdb::sequences_of_chart_as_fasta(acmacs::chart::Chart& chart)
{
    auto antigens = chart.antigens();
    const auto matches = get().match(*antigens, chart.info()->virus_type());
    std::string fasta;
    for (auto [ag_no, entry_seq] : acmacs::enumerate(matches)) {
        if (entry_seq) {
            try {
                const auto seq = entry_seq.seq().nucleotides(true);
                fasta += ">" + antigens->at(ag_no)->full_name() + '\n' + seq + '\n';
            }
            catch (seqdb::SequenceNotAligned& err) {
                std::cerr << "WARNING: " << err.what() << ' ' << entry_seq.entry().name() << '\n';
            }
        }
    }
    return fasta;

} // seqdb::sequences_of_chart_as_fasta

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
