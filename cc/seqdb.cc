#include <typeinfo>
#include <tuple>
#include <set>

#include "acmacs-base/read-file.hh"
#include "acmacs-base/timeit.hh"
#include "acmacs-base/string-matcher.hh"
#include "acmacs-base/virus-name.hh"
#include "acmacs-base/passage.hh"
#include "locationdb/locdb.hh"
#include "seqdb/seqdb.hh"
#include "clades.hh"
#include "seqdb-export.hh"
#include "seqdb-import.hh"

using namespace seqdb;

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
    if (!aPassage.empty() && std::find(mPassages.begin(), mPassages.end(), aPassage) == mPassages.end())
        mPassages.push_back(aPassage);

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
                aMessages.warning() << "[SAMESEQ] different genes " << mGene << " vs. " << aGene << std::endl;
        }
    }

} // SeqdbSeq::update_gene

// ----------------------------------------------------------------------

void SeqdbSeq::add_reassortant(std::string aReassortant)
{
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
              aMessages.warning() << "Nucs not translated/aligned" /* << mNucleotides */ << std::endl;
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
              aMessages.warning() << "AA not aligned" /* << mAminoAcids */ << std::endl;
          }
          break;
    }

    return align_data;

} // SeqdbSeq::align

// ----------------------------------------------------------------------

void SeqdbSeq::update_clades(std::string aVirusType, std::string aLineage)
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
        //     std::cerr << "Cannot update clades for virus type " << aVirusType << std::endl;
        // }
    }

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
      // std::cerr << "Lineage " << mName << " " << (aLineage.empty() ? std::string("?") : aLineage) << std::endl;
    if (!aLineage.empty()) {
        if (mLineage.empty())
            mLineage = aLineage;
        else if (aLineage != mLineage)
            aMessages.warning() << "Different lineages " << mLineage << " (stored) vs. " << aLineage << " (ignored)" << std::endl;
    }

} // SeqdbEntry::update_lineage

// ----------------------------------------------------------------------

void SeqdbEntry::update_subtype_name(std::string aSubtype, Messages& aMessages)
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
                aMessages.warning() << "Different subtypes " << mVirusType << " (stored) vs. " << aSubtype << " (ignored)" << std::endl;
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

std::string Seqdb::add_sequence(std::string aName, std::string aVirusType, std::string aLineage, std::string aLab, std::string aDate, std::string aLabId, std::string aPassage, std::string aReassortant, std::string aSequence, std::string aGene)
{
    Messages messages;
    std::string name = virus_name::normalize(aName);
    SeqdbEntry entry(name, aVirusType, aLineage);

    SeqdbSeq new_seq(aSequence, aGene);
    auto align_data = new_seq.align(false, messages);
    entry.update_subtype_name(align_data.subtype, messages); // updates entry.mName!
    // std::cerr << "add " << align_data.subtype << ' ' << entry.name() << std::endl;

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
        std::for_each(mEntries.begin(), mEntries.end(), std::mem_fn(&SeqdbEntry::remove_short_sequences));
    }

    std::for_each(mEntries.begin(), mEntries.end(), std::mem_fn(&SeqdbEntry::remove_not_translated_sequences));

      // Note empty passages must not be removed, this is just for testing purposes
      // std::for_each(mEntries.begin(), mEntries.end(), std::mem_fn(&SeqdbEntry::remove_empty_passages));

      // remove empty entries
    auto const num_entries_before = mEntries.size();
    mEntries.erase(std::remove_if(mEntries.begin(), mEntries.end(), std::mem_fn(&SeqdbEntry::empty)), mEntries.end());
    if (mEntries.size() != num_entries_before)
        messages.warning() << (num_entries_before - mEntries.size()) << " entries removed during cleanup" << std::endl;
    return messages;

} // Seqdb::cleanup

// ----------------------------------------------------------------------

std::string Seqdb::report() const
{
    std::ostringstream os;
    os << "Entries: " << mEntries.size() << std::endl;

    std::map<std::string, size_t> virus_types;
    for_each(mEntries.begin(), mEntries.end(), [&virus_types](auto& e) { ++virus_types[e.virus_type()]; });
    os << "Virus types: " << virus_types << std::endl;

    std::map<std::string, size_t> lineages;
    for_each(mEntries.begin(), mEntries.end(), [&lineages](auto& e) { ++lineages[e.lineage()]; });
    os << "Lineages: " << lineages << std::endl;

    std::map<std::string, size_t> aligned;
    for_each(mEntries.begin(), mEntries.end(), [&aligned](auto& e) { if (any_of(e.mSeq.begin(), e.mSeq.end(), std::mem_fn(&SeqdbSeq::aligned))) ++aligned[e.virus_type()]; });
    os << "Aligned: " << aligned << std::endl;

    std::map<std::string, size_t> matched;
    for_each(mEntries.begin(), mEntries.end(), [&matched](auto& e) { if (any_of(e.mSeq.begin(), e.mSeq.end(), std::mem_fn(&SeqdbSeq::matched))) ++matched[e.virus_type()]; });
    os << "Matched: " << matched << std::endl;

    std::map<std::string, size_t> have_dates;
    for_each(mEntries.begin(), mEntries.end(), [&have_dates](auto& e) { if (!e.mDates.empty()) ++have_dates[e.virus_type()]; });
    os << "Have dates: " << have_dates << std::endl;

    std::map<std::string, size_t> have_clades;
    for_each(mEntries.begin(), mEntries.end(), [&have_clades](auto& e) { if (any_of(e.mSeq.begin(), e.mSeq.end(), [](auto& seq) { return !seq.mClades.empty(); })) ++have_clades[e.virus_type()]; });
    os << "Have clades: " << have_clades << std::endl;

    return os.str();

} // Seqdb::report

// ----------------------------------------------------------------------

std::string Seqdb::report_identical() const
{
    std::ostringstream os;

    auto report = [&os](std::string prefix, const auto& groups) {
        if (!groups.empty()) {
            os << prefix << std::endl;
            for (auto const& group: groups) {
                for (auto const& entry: group) {
                    os << entry.make_name() << ' ';
                }
                os << std::endl;
            }
        }
    };

    try {
        report("Identical nucleotides:", find_identical_sequences([](const SeqdbEntrySeq& e) -> std::string { try { return e.seq().nucleotides(true); } catch (SequenceNotAligned&) { return std::string(); } }));
        os << std::endl;
    }
    catch (std::exception& err) {
        os << "Cannot report identical nucleotides: " << typeid(err).name() << ": " << err.what() << std::endl;
    }

    try {
        report("Identical amino-acids:", find_identical_sequences([](const SeqdbEntrySeq& e) -> std::string { try { return e.seq().amino_acids(true); } catch (SequenceNotAligned&) { return std::string(); } }));
        os << std::endl;
    }
    catch (std::exception& err) {
        os << "Cannot report identical amino-acids: " << typeid(err).name() << ": " << err.what() << std::endl;
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
    os << "Prefixes of not aligned sequences of length " << prefix_size << ": " << p_end - prefixes.begin() << std::endl;
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

void Seqdb::find_in_hidb_update_country_lineage_date(std::vector<const hidb::AntigenData*>& found, SeqdbEntry& entry) const
{
    try {
        const hidb::HiDb& hidb = mHiDbSet.get(entry.virus_type());
        const auto cdcids = entry.cdcids();
        if (!cdcids.empty()) {
            for (const auto& cdcid: cdcids) {
                const auto f_cdcid = hidb.find_antigens_by_cdcid(cdcid);
                std::copy(f_cdcid.begin(), f_cdcid.end(), std::back_inserter(found));
            }
        }
        const auto f_name = hidb.find_antigens_by_name(entry.name());
        std::copy(f_name.begin(), f_name.end(), std::back_inserter(found));

        std::sort(found.begin(), found.end());
        found.erase(std::unique(found.begin(), found.end()), found.end());

          // update country and continent
        if (entry.country().empty()) {
            try {
                const std::string country = hidb.locdb().find(virus_name::location(entry.name())).country();
                entry.country(country);
                entry.continent(hidb.locdb().continent_of_country(country));
            }
            catch (LocationNotFound&) {
            }
            catch (virus_name::Unrecognized&) {
            }
        }

          // update lineage
        if (entry.virus_type() == "B" && !found.empty()) {
            Messages messages;
            entry.update_lineage(found.front()->data().lineage(), messages);
            if (messages)
                std::cerr << messages << std::endl;
        }

          // update date
        if (!found.empty()) {
            for (const auto& e: found) {
                entry.add_date(e->date());
            }
        }
    }
    catch (hidb::NoHiDb&) {
    }

} // Seqdb::find_in_hidb_update_country_lineage_date

// ----------------------------------------------------------------------

void seqdb::Seqdb::match_hidb(bool aVerbose, bool aGreedy)
{
    std::ostream& report_stream = std::cerr;

    struct score_size_t
    {
        string_match::score_t score;
        size_t len;

          //inline score_size_t() = default;
        inline score_size_t(string_match::score_t s, size_t l) : score(s), len(l) {}
        inline bool operator<(const score_size_t& a) const { return score < a.score; }
    };

    struct score_seq_found_t : public score_size_t
    {
        size_t seq_no;
        size_t found_no;

          //inline score_seq_found_t() = default;
        inline score_seq_found_t(const score_size_t& ss, size_t sn, size_t fn) : score_size_t(ss.score, ss.len), seq_no(sn), found_no(fn) {}
        inline bool operator<(const score_seq_found_t& a) const { return score > a.score; }
    };

    auto report_found = [](std::ostream& out, const auto& found) {
        size_t found_no = 0;
        for (const auto& e: found) {
            out << "  >> " << found_no << ' ' << e->data().full_name() << std::endl;
            ++found_no;
        }
    };

    auto report_matching = [](std::ostream& out, const auto& matching) {
        for (const auto& m: matching) {
            out << "    ";
            for (const auto& sf: m) {
                out << " [" << sf.score << " " << sf.len << " S:" << sf.seq_no << " F:" << sf.found_no << ']';
            }
            out << std::endl;
        }
    };

    std::vector<const SeqdbEntry*> not_matched;
    for (auto& entry: mEntries) {
        std::vector<const hidb::AntigenData*> found;
        find_in_hidb_update_country_lineage_date(found, entry);

        if (aVerbose)
            report_stream << std::endl << entry << std::endl;
        if (!found.empty()) {

            std::vector<std::vector<score_seq_found_t>> matching; // for each seq list of matching [[score, min passage len], found_no] - sorted by score desc
            size_t seq_no = 0;
            for (auto& seq: entry.seqs()) {
                std::vector<score_seq_found_t> matching_for_seq;
                size_t found_no = 0;
                const passage::CellOrEgg seq_cell_or_egg = passage::cell_or_egg(seq.passages());
                for (const auto& f: found) {
                    const auto& f_passage = f->data().passage();
                    if (seq.reassortant_match(f->data().reassortant()) && passage::match_cell_egg(passage::cell_or_egg(f_passage), seq_cell_or_egg)) {
                        std::vector<score_size_t> scores; // score and min passage length (to avoid too incomplete matches)
                        if (!seq.passages().empty())
                            std::transform(seq.passages().begin(), seq.passages().end(), std::back_inserter(scores),
                                           [&f_passage](const auto& passage) -> score_size_t { return {string_match::match(passage, f_passage), std::min(passage.size(), f_passage.size())}; });
                        else
                            scores.emplace_back(string_match::match(std::string{}, f_passage), 0);
                        matching_for_seq.emplace_back(*std::max_element(scores.begin(), scores.end() /*$, [](const auto& a, const auto& b) { return a.first < b.first; }*/), seq_no, found_no);
                          // report_stream << "  @" << seq.passages() << " @ " << f_passage << " " << score_size->first << " " << score_size->second << std::endl;
                    }
                    ++found_no;
                }
                std::sort(matching_for_seq.begin(), matching_for_seq.end());
                matching.push_back(std::move(matching_for_seq));
                ++seq_no;
            }
            std::sort(matching.begin(), matching.end(), [](const auto& a, const auto& b) -> bool { return a.empty() ? false : (b.empty() || a[0] < b[0]); });

            if (aGreedy) {
                  // greedy matching: add all hi-names having matching reassortant and passage type (egg/cell) regardless of score
                if (aVerbose) {
                    report_found(report_stream, found);
                    report_matching(report_stream, matching);
                }
                std::set<size_t> found_assigned;
                for (const auto& m: matching) {
                    for (const auto& sf: m) {
                        const auto& antigen = found[sf.found_no]->data();
                        if (found_assigned.count(sf.found_no) == 0) {
                            const auto name = antigen.full_name();
                            entry.seqs()[sf.seq_no].add_hi_name(name);
                            if (aVerbose)
                                report_stream << "    +" << sf.seq_no << " " << name << std::endl;
                            found_assigned.insert(sf.found_no);
                        }
                    }
                }
            }
            else {
                if (matching.size() == 1) {
                    for (const auto& ms: matching[0]) {
                        if (ms.score == matching[0][0].score) {
                            const auto name = found[ms.found_no]->data().full_name();
                            entry.seqs()[0].add_hi_name(name);
                            if (aVerbose)
                                report_stream << "    + " << name << std::endl;
                        }
                    }
                }
                else {
                    if (aVerbose) {
                        report_found(report_stream, found);
                        report_matching(report_stream, matching);
                    }
                    std::set<size_t> found_assigned;
                    for (const auto& m: matching) {
                        for (const auto& sf: m) {
                            if (sf.score == m[0].score && found_assigned.count(sf.found_no) == 0) {
                                const auto name = found[sf.found_no]->data().full_name();
                                entry.seqs()[sf.seq_no].add_hi_name(name);
                                if (aVerbose)
                                    report_stream << "    +" << sf.seq_no << " " << name << std::endl;
                                found_assigned.insert(sf.found_no);
                            }
                        }
                    }
                }
            }
        }
        else {
            if (aVerbose)
                report_stream << "  ?? " << entry.name() << std::endl;
            not_matched.push_back(&entry);
        }
        if (aVerbose)
            report_stream << std::endl;
    }

    std::cout << "Matched " << (mEntries.size() - not_matched.size()) << " of " << mEntries.size() << "  " << ((mEntries.size() - not_matched.size()) * 100.0 / mEntries.size()) << '%' << std::endl;

    if (aVerbose && !not_matched.empty()) {
        report_stream << "Not matched " << not_matched.size() << std::endl;
        for (const auto& nm: not_matched) {
            report_stream << "  " << *nm << std::endl;
        }
    }

} // Seqdb::match_hidb

// ----------------------------------------------------------------------

SeqdbEntrySeq Seqdb::find_by_seq_id(std::string aSeqId) const
{
    SeqdbEntrySeq result;
    const std::string seq_id = name_decode(aSeqId);
    auto passage_separator = seq_id.find("__");
    if (passage_separator != std::string::npos) { // seq_id
        const auto entry = find_by_name(std::string(seq_id, 0, passage_separator));
        if (entry != nullptr) {
            const auto passage_distinct = string::split(std::string(seq_id, passage_separator + 2), "__", string::Split::KeepEmpty);
            auto index = passage_distinct.size() == 1 ? 0 : std::stoi(passage_distinct[1]);
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
            std::cerr << "Error: no entry for \"" << std::string(seq_id, 0, passage_separator) << "\" in seqdb [" << __FILE__ << ":" << __LINE__ << ']' << std::endl;
        }
    }
    else {
        std::smatch year_space;
        const auto year_space_present = std::regex_search(seq_id, year_space, sReYearSpace);
        const std::string look_for = year_space_present ? std::string(seq_id, 0, static_cast<std::string::size_type>(year_space.position(0) + year_space.length(0)) - 1) : seq_id;
        const auto entry = find_by_name(look_for);
        if (entry != nullptr) {
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
        std::cerr << "Error: \"" << seq_id << "\" not in seqdb" << std::endl;
    }
    return result;

} // Seqdb::find_by_seq_id

// ----------------------------------------------------------------------

void Seqdb::load(std::string filename)
{
    seqdb_import(filename, *this);

} // Seqdb::from_json_file

// ----------------------------------------------------------------------

void Seqdb::save(std::string filename, size_t indent) const
{
    seqdb_export(filename, *this, indent);

} // Seqdb::save

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
