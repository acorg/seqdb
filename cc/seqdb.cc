#include <typeinfo>
#include <tuple>
#include <set>

#include "acmacs-base/read-file.hh"
#include "acmacs-base/timeit.hh"
#include "acmacs-base/string-matcher.hh"
#include "seqdb.hh"
#include "clades.hh"
#include "seqdb-export.hh"

// ----------------------------------------------------------------------

bool SeqdbSeq::match_update_nucleotides(std::string aNucleotides)
{
    bool matches = false;
    if (mNucleotides == aNucleotides) {
        matches = true;
    }
    else if (mNucleotides.find(aNucleotides) != std::string::npos) { // sub
        matches = true;
    }
    else if (aNucleotides.find(mNucleotides) != std::string::npos) { // super
        matches = true;
        mNucleotides = aNucleotides;
        mNucleotidesShift.reset();
        mAminoAcidsShift.reset();
        mAminoAcids.clear();
    }
    return matches;

} // SeqdbSeq::match_update_nucleotides

// ----------------------------------------------------------------------

bool SeqdbSeq::match_update_amino_acids(std::string aAminoAcids)
{
    bool matches = false;
    if (mAminoAcids == aAminoAcids) {
        matches = true;
    }
    else if (mAminoAcids.find(aAminoAcids) != std::string::npos) { // sub
        matches = true;
    }
    else if (aAminoAcids.find(mAminoAcids) != std::string::npos) { // super
        matches = true;
        mNucleotides.clear();
        mNucleotidesShift.reset();
        mAminoAcidsShift.reset();
        mAminoAcids = aAminoAcids;
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
        if (aVirusType == "B" && aLineage == "YAMAGATA") {
            mClades = clades_b_yamagata(mAminoAcids, mAminoAcidsShift);
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
    auto insertion_pos = std::lower_bound(mDates.begin(), mDates.end(), aDate);
    if (insertion_pos == mDates.end() || aDate != *insertion_pos) {
        mDates.insert(insertion_pos, aDate);
    }

} // SeqdbEntry::add_date

// ----------------------------------------------------------------------

void SeqdbEntry::update_lineage(std::string aLineage, Messages& aMessages)
{
    if (!aLineage.empty()) {
        if (mLineage.empty())
            mLineage = aLineage;
        else if (aLineage != mLineage)
            aMessages.warning() << "Different lineages " << mLineage << " (stored) vs. " << aLineage << " (ignored)" << std::endl;
    }

} // SeqdbEntry::update_lineage

// ----------------------------------------------------------------------

void SeqdbEntry::update_subtype(std::string aSubtype, Messages& aMessages)
{
    if (!aSubtype.empty()) {
        if (mVirusType.empty())
            mVirusType = aSubtype;
        else if (aSubtype != mVirusType)
            aMessages.warning() << "Different subtypes " << mVirusType << " (stored) vs. " << aSubtype << " (ignored)" << std::endl;
    }

} // SeqdbEntry::update_subtype

// ----------------------------------------------------------------------

std::string SeqdbEntry::add_or_update_sequence(std::string aSequence, std::string aPassage, std::string aReassortant, std::string aLab, std::string aLabId, std::string aGene)
{
    Messages messages;
    const bool nucs = is_nucleotides(aSequence);
    decltype(mSeq.begin()) found;
    if (nucs)
        found = std::find_if(mSeq.begin(), mSeq.end(), [&aSequence](SeqdbSeq& seq) { return seq.match_update_nucleotides(aSequence); });
    else
        found = std::find_if(mSeq.begin(), mSeq.end(), [&aSequence](SeqdbSeq& seq) { return seq.match_update_amino_acids(aSequence); });
    if (found != mSeq.end()) {  // update
        found->update_gene(aGene, messages);
    }
    else { // add
        if (nucs)
            mSeq.push_back(SeqdbSeq(aSequence, aGene));
        else
            mSeq.push_back(SeqdbSeq(std::string(), aSequence, aGene));
        found = mSeq.end() - 1;
    }
    if (found != mSeq.end()) {
        found->add_passage(aPassage);
        found->add_reassortant(aReassortant);
        found->add_lab_id(aLab, aLabId);
        auto align_data = found->align(false, messages);
        update_subtype(align_data.subtype, messages);
        update_lineage(align_data.lineage, messages);
    }
    return messages;

} // SeqdbEntry::add_or_update_sequence

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

SeqdbEntry* Seqdb::new_entry(std::string aName)
{
    auto const first = find_insertion_place(aName);
    if (first != mEntries.end() && aName == first->name())
        throw std::runtime_error(std::string("Entry for \"") + aName + "\" already exists");
    auto inserted = mEntries.insert(first, SeqdbEntry(aName));
    return &*inserted;

} // Seqdb::new_entry

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
      //std::remove_if(mEntries.begin(), mEntries.end(), [](auto entry) { return entry.empty(); });
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
    // os << "Virus types: " << json::dump(virus_types) << std::endl;

    std::map<std::string, size_t> lineages;
    for_each(mEntries.begin(), mEntries.end(), [&lineages](auto& e) { ++lineages[e.lineage()]; });
    // os << "Lineages: " << json::dump(lineages) << std::endl;

    std::map<std::string, size_t> aligned;
    for_each(mEntries.begin(), mEntries.end(), [&aligned](auto& e) { if (any_of(e.mSeq.begin(), e.mSeq.end(), std::mem_fn(&SeqdbSeq::aligned))) ++aligned[e.virus_type()]; });
    // os << "Aligned: " << json::dump(aligned) << std::endl;

    std::map<std::string, size_t> matched;
    for_each(mEntries.begin(), mEntries.end(), [&matched](auto& e) { if (any_of(e.mSeq.begin(), e.mSeq.end(), std::mem_fn(&SeqdbSeq::matched))) ++matched[e.virus_type()]; });
    // os << "Matched: " << json::dump(matched) << std::endl;

    std::map<std::string, size_t> have_dates;
    for_each(mEntries.begin(), mEntries.end(), [&have_dates](auto& e) { if (!e.mDates.empty()) ++have_dates[e.virus_type()]; });
    // os << "Have dates: " << json::dump(have_dates) << std::endl;

    std::map<std::string, size_t> have_clades;
    for_each(mEntries.begin(), mEntries.end(), [&have_clades](auto& e) { if (any_of(e.mSeq.begin(), e.mSeq.end(), [](auto& seq) { return !seq.mClades.empty(); })) ++have_clades[e.virus_type()]; });
    // os << "Have clades: " << json::dump(have_clades) << std::endl;

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

void Seqdb::find_in_hidb(std::vector<const hidb::AntigenData*>& found, const SeqdbEntry& entry, HiDbPtrs& hidb_ptrs, std::string aHiDbDir) const
{
    try {
        const hidb::HiDb& hidb = get_hidb(entry.virus_type(), hidb_ptrs, aHiDbDir);
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
    }
    catch (NoHiDb&) {
    }

} // Seqdb::find_in_hidb

// ----------------------------------------------------------------------

void Seqdb::match_hidb(std::string aHiDbDir, bool aVerbose)
{
    HiDbPtrs hidb_ptrs;
    std::ostream& report_stream = std::cerr;

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
                out << " [" << std::get<0>(sf).first << " " << std::get<0>(sf).second << " S:" << std::get<1>(sf) << " F:" << std::get<2>(sf) << ']';
            }
            out << std::endl;
        }
    };

    std::vector<const SeqdbEntry*> not_matched;
    for (auto& entry: mEntries) {
        std::vector<const hidb::AntigenData*> found;
        find_in_hidb(found, entry, hidb_ptrs, aHiDbDir);

        if (aVerbose)
            report_stream << std::endl << entry << std::endl;
        if (!found.empty()) {
            typedef std::pair<string_match::score_t, size_t> score_size_t;
            typedef std::tuple<score_size_t, size_t, size_t> score_seq_found_t; // [score, len], seq_no, found_no

            std::vector<std::vector<score_seq_found_t>> matching; // for each seq list of matching [[score, min passage len], found_no] - sorted by score desc
            size_t seq_no = 0;
            for (auto& seq: entry.seqs()) {
                std::vector<score_seq_found_t> matching_for_seq;
                size_t found_no = 0;
                for (const auto& f: found) {
                    if (seq.reassortant_match(f->data().reassortant())) {
                        std::vector<score_size_t> scores; // score and min passage length (to avoid too incomplete matches)
                        const auto& f_passage = f->data().passage();
                        if (!seq.passages().empty())
                            std::transform(seq.passages().begin(), seq.passages().end(), std::back_inserter(scores),
                                           [&f_passage](const auto& passage) { return std::make_pair(string_match::match(passage, f_passage), std::min(passage.size(), f_passage.size())); });
                        else
                            scores.emplace_back(string_match::match(std::string{}, f_passage), 0);
                        matching_for_seq.emplace_back(*std::max_element(scores.begin(), scores.end(), [](const auto& a, const auto& b) { return a.first < b.first; }), seq_no, found_no);
                          // report_stream << "  @" << seq.passages() << " @ " << f_passage << " " << score_size->first << " " << score_size->second << std::endl;
                    }
                    ++found_no;
                }
                std::sort(matching_for_seq.begin(), matching_for_seq.end(), [](const auto& a, const auto& b) { return std::get<0>(a).first > std::get<0>(b).first; });
                matching.push_back(std::move(matching_for_seq));
                ++seq_no;
            }
            std::sort(matching.begin(), matching.end(), [](const auto& a, const auto& b) { return std::get<0>(a[0]).first > std::get<0>(b[0]).first; });

            if (matching.size() == 1) {
                for (const auto& ms: matching[0]) {
                    if (std::get<0>(ms).first == std::get<0>(matching[0][0]).first) {
                        const auto name = found[std::get<2>(ms)]->data().full_name();
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
                        if (std::get<0>(sf).first == std::get<0>(m[0]).first && found_assigned.count(std::get<2>(sf)) == 0) {
                            const auto name = found[std::get<2>(sf)]->data().full_name();
                            entry.seqs()[std::get<1>(sf)].add_hi_name(name);
                            if (aVerbose)
                                report_stream << "    +" << std::get<1>(sf) << " " << name << std::endl;
                            found_assigned.insert(std::get<2>(sf));
                        }
                    }
                }
            }
        }
        else {
            if (aVerbose)
                report_stream << "  ??" << std::endl;
            not_matched.push_back(&entry);
        }
    }

    std::cout << "Matched " << (mEntries.size() - not_matched.size()) << " of " << mEntries.size() << "  " << ((mEntries.size() - not_matched.size()) * 100.0 / mEntries.size()) << '%' << std::endl;

    std::cerr << "Not matched " << not_matched.size() << std::endl;
    if (!not_matched.empty()) {
        std::cerr << "  ";
        std::transform(not_matched.begin(), not_matched.end(), std::ostream_iterator<SeqdbEntry>(std::cerr, "\n  "), [](const auto* p) -> const SeqdbEntry& { return *p; });
    }

} // Seqdb::match_hidb

// ----------------------------------------------------------------------

const hidb::HiDb& Seqdb::get_hidb(std::string aVirusType, HiDbPtrs& aPtrs, std::string aHiDbDir) const
{
    auto h = aPtrs.find(aVirusType);
    if (h == aPtrs.end()) {
        std::string filename;
        if (aVirusType == "A(H1N1)")
            filename = aHiDbDir + "/hidb4.h1.json.xz";
        else if (aVirusType == "A(H3N2)")
            filename = aHiDbDir + "/hidb4.h3.json.xz";
        else if (aVirusType == "B")
            filename = aHiDbDir + "/hidb4.b.json.xz";
        else
            throw NoHiDb{};
          //throw std::runtime_error("No HiDb for " + aVirusType);

        std::unique_ptr<hidb::HiDb> hidb{new hidb::HiDb{}};
          // std::cout << "opening " << filename << std::endl;
        hidb->importFrom(filename);
        hidb->importLocDb(aHiDbDir + "/locationdb.json.xz");
        h = aPtrs.emplace(aVirusType, std::move(hidb)).first;
    }
    return *h->second;

} // Seqdb::get_hidb

// ----------------------------------------------------------------------

SeqdbEntrySeq Seqdb::find_by_seq_id(std::string aSeqId) const
{
    SeqdbEntrySeq result;
    const auto passage_separator = aSeqId.find("__");
    if (passage_separator != std::string::npos) { // seq_id
        const auto entry = find_by_name(std::string(aSeqId, 0, passage_separator));
        if (entry != nullptr) {
            const auto passage_distinct = string::split(std::string(aSeqId, passage_separator + 2), "__", string::Split::KeepEmpty);
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
    }
    else {
        std::smatch year_space;
        const auto year_space_present = std::regex_search(aSeqId, year_space, sReYearSpace);
        const std::string look_for = year_space_present ? std::string(aSeqId, 0, static_cast<std::string::size_type>(year_space.position(0) + year_space.length(0)) - 1) : aSeqId;
        const auto entry = find_by_name(look_for);
        if (entry != nullptr) {
            auto found = std::find_if(entry->begin_seq(), entry->end_seq(), [&aSeqId](const auto& seq) -> bool { return seq.hi_name_present(aSeqId); });
            if (found == entry->end_seq()) { // not found by hi_name, look by passage (or empty passage)
                const std::string passage = year_space_present ? std::string(aSeqId, static_cast<std::string::size_type>(year_space.position(0) + year_space.length(0))) : std::string();
                found = std::find_if(entry->begin_seq(), entry->end_seq(), [&passage](const auto& seq) -> bool { return seq.passage_present(passage); });
            }
            if (found != entry->end_seq()) {
                result.assign(*entry, *found);
            }
        }
    }

    if (!result) {
        std::cerr << "Error: \"" << aSeqId << "\" not in seqdb" << std::endl;
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
