#include <iostream>

#include "acmacs-base/string-matcher.hh"
#include "acmacs-base/passage.hh"

#include "seqdb/seqdb.hh"

using namespace seqdb;

// ----------------------------------------------------------------------

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

using Found = std::vector<const hidb::AntigenData*>;
using Matching = std::vector<std::vector<score_seq_found_t>>;

// ----------------------------------------------------------------------

static void report_found(std::ostream& out, const Found& found)
{
    size_t found_no = 0;
    for (const auto& e: found) {
        out << "  >> " << found_no << ' ' << e->data().full_name() << std::endl;
        ++found_no;
    }
}

static void report_matching(std::ostream& out, const Matching& matching)
{
    for (const auto& m: matching) {
        out << "    matching: ";
        for (const auto& sf: m) {
            out << " [" << sf.score << " " << sf.len << " S:" << sf.seq_no << " F:" << sf.found_no << ']';
        }
        out << std::endl;
    }
}

// ----------------------------------------------------------------------

static void make_matching(SeqdbEntry& entry, const Found& found, Matching& matching)
{
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
}

// ----------------------------------------------------------------------

// greedy matching: add all hi-names having matching reassortant and passage type (egg/cell) regardless of score
// if antigen is in multiple matching entries, use the one with the highest score
static void match_greedy(SeqdbEntry& entry, const Found& found, const Matching& matching, bool aVerbose, std::ostream& report_stream)
{
    if (aVerbose) {
        report_found(report_stream, found);
        report_matching(report_stream, matching);
    }

    std::map<size_t, std::pair<decltype(matching.begin() - matching.begin()), string_match::score_t>> antigen_to_matching; // antigen index in found to (matching index and score)
    for (auto mp = matching.begin(); mp != matching.end(); ++mp) {
        for (const auto& sf: *mp) {
            const auto ampi = antigen_to_matching.emplace(sf.found_no, std::make_pair(mp - matching.begin(), sf.score));
            if (!ampi.second && ampi.first->second.second < sf.score) {        // already present, replace if sf has hi higher score
                ampi.first->second.first = mp - matching.begin();
                ampi.first->second.second = sf.score;
            }
        }
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

// ----------------------------------------------------------------------

static void match_normal(SeqdbEntry& entry, const Found& found, const Matching& matching, bool aVerbose, std::ostream& report_stream)
{
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

// ----------------------------------------------------------------------

void Seqdb::match_hidb(bool aVerbose, bool aGreedy)
{
    std::ostream& report_stream = std::cerr;

    std::vector<const SeqdbEntry*> not_matched;
    for (auto& entry: mEntries) {
        Found found;
        find_in_hidb_update_country_lineage_date(found, entry);

        if (aVerbose)
            report_stream << std::endl << entry << std::endl;
        if (!found.empty()) {
            Matching matching; // for each seq list of matching [[score, min passage len], found_no] - sorted by score desc
            make_matching(entry, found, matching);
            if (aGreedy)
                match_greedy(entry, found, matching, aVerbose, report_stream);
            else
                match_normal(entry, found, matching, aVerbose, report_stream);
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
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
