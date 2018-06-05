#include <iostream>

#include "acmacs-base/string-matcher.hh"
#include "acmacs-base/passage.hh"
#include "locationdb/locdb.hh"
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

using Found = std::vector<hidb::AntigenP>;
using Matching = std::vector<std::vector<score_seq_found_t>>;

// ----------------------------------------------------------------------

static void report_found(std::ostream& out, const Found& found)
{
    size_t found_no = 0;
    for (const auto& e: found) {
        out << "  >> " << found_no << ' ' << e->full_name() << '\n';
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
        out << '\n';
    }
}

// ----------------------------------------------------------------------

static void make_matching(SeqdbEntry& entry, const Found& found, Matching& matching)
{
    size_t seq_no = 0;
    for (auto& seq: entry.seqs()) {
        std::vector<score_seq_found_t> matching_for_seq;
        size_t found_no = 0;
        const auto seq_cell_or_egg = acmacs::passage::cell_or_egg(seq.passages());
        for (auto f: found) {
            const auto f_passage = f->passage();
              // std::cerr << "match_cell_egg: " << acmacs::passage::match_cell_egg(acmacs::passage::cell_or_egg(f_passage), seq_cell_or_egg) << " -- " << f_passage << ':' << static_cast<int>(passage::cell_or_egg(f_passage)) << " " << seq.passages() << ':' << static_cast<int>(seq_cell_or_egg) << '\n';
            if (seq.reassortant_match(f->reassortant()) && acmacs::passage::match_cell_egg(acmacs::passage::cell_or_egg(f_passage), seq_cell_or_egg)) {
                std::vector<score_size_t> scores; // score and min passage length (to avoid too incomplete matches)
                if (!seq.passages().empty())
                    std::transform(seq.passages().begin(), seq.passages().end(), std::back_inserter(scores),
                                   [&f_passage](const auto& passage) -> score_size_t { return {string_match::match(passage, f_passage), std::min(passage.size(), f_passage.size())}; });
                else
                    scores.emplace_back(string_match::match(std::string{}, f_passage), 0);
                matching_for_seq.emplace_back(*std::max_element(scores.begin(), scores.end() /*$, [](const auto& a, const auto& b) { return a.first < b.first; }*/), seq_no, found_no);
                  // report_stream << "  @" << seq.passages() << " @ " << f_passage << " " << score_size->first << " " << score_size->second << '\n';
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
// returns if at least one seq matched
static bool match_greedy(SeqdbEntry& entry, const Found& found, const Matching& matching, seqdb::report aReport, std::ostream& report_stream)
{
    if (aReport == seqdb::report::yes) {
        report_found(report_stream, found);
        report_matching(report_stream, matching);
    }

    std::map<size_t, score_seq_found_t> antigen_to_matching; // antigen index in found to (matching index and score)
    for (auto mp = matching.begin(); mp != matching.end(); ++mp) {
        for (const auto& sf: *mp) {
            const auto ampi = antigen_to_matching.emplace(sf.found_no, sf);
            if (!ampi.second && ampi.first->second.score < sf.score) {        // already present, replace if sf has hi higher score
                ampi.first->second = sf;
            }
        }
    }

    bool matched = false;
    for (const auto& e: antigen_to_matching) {
        const auto name = found[e.first]->full_name();
        entry.seqs()[e.second.seq_no].add_hi_name(name);
        matched = true;
        if (aReport == seqdb::report::yes)
            report_stream << "    +" << e.second.seq_no << " " << name << '\n';
    }
    return matched;
}

// ----------------------------------------------------------------------

static bool match_normal(SeqdbEntry& entry, const Found& found, const Matching& matching, seqdb::report aReport, std::ostream& report_stream)
{
    bool matched = false;
    if (matching.size() == 1) {
        for (const auto& ms: matching[0]) {
            if (ms.score == matching[0][0].score) {
                const auto name = found[ms.found_no]->full_name();
                entry.seqs()[0].add_hi_name(name);
                matched = true;
                if (aReport == seqdb::report::yes)
                    report_stream << "    + " << name << '\n';
            }
        }
    }
    else {
        if (aReport == seqdb::report::yes) {
            report_found(report_stream, found);
            report_matching(report_stream, matching);
        }
        std::set<size_t> found_assigned;
        for (const auto& m: matching) {
            for (const auto& sf: m) {
                if (sf.score == m[0].score && found_assigned.count(sf.found_no) == 0) {
                    const auto name = found[sf.found_no]->full_name();
                    entry.seqs()[sf.seq_no].add_hi_name(name);
                    matched = true;
                    if (aReport == seqdb::report::yes)
                        report_stream << "    +" << sf.seq_no << " " << name << '\n';
                    found_assigned.insert(sf.found_no);
                }
            }
        }
    }
    return matched;
}

// ----------------------------------------------------------------------

std::vector<std::string> Seqdb::match_hidb(seqdb::report aReport, bool aGreedy)
{
    using namespace std::string_literals;
    std::vector<std::string> not_found_locations;
    std::ostream& report_stream = std::cerr;

    std::vector<const SeqdbEntry*> not_matched;
    for (auto& entry: mEntries) {
        if (aReport == seqdb::report::yes)
            report_stream << '\n' << entry << '\n';

        Found found;
        try {
            find_in_hidb_update_country_lineage_date(found, entry);
        }
        catch (LocationNotFound& err) {
            not_found_locations.push_back(err.what() + " in "s + entry.name());
        }

        if (!found.empty()) {
            Matching matching; // for each seq list of matching [[score, min passage len], found_no] - sorted by score desc
            make_matching(entry, found, matching);
            bool matched = false;
            if (aGreedy)
                matched = match_greedy(entry, found, matching, aReport, report_stream);
            else
                matched = match_normal(entry, found, matching, aReport, report_stream);
            if (!matched) {
                if (aReport == seqdb::report::yes)
                    report_stream << "  <no-passage-matches-in-hidb>?? " << entry.name() << '\n';
                not_matched.push_back(&entry);
            }
        }
        else {
            if (aReport == seqdb::report::yes)
                report_stream << "  <no-name-matches-in-hidb>?? " << entry.name() << '\n';
            not_matched.push_back(&entry);
        }
        if (aReport == seqdb::report::yes)
            report_stream << '\n';
    }

    std::cout << "Matched " << (mEntries.size() - not_matched.size()) << " of " << mEntries.size() << "  " << ((mEntries.size() - not_matched.size()) * 100.0 / mEntries.size()) << '%' << '\n';

    if (aReport == seqdb::report::yes && !not_matched.empty()) {
        report_stream << "Not matched " << not_matched.size() << '\n';
        for (const auto& nm: not_matched) {
            report_stream << "  " << *nm << '\n';
        }
    }
    return not_found_locations;

} // Seqdb::match_hidb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
