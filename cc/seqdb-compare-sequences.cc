#include <iostream>

#include "acmacs-base/argc-argv.hh"
#include "acmacs-base/range.hh"
#include "acmacs-base/enumerate.hh"
#include "acmacs-base/stream.hh"
#include "seqdb.hh"

using namespace std::string_literals;

// ----------------------------------------------------------------------

class Comparer
{
 public:
    Comparer() = default;

    void add(const seqdb::SeqdbEntrySeq& entry_seq) { entries_.push_back(entry_seq); }
    void analyse();
    void report(std::ostream& output, bool all_pos, bool vertical) const;
    void report_old() const;

 private:
    // struct Entry
    // {
    //     std::string name;
    //     std::string aa;
    //     std::string nuc;
    // };

    struct EntryPos
    {
        std::map<char, size_t> aa_count;
    };

    std::vector<seqdb::SeqdbEntrySeq> entries_;

}; //  class Comparer

// ----------------------------------------------------------------------

int main(int argc, char* const argv[])
{
    try {
        argc_argv args(argc, argv, {
                {"--all-pos", false, "show all positions"},
                {"--vertical", false, "show nucs vertically"},
                {"--old", false, "report in the old style"},

                {"--db-dir", ""},
                // {"--lab", ""},
                // {"--flu", "", "A(H1N1), A(H3N2), B"},
                {"-v", false},
                {"--verbose", false},
                {"-h", false},
                {"--help", false},
            });
        if (args["-h"] || args["--help"] || args.number_of_arguments() < 2) {
            throw std::runtime_error("Usage: "s + args.program() + " [options] <name> ...\n" + args.usage_options());
        }
        const bool verbose = args["-v"] || args["--verbose"];
        seqdb::setup_dbs(args["--db-dir"], verbose ? seqdb::report::yes : seqdb::report::no);
        Comparer comparer;
        const auto& seqdb = seqdb::get();
        for (auto arg_no : acmacs::range(args.number_of_arguments())) {
            if (const auto* entry_seq_1 = seqdb.find_hi_name(args[arg_no]); entry_seq_1 && entry_seq_1->seq().aligned()) {
                comparer.add(*entry_seq_1);
                // std::cout << std::setw(60) << std::left << entry_seq_1->make_name() << entry_seq_1->seq().amino_acids(true) << '\n';
            }
            else if (const auto entry_seq_2 = seqdb.find_by_seq_id(args[arg_no]); entry_seq_2 && entry_seq_2.seq().aligned()) {
                comparer.add(entry_seq_2);
            }
            else if (const auto* entry = seqdb.find_by_name(args[arg_no]); entry) {
                for (const auto& seq : entry->seqs()) {
                    if (seq.aligned())
                        comparer.add({*entry, seq});
                }
            }
            else {
                std::cerr << "WARNING: no found or not aligned: " << args[arg_no] << '\n';
            }
        }
        comparer.analyse();
        if (args["--old"])
            comparer.report_old();
        else
            comparer.report(std::cout, args["--all-pos"], args["--vertical"]);
        return 0;
    }
    catch (std::exception& err) {
        std::cerr << "ERROR: " << err.what() << '\n';
        return 1;
    }
}

// ----------------------------------------------------------------------

void Comparer::analyse()
{

} // Comparer::analyse

// ----------------------------------------------------------------------

void Comparer::report(std::ostream& output, bool all_pos, bool vertical) const
{
    std::vector<std::string> aas(entries_.size()), nucs(entries_.size());
    size_t max_aa = 0;
    for (const auto [no, entry] : acmacs::enumerate(entries_)) {
        output << std::setw(2) << std::right << no + 1 << ' ' << entry.make_name() << ' ' << entry.entry().lineage() << ' ' << entry.seq().clades() << '\n';
        aas[no] = entry.seq().amino_acids(true);
        nucs[no] = entry.seq().nucleotides(true);
        max_aa = std::max(max_aa, aas[no].size());
    }
    output << '\n';

    if (vertical) {
        output << "pos ";
        for (size_t no = 1; no <= entries_.size(); ++no)
            output << ' ' << std::setw(2) << std::right << no << ' ';
    }
    else {
        output << "pos  ";
        for (size_t no = 1; no <= entries_.size(); ++no)
            output << ' ' << std::setw(2) << std::right << no << "    ";
    }
    output << '\n';

    auto nucs_horizontal = [&output, &aas, &nucs, this](size_t pos, const auto& aa_count, const auto& nuc_count) {
        output << std::setw(3) << std::right << pos + 1 << "  ";
        for (size_t no = 0; no < this->entries_.size(); ++no) {
            if (pos < aas[no].size())
                output << aas[no][pos] << ' ' << std::setw(3) << std::left << nucs[no].substr(pos * 3, 3);
            else
                output << "     ";
            output << "  ";
        }
        if (nuc_count.size() > 1)
            output << nuc_count;
        if (aa_count.size() > 1)
            output << ' ' << aa_count;
        output << '\n';
    };

    auto nucs_vertical = [&output, &aas, &nucs, this](size_t pos, const auto& aa_count, const auto& nuc_count) {
        output << std::setw(3) << std::right << pos + 1 << "  ";
        std::vector<std::string> nucs_at(this->entries_.size());
        for (size_t no = 0; no < this->entries_.size(); ++no) {
            nucs_at[no] = nucs[no].substr(pos * 3, 3);
            if (pos < aas[no].size())
                output << aas[no][pos] << nucs_at[no][0];
            else
                output << "  ";
            output << "  ";
        }
        if (nuc_count.size() > 1)
            output << nuc_count;
        if (aa_count.size() > 1)
            output << ' ' << aa_count;
        for (size_t nn = 1; nn < 3; ++nn) {
            output << "\n     ";
            for (size_t no = 0; no < this->entries_.size(); ++no) {
                output << ' ';
                if (nucs_at[no].size() > nn)
                    output << nucs_at[no][nn];
                else
                    output << ' ';
                output << "  ";
            }
        }
        output << '\n';
    };

    for (size_t pos = 0; pos < max_aa; ++pos) {
        std::map<char, size_t> aa_count;
        std::map<std::string, size_t> nuc_count;
        for (size_t no = 0; no < entries_.size(); ++no) {
            if (pos < aas[no].size()) {
                if (const std::string nucs_at = nucs[no].substr(pos * 3, 3); !nucs_at.empty())
                    ++nuc_count[nucs_at];
                ++aa_count[aas[no][pos]];
            }
        }
        if (all_pos || nuc_count.size() > 1) {
            if (vertical)
                nucs_vertical(pos, aa_count, nuc_count);
            else
                nucs_horizontal(pos, aa_count, nuc_count);
        }
    }

} // Comparer::report

// ----------------------------------------------------------------------

void Comparer::report_old() const
{
    // const auto max_seq_size = std::accumulate(entries_.begin(), entries_.end(), 0, [](auto mns, const auto& entry) { return std::max(mns, static_cast<decltype(mns)>(entry.sequence.size())); });
    // std::cout << "max_seq_size " << max_seq_size << '\n';

    // const auto max_name_size = std::accumulate(entries_.begin(), entries_.end(), 0, [](auto mns, const auto& entry) { return std::max(mns, static_cast<decltype(mns)>(entry.name.size())); });
    // for (const auto& entry : entries_)
    //     std::cout << std::setw(max_name_size + 2) << std::left << entry.name << entry.sequence << '\n';

    std::vector<EntryPos> entries_pos{600};
    for (const auto& entry : entries_) {
        for (const auto [pos, aa] : acmacs::enumerate(entry.seq().amino_acids(true))) {
            if (aa != 'X')
                ++entries_pos[pos].aa_count[aa];
        }
    }

    for (const auto [pos, entry_pos] : acmacs::enumerate(entries_pos, 1)) {
        if (entry_pos.aa_count.size() > 1) {
            const auto max_occur = std::accumulate(entry_pos.aa_count.begin(), entry_pos.aa_count.end(), 0UL, [](auto max_count, const auto& elt) { return std::max(max_count, elt.second); });
            const auto max_occur_percent = max_occur * 100.0 / entries_.size();
            if (max_occur_percent < 95.0)
                std::cout << std::setw(3) << std::right << pos << ' ' << entry_pos.aa_count << ' ' << std::setprecision(3) << max_occur_percent << "%\n";
        }
    }

} // Comparer::report_old

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
