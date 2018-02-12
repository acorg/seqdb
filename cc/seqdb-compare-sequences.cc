#include <iostream>

#include "acmacs-base/argc-argv.hh"
#include "acmacs-base/range.hh"
#include "acmacs-base/enumerate.hh"
#include "seqdb.hh"

using namespace std::string_literals;

// ----------------------------------------------------------------------

class Comparer
{
 public:
    Comparer() = default;

    void add(const seqdb::SeqdbEntrySeq& entry_seq);
    void analyse();
    void report() const;

 private:
    struct Entry
    {
        std::string name;
        std::string sequence;
    };

    struct EntryPos
    {
        std::map<char, size_t> aa_count;
    };

    std::vector<Entry> entries_;
    std::vector<EntryPos> entries_pos_{600};

}; //  class Comparer

// ----------------------------------------------------------------------

int main(int argc, char* const argv[])
{
    try {
        argc_argv args(argc, argv, {
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
        seqdb::setup_dbs(args["--db-dir"], verbose);
        Comparer comparer;
        for (auto arg_no : acmacs::range(args.number_of_arguments())) {
            if (const auto* entry_seq = seqdb::get().find_hi_name(args[arg_no]); entry_seq && entry_seq->seq().aligned()) {
                comparer.add(*entry_seq);
                // std::cout << std::setw(60) << std::left << entry_seq->make_name() << entry_seq->seq().amino_acids(true) << '\n';
            }
            else {
                std::cerr << "WARNING: no found or not aligned: " << args[arg_no] << '\n';
            }
        }
        comparer.analyse();
        comparer.report();
        return 0;
    }
    catch (std::exception& err) {
        std::cerr << "ERROR: " << err.what() << '\n';
        return 1;
    }
}

// ----------------------------------------------------------------------

void Comparer::add(const seqdb::SeqdbEntrySeq& entry_seq)
{
    const Entry entry{entry_seq.make_name(), entry_seq.seq().amino_acids(true)};
    entries_.push_back(entry);
    for (const auto [pos, aa] : acmacs::enumerate(entry.sequence)) {
        if (aa != 'X')
            ++entries_pos_[pos].aa_count[aa];
    }

} // Comparer::add

// ----------------------------------------------------------------------

void Comparer::analyse()
{

} // Comparer::analyse

// ----------------------------------------------------------------------

void Comparer::report() const
{
    // const auto max_seq_size = std::accumulate(entries_.begin(), entries_.end(), 0, [](auto mns, const auto& entry) { return std::max(mns, static_cast<decltype(mns)>(entry.sequence.size())); });
    // std::cout << "max_seq_size " << max_seq_size << '\n';

    // const auto max_name_size = std::accumulate(entries_.begin(), entries_.end(), 0, [](auto mns, const auto& entry) { return std::max(mns, static_cast<decltype(mns)>(entry.name.size())); });
    // for (const auto& entry : entries_)
    //     std::cout << std::setw(max_name_size + 2) << std::left << entry.name << entry.sequence << '\n';

    for (const auto [pos, entry_pos] : acmacs::enumerate(entries_pos_, 1)) {
        if (entry_pos.aa_count.size() > 1) {
            const auto max_occur = std::accumulate(entry_pos.aa_count.begin(), entry_pos.aa_count.end(), 0UL, [](auto max_count, const auto& elt) { return std::max(max_count, elt.second); });
            const auto max_occur_percent = max_occur * 100.0 / entries_.size();
            if (max_occur_percent < 95.0)
                std::cout << std::setw(3) << std::right << pos << ' ' << entry_pos.aa_count << ' ' << max_occur_percent << '\n';
        }
    }

} // Comparer::report

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
