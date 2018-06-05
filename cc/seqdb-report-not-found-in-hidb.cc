#include <iostream>

#include "seqdb.hh"
using namespace seqdb;
using namespace std::string_literals;

// ----------------------------------------------------------------------

int main(int argc, char* const argv[])
{
    try {
        if (argc != 2)
            throw std::runtime_error("Usage: "s + argv[0] + " seqdb.json.xz");

        seqdb::setup(argv[1], seqdb::report::yes);
        const auto& seqdb = seqdb::get(seqdb::ignore_errors::no, report_time::Yes);

        size_t not_matched = 0;
        for (const auto& entry: seqdb.entries()) {
            size_t hi_matched = 0;
            for (const auto& seq: entry.seqs()) {
                if (!seq.hi_names().empty())
                    ++hi_matched;
            }
            if (hi_matched == 0) {
                std::cout << entry.name() << ' ' << entry.dates() << ' ' << entry.lab_ids() << '\n';
                ++not_matched;
            }
            else if (hi_matched != entry.seqs().size())
                std::cout << "   partially matched: " << entry.name() << '\n';
        }
        std::cout << "Total       entries: " << seqdb.entries().size() << '\n';
        std::cout << "Not matched entries: " << not_matched << '\n';
        return 0;
    }
    catch (std::exception& err) {
        std::cerr << "ERROR: " << err.what() << '\n';
        return 1;
    }
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
