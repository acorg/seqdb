#include <iostream>

#include "acmacs-base/timeit.hh"
#include "seqdb.hh"
using namespace seqdb;
using namespace std::string_literals;

// ----------------------------------------------------------------------

int main(int argc, char* const argv[])
{
    try {
        if (argc != 3)
            throw std::runtime_error("Usage: "s + argv[0] + " seqdb.json.xz <YYYY-MM>");

        const std::string date_to_match{argv[2]};
        if (date_to_match.size() != 4 && date_to_match.size() != 7 && date_to_match.size() != 10)
            throw std::runtime_error{"Invalid date format, expected either of YYYY, YYYY-MM, YYYY-MM-DD"};

        Timeit timeit("seqdb loading ");
        Seqdb seqdb;
        seqdb.load(argv[1]);
        timeit.report();

        for (const auto entry: seqdb) {
            const auto seq_date = entry.entry().date();
            if (seq_date.substr(0, date_to_match.size()) == date_to_match)
                std::cout << string::strip(fmt::format("{} {}", entry.entry().name(), entry.seq().passage())) << ' ' << entry.seq().hi_names() << ' ' << seq_date << '\n';
        }
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
