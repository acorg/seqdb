#include <iostream>
#include <algorithm>
#include <array>
#include <vector>

#include "acmacs-base/argv.hh"
#include "acmacs-base/string.hh"
#include "acmacs-base/normalize.hh"
#include "seqdb.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str>   flu{*this, "flu", desc{"filter by virus type: H1, H3, B"}};
    option<str>   lineage{*this, "lineage", desc{"filter by lineage"}};
    option<str>   clade{*this, "clade", desc{"filter by clade"}};
    option<str>   start_date{*this, "start-date", desc{"isolated on or after date: YYYY-MM-DD"}};
    option<str>   end_date{*this, "end-date", desc{"isolated before date: YYYY-MM-DD"}};
};

constexpr const size_t NUM_POSITIONS = 600;
// constexpr const size_t FIRST_AA = static_cast<size_t>('A');
// constexpr const size_t NUM_AAS = static_cast<size_t>('Z') + 1 - FIRST_AA;

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        // seqdb::setup(opt.seqdb_file, seqdb::report::yes);
        const auto& seqdb = seqdb::get(seqdb::ignore_errors::no, report_time::yes);
        size_t num_sequences = 0;
        auto iter = seqdb.begin();
        iter.filter_subtype(acmacs::normalize_virus_type(*opt.flu))
                .filter_lineage(acmacs::normalize_lineage(*opt.lineage))
                .filter_clade(string::upper(*opt.clade))
                .filter_date_range(opt.start_date, opt.end_date)
                ;

        std::array<std::vector<std::pair<char, unsigned>>, NUM_POSITIONS> amino_acids_data;
        for (; iter != seqdb.end(); ++iter) {
            const auto seq = (*iter).seq().amino_acids(true);
            for (size_t pos = 0; pos < seq.size(); ++pos) {
                auto& at = amino_acids_data[pos];
                if (const auto found = std::find_if(std::begin(at), std::end(at), [aa=seq[pos]](const auto& entry) { return entry.first == aa; }); found != std::end(at))
                    ++found->second;
                else
                    at.emplace_back(seq[pos], 1);
            }
            ++num_sequences;
        }
        std::cout << "sequences: " << num_sequences << '\n';

        for (size_t pos = 0; pos < amino_acids_data.size(); ++pos) {
            if (amino_acids_data[pos].size() > 1) {
                std::cout << std::setw(3) << (pos + 1);
                for (const auto& [aa, num] : amino_acids_data[pos])
                    std::cout << ' ' << aa << ':' << num;
                std::cout << '\n';
            }
        }
        
        // std::array<std::array<unsigned, NUM_AAS>, NUM_POSITIONS> amino_acids_data; // amino_acids_data[pos][aa]
        // std::for_each(std::begin(amino_acids_data), std::end(amino_acids_data), [](auto& entry) { entry.fill(0); }); // array requires initialization!

        // for (; iter != seqdb.end(); ++iter) {
        //     const auto seq = (*iter).seq().amino_acids(true);
        //     for (size_t pos = 0; pos < seq.size(); ++pos)
        //         ++amino_acids_data[pos][static_cast<size_t>(seq[pos]) - FIRST_AA];
        //     ++num_sequences;
        // }
        // std::cout << "sequences: " << num_sequences << '\n';

        // for (size_t pos = 0; pos < amino_acids_data.size(); ++pos) {
        //     std::cout << std::setw(3) << pos;
        //     for (size_t aa = 0; aa < NUM_AAS; ++aa) {
        //         if (amino_acids_data[pos][aa])
        //             std::cout << ' ' << static_cast<char>(aa + FIRST_AA + 1) << ':' << amino_acids_data[pos][aa];
        //     }
        //     std::cout << '\n';
        // }
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
