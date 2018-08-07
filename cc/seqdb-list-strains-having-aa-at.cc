#include <iostream>
#include <cstdlib>
#include <cctype>

#include "acmacs-base/argc-argv.hh"
#include "acmacs-base/stream.hh"
#include "seqdb/seqdb.hh"

using namespace std::string_literals;

// ----------------------------------------------------------------------

constexpr const char* sUsage = "<142N> <283Q> ... [options]\n";

struct PosAA
{
    PosAA(size_t a_pos, char a_aa) : pos(a_pos), aa(a_aa) {}
    size_t pos;
    char aa;
};
inline std::ostream& operator << (std::ostream& out, const PosAA& paa)
{
    return out << paa.pos << paa.aa;
}

int main(int argc, char* const argv[])
{
    try {
        argc_argv args(argc, argv,
                       {
                           {"--db-dir", ""},
                           {"-v", false},
                           {"--verbose", false},
                           {"-h", false},
                           {"--help", false},
                       });
        if (args["-h"] || args["--help"] || args.number_of_arguments() < 1) {
            throw std::runtime_error("Usage: "s + args.program() + sUsage + args.usage_options());
        }
        const bool verbose = args["-v"] || args["--verbose"];

        seqdb::setup_dbs(args["--db-dir"], verbose ? seqdb::report::yes : seqdb::report::no);
        const auto& seqdb = seqdb::get();

        std::vector<PosAA> aa_at;
        for (size_t arg_no = 0; arg_no < args.number_of_arguments(); ++arg_no) {
            const std::string arg = args[arg_no];
            aa_at.emplace_back(std::strtoul(arg.c_str(), nullptr, 10), std::toupper(arg[arg.size() - 1]));
        }
        // std::cout << aa_at << '\n';

        for (const auto entry_seq : seqdb) {
            try {
                const auto aa = entry_seq.seq().amino_acids(true);
                if (all_of(aa_at.begin(), aa_at.end(), [&aa](const auto& paa) -> bool { return aa[paa.pos - 1] == paa.aa; }))
                    std::cout << entry_seq.make_name() << '\n';
            }
            catch (seqdb::SequenceNotAligned&) {
            }
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
