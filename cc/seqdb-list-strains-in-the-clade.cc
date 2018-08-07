#include <iostream>

#include "acmacs-base/argc-argv.hh"
#include "seqdb/seqdb.hh"

using namespace std::string_literals;

// ----------------------------------------------------------------------

constexpr const char* sUsage = "<clade> [options]\n";

int main(int argc, char* const argv[])
{
    try {
        argc_argv args(argc, argv, {
                {"--db-dir", ""},
                {"-v", false},
                {"--verbose", false},
                {"-h", false},
                {"--help", false},
            });
        if (args["-h"] || args["--help"] || args.number_of_arguments() != 1) {
            throw std::runtime_error("Usage: "s + args.program() + sUsage + args.usage_options());
        }
        const bool verbose = args["-v"] || args["--verbose"];
        const std::string clade = args[0];
        seqdb::setup_dbs(args["--db-dir"], verbose ? seqdb::report::yes : seqdb::report::no);
        const auto& seqdb = seqdb::get();
        for (const auto entry_seq : seqdb) {
            if (entry_seq.seq().has_clade(clade)) {
                std::cout << entry_seq.make_name() << '\n';
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
