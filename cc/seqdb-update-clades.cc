#include <iostream>

#include "acmacs-base/argc-argv.hh"
#include "seqdb.hh"

using namespace std::string_literals;

// ----------------------------------------------------------------------

constexpr const char* sUsage = " [options]\n";

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
        if (args["-h"] || args["--help"] || args.number_of_arguments() != 0) {
            throw std::runtime_error("Usage: "s + args.program() + sUsage + args.usage_options());
        }
        const bool verbose = args["-v"] || args["--verbose"];
        seqdb::setup_dbs(args["--db-dir"], verbose ? seqdb::report::yes : seqdb::report::no);
        auto& seqdb = seqdb::get_for_updating();
        seqdb.update_clades(verbose ? seqdb::report::yes : seqdb::report::no);
        seqdb.save(std::string{}, 2);
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
