#include <iostream>

#include "acmacs-base/argc-argv.hh"
#include "seqdb.hh"

using namespace std::string_literals;

// ----------------------------------------------------------------------

constexpr const char* sUsage = " [options] <clade-name, e.g. DEL2017 or all>\n";

int main(int argc, char* const argv[])
{
    try {
        argc_argv args(argc, argv, {
                {"--db-dir", ""},
                {"--lab", ""},
                {"--flu", "", "A(H1N1), A(H3N2), B"},
                {"-v", false},
                {"--verbose", false},
                {"-h", false},
                {"--help", false},
            });
        if (args["-h"] || args["--help"] || args.number_of_arguments() != 1) {
            throw std::runtime_error("Usage: "s + args.program() + sUsage + args.usage_options());
        }
        const bool verbose = args["-v"] || args["--verbose"];
        seqdb::setup_dbs(std::string(args["--db-dir"]), verbose ? seqdb::report::yes : seqdb::report::no);
        if (args[0] == "all"s) {
            for (const auto entry_seq: seqdb::get()) {
                if ((!args["--lab"] || entry_seq.seq().has_lab(args["--lab"].str())) && (!args["--flu"] || entry_seq.entry().virus_type() == args["--flu"].str_view()))
                    std::cout << std::setw(60) << std::left << entry_seq.make_name() << '\t' << entry_seq.seq().clades()
                              << '\t' << entry_seq.entry().virus_type() << '\t' << entry_seq.entry().lineage() << '\t' << entry_seq.entry().dates()
                              << '\t' << entry_seq.seq().lab() << '\n';
            }
        }
        else {
            for (const auto entry_seq: seqdb::get()) {
                if (entry_seq.seq().has_clade(std::string(args[0])) && (!args["--lab"] || entry_seq.seq().has_lab(args["--lab"].str())) && (!args["--flu"] || entry_seq.entry().virus_type() == args["--flu"].str_view()))
                    std::cout << std::setw(60) << std::left << entry_seq.make_name()
                              << '\t' << entry_seq.entry().virus_type() << '\t' << entry_seq.entry().lineage() << '\t' << entry_seq.entry().dates()
                              << '\t' << entry_seq.seq().lab() << '\n';
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
