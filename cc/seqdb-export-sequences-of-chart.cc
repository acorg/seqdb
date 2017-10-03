#include <string>
#include <fstream>

#include "acmacs-base/argc-argv.hh"
#include "acmacs-base/string.hh"
#include "acmacs-base/read-file.hh"
#include "acmacs-chart/ace.hh"
#include "seqdb.hh"

using namespace std::string_literals;

// ----------------------------------------------------------------------

constexpr const char* sUsage = " [options] <chart.ace> <output.fasta>\n";

int main(int argc, char* const argv[])
{
    try {
        argc_argv args(argc, argv, {
                {"--db-dir", ""},
                {"--amino-acids", false},
                {"--aligned", false},
                {"--replace-spaces-in-names", false},
                {"-v", false},
                {"--verbose", false},
                {"-h", false},
                {"--help", false},
            });
        if (args["-h"] || args["--help"] || args.number_of_arguments() != 2) {
            throw std::runtime_error("Usage: "s + args.program() + sUsage + args.usage_options());
        }
        // const bool verbose = args["-v"] || args["--verbose"];
        seqdb::setup_dbs(args["--db-dir"]);
        const auto& seqdb = seqdb::get();
        std::unique_ptr<Chart> chart{import_chart(args[0])};
        const auto per_antigen = seqdb.match(chart->antigens());
        std::string output;
        for (const auto& entry: per_antigen) {
            if (entry) {
                auto name = entry.make_name();
                if (args["--replace-spaces-in-names"])
                    name = string::replace(name, " ", "_");
                output += ">" + name + "\n";
                if (args["--amino-acids"])
                    output += entry.seq().amino_acids(args["--aligned"]);
                else
                    output += entry.seq().nucleotides(args["--aligned"]);
                output += "\n";
            }
        }
        acmacs_base::write_file(args[1], output);
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
