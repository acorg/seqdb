#include <string>
#include <fstream>

#include "acmacs-base/argc-argv.hh"
#include "acmacs-chart/ace.hh"
#include "seqdb.hh"

using namespace std::string_literals;

// ----------------------------------------------------------------------

int main(int argc, char* const argv[])
{
    try {
        argc_argv args(argc, argv, {"--seqdb"});
        if (args["-h"] || args["--help"] || args.number_of_arguments() != 2)
            throw std::runtime_error("Usage: "s + args.program() + " [--seqdb <seqdb.json.xz>] [--amino-acids] [--aligned] <chart.ace> <output.fasta>");
        const auto& seqdb = seqdb::get(args.get("--seqdb", "/Users/eu/AD/data/seqdb.json.xz"));
        std::unique_ptr<Chart> chart{import_chart(args[0])};
        const auto per_antigen = seqdb.match(chart->antigens());
        std::string output;
        for (const auto& entry: per_antigen) {
            if (entry) {
                output += ">" + entry.make_name() + "\n";
                if (args["--amino-acids"])
                    output += entry.seq().amino_acids(args["--aligned"]);
                else
                    output += entry.seq().nucleotides(args["--aligned"]);
                output += "\n";
            }
        }
        std::ofstream{args[1]} << output;
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
