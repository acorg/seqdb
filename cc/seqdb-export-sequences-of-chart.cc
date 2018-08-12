#include <string>
#include <fstream>

#include "acmacs-base/argc-argv.hh"
#include "acmacs-base/enumerate.hh"
#include "acmacs-base/string.hh"
#include "acmacs-base/read-file.hh"
#include "acmacs-chart-2/factory-import.hh"
#include "acmacs-chart-2/chart.hh"
#include "seqdb/seqdb.hh"

using namespace std::string_literals;

// ----------------------------------------------------------------------

constexpr const char* sUsage = " [options] <chart> <output.fasta>\n";

int main(int argc, char* const argv[])
{
    try {
        argc_argv args(argc, argv, {
                {"--db-dir", ""},
                {"--amino-acids", false},
                {"--aligned", false},
                {"--replace-spaces-in-names", false},
                {"--name-from-chart", false},
                {"--time", false, "report time of loading chart"},
                {"-v", false},
                {"--verbose", false},
                {"-h", false},
                {"--help", false},
            });
        if (args["-h"] || args["--help"] || args.number_of_arguments() != 2) {
            throw std::runtime_error("Usage: "s + args.program() + sUsage + args.usage_options());
        }
        const bool verbose = args["-v"] || args["--verbose"];
        seqdb::setup_dbs(args["--db-dir"], verbose ? seqdb::report::yes : seqdb::report::no);
        const auto& seqdb = seqdb::get();
        auto chart = acmacs::chart::import_from_file(args[0], acmacs::chart::Verify::None, args["--time"] ? report_time::Yes : report_time::No);
        auto antigens = chart->antigens();
        const auto per_antigen = seqdb.match(*antigens, chart->info()->virus_type());
        std::string output;
        for (auto [ag_no, entry] : acmacs::enumerate(per_antigen)) {
            if (entry) {
                if (!entry.seq().hi_name_present(antigens->at(ag_no)->full_name()))
                    throw std::runtime_error("ERROR: internal: matched sequence " + entry.entry().name() + " has no matched HI name for " + antigens->at(ag_no)->full_name());
                auto name = entry.make_name();
                if (args["--replace-spaces-in-names"])
                    name = string::replace(name, " ", "_");
                if (args["--name-from-chart"])
                    name = antigens->at(ag_no)->full_name() + " ==> " + name;
                output += ">" + name + "\n";
                if (args["--amino-acids"])
                    output += entry.seq().amino_acids(args["--aligned"]);
                else
                    output += entry.seq().nucleotides(args["--aligned"]);
                output += "\n";
            }
        }
        acmacs::file::write(args[1], output);
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
