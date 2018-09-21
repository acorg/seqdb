#include <string>
#include <fstream>

#include "acmacs-base/argc-argv.hh"
#include "acmacs-base/enumerate.hh"
#include "acmacs-base/stream.hh"
#include "acmacs-base/string.hh"
#include "acmacs-chart-2/factory-import.hh"
#include "acmacs-chart-2/chart.hh"
#include "seqdb/seqdb.hh"

// ----------------------------------------------------------------------

constexpr const char* sUsage = " [options] <chart>\n";

int main(int argc, char* const argv[])
{
    try {
        using namespace std::string_literals;
        argc_argv args(argc, argv,
                       {
                           {"--chart-name", false},
                           {"--no-gly", false},
                           {"--no-unknown", false},
                           {"--clade-left", false},
                           {"--db-dir", ""},
                           {"--time", false, "report time of loading chart"},
                           {"-v", false},
                           {"--verbose", false},
                           {"-h", false},
                           {"--help", false},
                       });
        if (args["-h"] || args["--help"] || args.number_of_arguments() != 1) {
            throw std::runtime_error("Usage: "s + args.program() + sUsage + args.usage_options());
        }
        const bool verbose = args["-v"] || args["--verbose"];
        seqdb::setup_dbs(args["--db-dir"], verbose ? seqdb::report::yes : seqdb::report::no);
        const auto& seqdb = seqdb::get();
        auto chart = acmacs::chart::import_from_file(args[0], acmacs::chart::Verify::None, args["--time"] ? report_time::Yes : report_time::No);
        if (args["--chart-name"])
            std::cout << chart->make_name() << '\n';
        auto sera = chart->sera();
        chart->set_homologous(acmacs::chart::Chart::find_homologous_for_big_chart::yes, sera);
        auto antigens = chart->antigens();
        const auto per_antigen = seqdb.match(*antigens, chart->info()->virus_type());
        for (auto [sr_no, serum] : acmacs::enumerate(*sera)) {
            std::set<std::string> clades;
            for (auto ag_no : serum->homologous_antigens()) {
                if (const auto& entry_seq = per_antigen[ag_no]; entry_seq) {
                    for (const auto& clade : entry_seq.seq().clades()) {
                        if (!args["--no-gly"] || (clade != "GLY" && clade != "NO-GLY"))
                            clades.insert(clade);
                    }
                }
            }
            if (!args["--no-unknown"] || !clades.empty()) {
                if (args["--clade-left"]) {
                    std::cout << std::setw(20) << std::left << string::join(" ", clades.begin(), clades.end()) << ' ' << serum->full_name() << '\n';
                }
                else {
                    std::cout << "SR " << std::setw(5) << sr_no << ' ' << serum->full_name() << ' ' << clades << '\n';
                }
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
