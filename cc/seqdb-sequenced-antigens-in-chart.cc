#include <string>
#include <fstream>

#include "acmacs-base/argc-argv.hh"
#include "acmacs-base/enumerate.hh"
#include "acmacs-base/string.hh"
#include "acmacs-chart-2/factory-import.hh"
#include "acmacs-chart-2/chart.hh"
#include "seqdb/seqdb.hh"

using namespace std::string_literals;

// ----------------------------------------------------------------------

constexpr const char* sUsage = " [options] <chart>\n";

int main(int argc, char* const argv[])
{
    try {
        argc_argv args(argc, argv, {
                {"--db-dir", ""},
                {"--with-known-clade", false, "GLY and NO-GLY are not clades"},
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
        auto antigens = chart->antigens();
        const auto per_antigen = seqdb.match(*antigens, chart->info()->virus_type());
        for (auto [ag_no, entry] : acmacs::enumerate(per_antigen)) {
            if (entry) {
                if (!entry.seq().hi_name_present(antigens->at(ag_no)->full_name()))
                    throw std::runtime_error("ERROR: internal: matched sequence " + entry.entry().name() + " has no matched HI name for " + antigens->at(ag_no)->full_name());
                if (args["--with-known-clade"]) {
                    auto clades = entry.seq().clades(); // copy!
                    clades.erase(std::find_if(clades.begin(), clades.end(), [](const auto& clade) -> bool { return clade == "GLY" || clade == "NO-GLY"; }), clades.end());
                    if (clades.empty())
                        continue;
                }
                std::cout << std::setw(5) << ag_no << ' ' << antigens->at(ag_no)->full_name() << ' ' << entry.seq().clades() << '\n';
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
