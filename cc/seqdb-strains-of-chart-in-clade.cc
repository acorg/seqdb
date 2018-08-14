#include <iostream>

#include "acmacs-base/argc-argv.hh"
#include "acmacs-base/enumerate.hh"
#include "acmacs-chart-2/factory-import.hh"
#include "acmacs-chart-2/chart.hh"
#include "seqdb/seqdb.hh"

// ----------------------------------------------------------------------

constexpr const char* sUsage = "<clade> <chart> [options]\n";

int main(int argc, char* const argv[])
{
    using namespace std::string_literals;
    try {
        argc_argv args(argc, argv, {
                {"--db-dir", ""},
                {"-v", false},
                {"--verbose", false},
                {"-h", false},
                {"--help", false},
            });
        if (args["-h"] || args["--help"] || args.number_of_arguments() != 2) {
            throw std::runtime_error("Usage: "s + args.program() + sUsage + args.usage_options());
        }
        const bool verbose = args["-v"] || args["--verbose"];
        const std::string clade = args[0];
        seqdb::setup_dbs(args["--db-dir"], verbose ? seqdb::report::yes : seqdb::report::no);
        const auto& seqdb = seqdb::get();
        auto chart = acmacs::chart::import_from_file(args[1], acmacs::chart::Verify::None, report_time::No);
        auto antigens = chart->antigens();
        const auto per_antigen = seqdb.match(*antigens, chart->info()->virus_type());
        for (auto [ag_no, entry] : acmacs::enumerate(per_antigen)) {
            if (entry && entry.seq().has_clade(clade)) {
                if (!entry.seq().hi_name_present(antigens->at(ag_no)->full_name()))
                    throw std::runtime_error("ERROR: internal: matched sequence " + entry.entry().name() + " has no matched HI name for " + antigens->at(ag_no)->full_name());
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
