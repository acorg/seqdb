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

constexpr const char* sUsage = " [options] <chart> <output.json>\n";

int main(int argc, char* const argv[])
{
    try {
        argc_argv args(argc, argv,
                       {
                           {"--db-dir", ""},
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
        seqdb::setup_dbs(args["--db-dir"].str(), verbose ? seqdb::report::yes : seqdb::report::no);
        auto chart = acmacs::chart::import_from_file(args[0], acmacs::chart::Verify::None, do_report_time(args["--time"]));
        const auto output = seqdb::sequences_of_chart_for_ace_view_1(*chart);
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
