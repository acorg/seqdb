#include <string>
#include <fstream>

#include "acmacs-base/argc-argv.hh"
#include "acmacs-base/csv.hh"
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
                           {"--csv", false},
                           {"--gly", false},
                           {"--no-unknown", false},
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
        seqdb::setup_dbs(args["--db-dir"].str(), verbose ? seqdb::report::yes : seqdb::report::no);
        const auto& seqdb = seqdb::get();
        auto chart = acmacs::chart::import_from_file(args[0], acmacs::chart::Verify::None, do_report_time(args["--time"]));
        auto antigens = chart->antigens();
        const auto per_antigen = seqdb.match(*antigens, chart->info()->virus_type());
        const bool csv = args["--csv"];
        acmacs::CsvWriter csv_writer;
        for (auto [ag_no, antigen] : acmacs::enumerate(*antigens)) {
            const auto write_name = [&csv_writer,csv](auto ag_no, auto antigen) {
                if (csv) {
                    csv_writer.field(acmacs::to_string(ag_no));
                    csv_writer.field(antigen->full_name());
                }
                else
                    std::cout << ag_no << ' ' << antigen->full_name();
            };
            if (const auto& entry_seq = per_antigen[ag_no]; entry_seq) {
                write_name(ag_no, antigen);
                for (const auto& clade : entry_seq.seq().clades()) {
                    if (args["--gly"] || (clade != "GLY" && clade != "NO-GLY")) {
                        if (csv)
                            csv_writer.field(clade);
                        else
                            std::cout << ' ' << clade;
                    }
                }
            }
            else if (!args["--no-unknown"])
                write_name(ag_no, antigen);
            if (csv)
                csv_writer.new_row();
            else
                std::cout << '\n';
        }
        if (csv)
            std::cout << static_cast<std::string_view>(csv_writer) << '\n';
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
