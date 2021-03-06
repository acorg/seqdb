#include <string>
#include <fstream>

#include "acmacs-base/argv.hh"
#include "acmacs-base/csv.hh"
#include "acmacs-base/enumerate.hh"
#include "acmacs-base/stream.hh"
#include "acmacs-base/string.hh"
#include "acmacs-chart-2/factory-import.hh"
#include "acmacs-chart-2/chart.hh"
#include "seqdb/seqdb.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str>  only_clade{*this, "clade", desc{"report antigens/sera of that clade only"}};
    option<bool> sera_only{*this, "sera-only"};
    option<bool> antigens_only{*this, "antigens-only"};
    option<bool> indexes_only{*this, "indexes-only", desc{"requires --clade"}};
    option<bool> csv{*this, "csv"};
    option<bool> no_unknown{*this, "no-unknown"};
    option<bool> gly{*this, "csv"};
    option<str>  db_dir{*this, "db-dir"};
    option<bool> report_time{*this, "time", desc{"report time of loading chart"}};
    option<bool> verbose{*this, 'v', "verbose"};

    argument<str> chart{*this, arg_name{"chart"}, mandatory};
};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);
        seqdb::setup_dbs(opt.db_dir, opt.verbose ? seqdb::report::yes : seqdb::report::no);
        const auto& seqdb = seqdb::get();
        auto chart = acmacs::chart::import_from_file(opt.chart, acmacs::chart::Verify::None, do_report_time(opt.report_time));
        auto antigens = chart->antigens();
        auto sera = chart->sera();
        auto layout = chart->number_of_projections() > 0 ? chart->projection(0)->layout() : std::shared_ptr<acmacs::Layout>{};
        acmacs::CsvWriter csv_writer;

        const auto write_name = [&csv_writer, &layout, csv = *opt.csv](const char* ag_sr, auto ag_no, auto antigen) {
            if (csv) {
                csv_writer << ag_sr << ag_no << antigen->full_name();
                // csv_writer.field(acmacs::to_string(ag_no));
                // csv_writer.field(antigen->full_name());
            }
            else {
                std::cout << ag_sr << ' ' << ag_no << ' ' << antigen->full_name();
                if (layout && !layout->point_has_coordinates(ag_no))
                    std::cout << " <not-shown-on-map>";
            }
        };
        const auto write_clade = [&csv_writer, csv = *opt.csv](auto clade) {
            if (csv)
                csv_writer.field(clade);
            else
                std::cout << ' ' << clade;
        };
        const auto endl = [&csv_writer, csv = *opt.csv]() {
            if (csv)
                csv_writer.new_row();
            else
                std::cout << '\n';
        };
        const auto show = [&](const auto& ag_sr, const char* prefix) {
            std::vector<size_t> indexes;
            for (auto [ag_no, antigen] : acmacs::enumerate(ag_sr)) {
                bool new_row = false;
                const auto clades = seqdb.clades_for_name(antigen->name());
                if (opt.only_clade.has_value()) {
                    if (std::find(std::begin(clades), std::end(clades), *opt.only_clade) != std::end(clades)) {
                        indexes.push_back(ag_no);
                        if (!opt.indexes_only) {
                            write_name(prefix, ag_no, antigen);
                            new_row = true;
                        }
                    }
                }
                else {
                    write_name(prefix, ag_no, antigen);
                    for (const auto& clade : clades) {
                        if (opt.gly || (clade != "GLY" && clade != "NO-GLY"))
                            write_clade(clade);
                    }
                    new_row = true;
                }
                if (new_row)
                    endl();
            }
            if (!indexes.empty()) {
                if (!opt.indexes_only)
                    std::cout << prefix << " (" << indexes.size() << ") ";
                std::cout << string::join(",", indexes.begin(), indexes.end()) << '\n';
            }
        };

        if (!opt.sera_only)
            show(*antigens, "AG");
        if (!opt.antigens_only)
            show(*sera, "SR");

        if (opt.csv)
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
