#include <string>
#include <fstream>

#include "acmacs-base/argv.hh"
#include "acmacs-base/enumerate.hh"
#include "acmacs-base/string.hh"
#include "acmacs-chart-2/factory-import.hh"
#include "acmacs-chart-2/chart.hh"
#include "seqdb/seqdb.hh"

using namespace std::string_literals;

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str>  db_dir{*this, "db-dir"};
    option<bool> with_known_clade{*this, "with-known-clade"};

    argument<str> chart{*this, arg_name{"chart"}, mandatory};
};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);
        seqdb::setup_dbs(opt.db_dir, seqdb::report::no);
        const auto& seqdb = seqdb::get();
        auto chart = acmacs::chart::import_from_file(opt.chart, acmacs::chart::Verify::None, do_report_time(false));
        auto antigens = chart->antigens();
        const auto per_antigen = seqdb.match(*antigens, chart->info()->virus_type());
        for (auto [ag_no, entry] : acmacs::enumerate(per_antigen)) {
            if (entry) {
                if (!entry.seq().hi_name_present(antigens->at(ag_no)->full_name()))
                    throw std::runtime_error(fmt::format("ERROR: internal: matched sequence {} has no matched HI name for {}", entry.entry().name(), antigens->at(ag_no)->full_name()));
                if (opt.with_known_clade) {
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
