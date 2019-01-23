#include <string>
#include <fstream>

#include "acmacs-base/argv.hh"
#include "acmacs-base/enumerate.hh"
#include "acmacs-base/string.hh"
#include "acmacs-base/csv.hh"
#include "acmacs-base/read-file.hh"
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

    argument<str> chart{*this, arg_name{"chart"}, mandatory};
    argument<str> output{*this, arg_name{"chart-seqids.csv"}, mandatory};
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
        acmacs::CsvWriter csv;
        csv << "No" << "Name" << "Reassortant" << "Annotations" << "Passage" << "Seq Id" << acmacs::CsvWriter::end_of_row;
        for (auto [ag_no, entry] : acmacs::enumerate(per_antigen)) {
            if (entry) {
                if (!entry.seq().hi_name_present(antigens->at(ag_no)->full_name()))
                    std::cerr << "WARNING: matched sequence " << entry.seq_id(seqdb::SeqdbEntrySeq::encoded_t::yes) << " has no matched HI name for " << antigens->at(ag_no)->full_name() << '\n';
                auto antigen = antigens->at(ag_no);
                csv << ag_no << antigen->name() << antigen->reassortant() << string::join(" ", antigen->annotations()) << antigen->passage()
                    << entry.seq_id(seqdb::SeqdbEntrySeq::encoded_t::yes) << acmacs::CsvWriter::end_of_row;
            }
        }
        acmacs::file::write(opt.output, csv);
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
