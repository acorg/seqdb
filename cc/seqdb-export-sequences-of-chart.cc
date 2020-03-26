#include <string>
#include <fstream>

#include "acmacs-base/argv.hh"
#include "acmacs-base/enumerate.hh"
#include "acmacs-base/string.hh"
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
    option<bool> amino_acids{*this, "amino-acids"};
    option<bool> aligned{*this, "aligned"};
    option<str>  clade{*this, "clade", desc{"only sequences in clade"}};
    option<bool> replace_spaces_in_names{*this, "replace-spaces-in-names"};
    option<bool> name_from_chart{*this, "name-from-chart"};
    option<bool> report_time{*this, "time", desc{"report time of loading chart"}};
    option<bool> verbose{*this, 'v', "verbose"};

    argument<str> chart{*this, arg_name{"chart"}, mandatory};
    argument<str> output_fasta{*this, arg_name{"output.fasta"}, mandatory};
};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);
        seqdb::setup_dbs(*opt.db_dir, opt.verbose ? seqdb::report::yes : seqdb::report::no);
        const auto& seqdb = seqdb::get();
        auto chart = acmacs::chart::import_from_file(opt.chart, acmacs::chart::Verify::None, do_report_time(opt.report_time));
        auto antigens = chart->antigens();
        const auto per_antigen = seqdb.match(*antigens, chart->info()->virus_type());
        std::string output;
        size_t seq_found = 0;
        for (auto [ag_no, entry] : acmacs::enumerate(per_antigen)) {
            if (entry && (!opt.clade || entry.seq().has_clade(opt.clade))) {
                ++seq_found;
                if (!entry.seq().hi_name_present(antigens->at(ag_no)->full_name()))
                    throw std::runtime_error(fmt::format("ERROR: internal: matched sequence {} has no matched HI name for {}", entry.entry().name(), antigens->at(ag_no)->full_name()));
                auto name = entry.make_name();
                if (opt.replace_spaces_in_names)
                    name = string::replace(name, " ", "_");
                if (opt.name_from_chart)
                    name = antigens->at(ag_no)->full_name(); // + " ==> " + name;
                output += ">" + name + "\n";
                try {
                    if (opt.amino_acids)
                        output += entry.seq().amino_acids(opt.aligned);
                    else
                        output += entry.seq().nucleotides(opt.aligned);
                }
                catch (seqdb::SequenceNotAligned& err) {
                    std::cerr << "WARNING: " << err.what() << ' ' << entry.entry().name() << '\n';
                }
                output += "\n";
            }
        }
        std::cerr << "INFO: " << seq_found << " sequences exported\n";
        acmacs::file::write(opt.output_fasta, output);
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
