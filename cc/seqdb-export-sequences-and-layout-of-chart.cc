#include <string>
#include <fstream>

#include "acmacs-base/argc-argv.hh"
#include "acmacs-base/string.hh"
#include "acmacs-base/read-file.hh"
#include "acmacs-base/csv.hh"
#include "acmacs-chart-2/factory-import.hh"
#include "acmacs-chart-2/chart.hh"
#include "seqdb/seqdb.hh"

using namespace std::string_literals;

// ----------------------------------------------------------------------

constexpr const char* sUsage = " [options] <chart> <output.csv>\n";

int main(int argc, char* const argv[])
{
    try {
        argc_argv args(argc, argv, {
                {"--db-dir", ""},
                {"--amino-acids", false},
                {"--aligned", false},
                {"--projection", 0L},
                {"-v", false},
                {"--verbose", false},
                {"-h", false},
                {"--help", false},
            });
        if (args["-h"] || args["--help"] || args.number_of_arguments() != 2) {
            throw std::runtime_error("Usage: "s + args.program() + sUsage + args.usage_options());
        }
        const bool verbose = args["-v"] || args["--verbose"];
        const bool aligned = args["--aligned"];
        const bool amino_acids = args["--amino-acids"];
        seqdb::setup_dbs(args["--db-dir"], verbose);
        const auto& seqdb = seqdb::get();
        auto chart = acmacs::chart::import_factory(args[0], acmacs::chart::Verify::None);
        auto antigens = chart->antigens();
        const auto per_antigen = seqdb.match(*antigens, chart->info()->virus_type());

        acmacs::CsvWriter writer;
        auto layout = chart->projection(args["--projection"])->layout();
        const auto number_of_dimensions = layout->number_of_dimensions();
        size_t have_sequences = 0;
        for (size_t ag_no = 0; ag_no < antigens->size(); ++ag_no) {
            writer.add_field((*antigens)[ag_no]->full_name());
            for (size_t dim = 0; dim < number_of_dimensions; ++dim)
                writer.add_field(acmacs::to_string(layout->coordinate(ag_no, dim)));
            const auto& entry = per_antigen[ag_no];
            if (entry) {
                ++have_sequences;
                if (amino_acids)
                    writer.add_field(entry.seq().amino_acids(aligned));
                else
                    writer.add_field(entry.seq().nucleotides(aligned));
            }
            else {
                writer.add_empty_field();
            }
            writer.new_row();
        }
        acmacs::file::write(args[1], writer);
        std::cout << "Antigens written:  " << antigens->size() << '\n';
        std::cout << "Sequences written: " << have_sequences << '\n';
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