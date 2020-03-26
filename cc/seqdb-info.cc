#include <iostream>
#include <string>
using namespace std::string_literals;

#include "acmacs-base/argv.hh"
#include "acmacs-base/string.hh"
#include "seqdb.hh"

// ----------------------------------------------------------------------

struct Info
{
    size_t entries = 0;
    size_t sequences = 0;
    size_t with_hi_names = 0;
    std::map<std::string, size_t> all_by_month;
    std::map<std::string, size_t> with_hi_names_by_month;
    std::map<std::string, size_t> clades;
    std::map<std::string, size_t> clade_set;
};

inline std::ostream& operator<<(std::ostream& s, const Info& c)
{
    s << "Entries:    " << c.entries << '\n' << "Seqs:       " << c.sequences << '\n';
    for (auto [month, count] : c.all_by_month)
        s << "    " << month << ' ' << std::setw(5) << std::right << count << '\n';
    s << "HI matched: " << c.with_hi_names << '\n';
    for (auto [month, count] : c.with_hi_names_by_month)
        s << "    " << month << ' ' << std::setw(5) << std::right << count << '\n';
    s << "By clade:\n";
    for (auto [clade, count] : c.clades)
        s << std::setw(10) << clade << ": " << count << '\n';
    s << "By clade combination:\n";
    for (auto [clade_set, count] : c.clade_set)
        s << std::setw(10) << clade_set << ": " << count << '\n';
    return s;
}

// ----------------------------------------------------------------------

using namespace acmacs::argv;

struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str>   lab{*this, "lab", dflt{"all"}, desc{"filter by lab that submitted the sequence"}};
    argument<str> seqdb_file{*this, arg_name{"~/AD/data/seqdb.json.xz"}, mandatory};
};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        seqdb::setup(opt.seqdb_file, seqdb::report::yes);
        const auto& seqdb = seqdb::get(seqdb::ignore_errors::no, report_time::yes);

        auto update = [](Info& target, const auto& entry) {
            ++target.entries;
            target.sequences += entry.number_of_seqs();
            for (const auto& seq : entry.seqs()) {
                if (!seq.hi_names().empty()) {
                    ++target.with_hi_names;
                    if (const auto date = entry.date(); !date.empty()) {
                        ++target.with_hi_names_by_month[std::string{date.substr(0, 7)}];
                        ++target.with_hi_names_by_month[std::string{date.substr(0, 4)} + "all"];
                    }
                }
                if (const auto date = entry.date(); !date.empty()) {
                    ++target.all_by_month[std::string{date.substr(0, 7)}];
                    ++target.all_by_month[std::string{date.substr(0, 4)} + "all"];
                }
                if (const auto& clades = seq.clades(); !clades.empty()) {
                    for (const auto& clade : clades)
                        ++target.clades[clade];
                    ++target.clade_set[string::join(" ", clades)];
                }
            }
        };

        auto filter_lab = [lab = *opt.lab](const auto& entry) -> bool {
            if (lab.empty() || lab == "all")
                return true;
            return entry.has_lab(lab);
        };

        std::map<std::string, Info> by_virus_type;
        Info total;
        for (const auto& entry : seqdb.entries()) {
            if (filter_lab(entry)) {
                update(total, entry);
                const auto vt = entry.virus_type();
                update(by_virus_type[std::string{vt}], entry);
                if (vt.empty())
                    std::cerr << "No virus_type for " << entry.name() << '\n';
                else if (vt == "B") {
                    const auto lineage = entry.lineage();
                    if (lineage.empty())
                        std::cerr << "No lineage for " << entry.name() << '\n';
                    else
                        update(by_virus_type[std::string{vt} + std::string{lineage}], entry);
                }
            }
        }

        std::cout << "Total:\n" << total << '\n';
        for (auto [virus_type, for_vt] : by_virus_type) {
            std::cout << virus_type << ":\n" << for_vt << '\n';
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
