#include <iostream>
#include <string>
using namespace std::string_literals;

#include "acmacs-base/argc-argv.hh"
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
    return s;
}

// ----------------------------------------------------------------------

int main(int argc, char* const argv[])
{
    try {
        argc_argv args(argc, argv,
                       {
                           {"--lab", "all", "filter by lab that submitted the sequence"},
                           // {"--virus-type", "all", "filter by subtype: H1, H3, BV, BY"},
                           {"-h", false},
                           {"--help", false},
                       });
        if (args["-h"] || args["--help"] || args.number_of_arguments() != 1)
            throw std::runtime_error("Usage: "s + args.program() + " <seqdb.json.xz>");

        seqdb::setup(std::string(args[0]), seqdb::report::yes);
        const auto& seqdb = seqdb::get(seqdb::ignore_errors::no, report_time::Yes);

        auto update = [](Info& target, const auto& entry) {
            ++target.entries;
            target.sequences += entry.number_of_seqs();
            for (const auto& seq : entry.seqs()) {
                if (!seq.hi_names().empty()) {
                    ++target.with_hi_names;
                    if (const auto date = entry.date(); !date.empty()) {
                        ++target.with_hi_names_by_month[date.substr(0, 7)];
                        ++target.with_hi_names_by_month[date.substr(0, 4) + "all"];
                    }
                }
                if (const auto date = entry.date(); !date.empty()) {
                    ++target.all_by_month[date.substr(0, 7)];
                    ++target.all_by_month[date.substr(0, 4) + "all"];
                }
                for (const auto& clade : seq.clades())
                    ++target.clades[clade];
            }
        };

        auto filter_lab = [lab = args["--lab"].str()](const auto& entry) -> bool {
            if (lab.empty() || lab == "all")
                return true;
            return entry.has_lab(lab);
        };

        std::map<std::string, Info> by_virus_type;
        Info total;
        for (const auto& entry : seqdb.entries()) {
            if (filter_lab(entry)) {
                update(total, entry);
                const auto& vt = entry.virus_type();
                update(by_virus_type[vt], entry);
                if (vt.empty())
                    std::cerr << "No virus_type for " << entry.name() << '\n';
                else if (vt == "B") {
                    const auto& lineage = entry.lineage();
                    if (lineage.empty())
                        std::cerr << "No lineage for " << entry.name() << '\n';
                    else
                        update(by_virus_type[vt + lineage], entry);
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
