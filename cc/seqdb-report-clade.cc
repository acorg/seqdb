#include <iostream>

#include "acmacs-base/argv.hh"
#include "acmacs-base/normalize.hh"
#include "seqdb.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str>  db_dir{*this, "db-dir"};
    option<str>  lab{*this, "lab"};
    option<str>  flu{*this, "flu", desc{"A(H1N1), A(H3N2), H1, H3 B"}};
    option<bool> sort_by_date{*this, "sort-by-date"};

    argument<str> clade{*this, arg_name{"clade"}, mandatory};
};

struct Entry
{
    std::string name;
    std::string date;
    std::string lab;
    std::string seq_id;
};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);
        seqdb::setup_dbs(opt.db_dir, seqdb::report::no);
        std::string flu{acmacs::normalize_virus_type(opt.flu)};
        if (*opt.clade == "all") {
            for (const auto entry_seq : seqdb::get()) {
                if ((opt.lab->empty() || entry_seq.seq().has_lab(*opt.lab)) && (flu.empty() || entry_seq.entry().virus_type() == flu))
                    std::cout << std::setw(60) << std::left << entry_seq.make_name() << '\t' << entry_seq.seq().clades() << '\t' << entry_seq.entry().virus_type() << '\t'
                              << entry_seq.entry().lineage() << '\t' << entry_seq.entry().dates() << '\t' << entry_seq.seq().lab() << '\n';
            }
        }
        else {
            const auto& seqdb = seqdb::get();
            std::vector<Entry> seqs;

            for (const auto entry_seq : seqdb) {
                if (entry_seq.seq().has_clade(*opt.clade) && (opt.lab->empty() || entry_seq.seq().has_lab(*opt.lab)) && (flu.empty() || entry_seq.entry().virus_type() == flu))
                    seqs.push_back({entry_seq.make_name(), std::string{entry_seq.entry().date()}, std::string{entry_seq.seq().lab()}, entry_seq.seq_id(seqdb::SeqdbEntrySeq::encoded_t::yes)});
            }
            if (opt.sort_by_date)
                std::sort(std::begin(seqs), std::end(seqs), [](const auto& e1, const auto& e2) { return e1.date < e2.date; });
            for (const auto& entry : seqs) {
                std::cout << std::setw(60) << std::left << entry.name << '\t' << entry.date << '\t' << entry.lab << '\t' << entry.seq_id << '\n';
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
