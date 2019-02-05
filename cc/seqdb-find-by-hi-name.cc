#include "acmacs-base/argv.hh"
#include "seqdb/seqdb.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str> db_dir{*this, "db-dir"};

    argument<str_array> names{*this, arg_name{"name"}, mandatory};
};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        seqdb::setup_dbs(opt.db_dir, seqdb::report::no);
        const auto& seqdb = seqdb::get();

        for (const auto& name : *opt.names) {
            if (const auto* entry_seq = seqdb.find_hi_name(name); entry_seq) {
                std::cout << entry_seq->seq_id(seqdb::SeqdbEntrySeq::encoded_t::yes) << '\n';
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
