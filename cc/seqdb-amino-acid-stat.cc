#include <iostream>
#include "acmacs-base/argv.hh"
#include "acmacs-base/string.hh"
#include "acmacs-base/normalize.hh"
#include "seqdb.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str>   flu{*this, "flu", desc{"filter by virus type: H1, H3, B"}};
    option<str>   lineage{*this, "lineage", desc{"filter by lineage"}};
    option<str>   clade{*this, "clade", desc{"filter by clade"}};
    option<str>   start_date{*this, "start-date", desc{"isolated on or after date: YYYY-MM-DD"}};
    option<str>   end_date{*this, "end-date", desc{"isolated before date: YYYY-MM-DD"}};
};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        // seqdb::setup(opt.seqdb_file, seqdb::report::yes);
        const auto& seqdb = seqdb::get(seqdb::ignore_errors::no, report_time::yes);
        size_t num_sequences = 0;
        auto iter = seqdb.begin();
        iter.filter_subtype(acmacs::normalize_virus_type(*opt.flu))
                .filter_lineage(acmacs::normalize_lineage(*opt.lineage))
                .filter_clade(string::upper(*opt.clade))
                .filter_date_range(opt.start_date, opt.end_date)
                ;
        for (; iter != seqdb.end(); ++iter) {
            ++num_sequences;
        }
        std::cout << "sequences: " << num_sequences << '\n';
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
