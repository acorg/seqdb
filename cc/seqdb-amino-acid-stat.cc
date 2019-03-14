#include <iostream>
#include <sstream>
#include <algorithm>
#include <array>
#include <vector>

#include "acmacs-base/argv.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/normalize.hh"
#include "acmacs-base/range.hh"
#include "seqdb.hh"

// ----------------------------------------------------------------------

constexpr const size_t MAX_NUM_POSITIONS = 1000;
// constexpr const size_t FIRST_AA = static_cast<size_t>('A');
// constexpr const size_t NUM_AAS = static_cast<size_t>('Z') + 1 - FIRST_AA;

using position_entry_t = std::vector<std::pair<char, unsigned>>;
using amino_acids_data_t = std::vector<position_entry_t>;

struct Options;

static std::pair<std::string_view, std::string_view> extract_date_range(std::string_view date_range);
static amino_acids_data_t collect(const seqdb::Seqdb& seqdb, std::string_view date_range, const Options& opt);
static void report(const amino_acids_data_t& amino_acids_data, const Options& opt);
static void report(const amino_acids_data_t& amino_acids_data_1, const amino_acids_data_t& amino_acids_data_2, const Options& opt);

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str>   flu{*this, "flu", desc{"filter by virus type: H1, H3, B"}};
    option<str>   lineage{*this, "lineage", desc{"filter by lineage"}};
    option<str>   clade{*this, "clade", desc{"filter by clade"}};
    option<str>   date_range_1{*this, "date-range", desc{"YYYY-MM-DD:YYYY-MM-DD, parts optional, semicolon required"}};
    option<str>   date_range_2{*this, "date-range-2", desc{"YYYY-MM-DD:YYYY-MM-DD, parts optional, semicolon required"}};
    option<bool>  report_all_pos{*this, "report-all-pos", desc{"report all positions"}};
};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        // seqdb::setup(opt.seqdb_file, seqdb::report::yes);
        const auto& seqdb = seqdb::get(seqdb::ignore_errors::no, report_time::yes);
        const amino_acids_data_t amino_acids_data = collect(seqdb, opt.date_range_1, opt);
        if (opt.date_range_2.has_value())
            report(amino_acids_data, collect(seqdb, opt.date_range_2, opt), opt);
        else
            report(amino_acids_data, opt);
    }
    catch (std::exception& err) {
        std::cerr << "ERROR: " << err.what() << '\n';
        return 1;
    }
}

// ----------------------------------------------------------------------

amino_acids_data_t collect(const seqdb::Seqdb& seqdb, std::string_view date_range, const Options& opt)
{
    size_t num_sequences = 0;
    const auto [start, end] = extract_date_range(date_range);
    auto iter = seqdb.begin();
    iter.filter_subtype(acmacs::normalize_virus_type(*opt.flu))
            .filter_lineage(acmacs::normalize_lineage(*opt.lineage))
            .filter_clade(string::upper(*opt.clade))
            .filter_date_range(start, end)
            ;

    amino_acids_data_t amino_acids_data(MAX_NUM_POSITIONS);
    size_t max_pos = 0;
    for (; iter != seqdb.end(); ++iter) {
        const auto seq = (*iter).seq().amino_acids(true);
        for (size_t pos = 0; pos < seq.size(); ++pos) {
            auto& at = amino_acids_data[pos];
            if (const auto found = std::find_if(std::begin(at), std::end(at), [aa = seq[pos]](const auto& entry) { return entry.first == aa; }); found != std::end(at))
                ++found->second;
            else
                at.emplace_back(seq[pos], 1);
        }
        max_pos = std::max(max_pos, seq.size());
        ++num_sequences;
    }
    std::cout << "sequences: " << num_sequences << " max-pos: " << max_pos << '\n';
    amino_acids_data.resize(max_pos);

    std::for_each(std::begin(amino_acids_data), std::end(amino_acids_data),
                  [](auto& entry) { std::sort(std::begin(entry), std::end(entry), [](const auto& e1, const auto& e2) { return e1.second > e2.second; }); });

    return amino_acids_data;

} // collect

// ----------------------------------------------------------------------

static std::pair<std::string_view, std::string_view> extract_date_range(std::string_view date_range)
{
    if (!date_range.empty()) {
        const auto parts = acmacs::string::split(date_range, ":");
        if (parts.size() != 2)
            throw std::runtime_error(string::concat("Invalid date range: ", date_range));
        return {parts[0], parts[1]};
    }
    return {std::string_view{}, std::string_view{}};

} // extract_date_range

// ----------------------------------------------------------------------

void report(const amino_acids_data_t& amino_acids_data, const Options& opt)
{
    auto positions_sorted = acmacs::filled_with_indexes(amino_acids_data.size());
    std::sort(std::begin(positions_sorted), std::end(positions_sorted),
              [&amino_acids_data](size_t p1, size_t p2) { return amino_acids_data[p1].front().second < amino_acids_data[p2].front().second; });

    for (auto pos : positions_sorted) {
        if (opt.report_all_pos || amino_acids_data[pos].size() > 1) {
            std::cout << std::setw(3) << (pos + 1);
            for (const auto& [aa, num] : amino_acids_data[pos])
                std::cout << ' ' << aa << ':' << num;
            std::cout << '\n';
        }
    }

} // report

// ----------------------------------------------------------------------

void report(const amino_acids_data_t& amino_acids_data_1, const amino_acids_data_t& amino_acids_data_2, const Options& opt)
{
    const auto write_aas = [](const auto& data) {
        std::stringstream s1;
        for (const auto& [aa, num] : data)
            s1 << ' ' << aa << ':' << num;
        return s1.str();
    };

    for (size_t pos = 0; pos < std::min(amino_acids_data_1.size(), amino_acids_data_2.size()); ++pos) {
        const auto& aad1 = amino_acids_data_1[pos];
        const auto& aad2 = amino_acids_data_2[pos];
        const char aa1 = aad1.front().first, aa2 = aad2.front().first;
        if (opt.report_all_pos || aad1.size() > 1 || aad2.size() > 1 || aa1 != aa2) {
            std::cout << std::setw(3) << std::right << (pos + 1) << "  " << aa1;
            if (aa1 != aa2)
                std::cout << aa2;
            else
                std::cout << ' ';
            std::cout << "  " << std::setw(30) << std::left << write_aas(aad1) << write_aas(aad2) << '\n';
        }
    }

} // report

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
