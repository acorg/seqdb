#include <string>
#include <fstream>

#include "acmacs-base/argc-argv.hh"
#include "acmacs-base/string.hh"
#include "acmacs-base/stream.hh"
#include "acmacs-base/counter.hh"
#include "acmacs-base/read-file.hh"
#include "seqdb/seqdb.hh"
#include "seqdb/hamming-distance.hh"

// ----------------------------------------------------------------------

using entry_seq_iter_t = std::vector<seqdb::SeqdbEntrySeq>::const_iterator;

static std::pair<std::string, std::string> parse_flu(std::string flu);
static seqdb::SeqdbEntrySeq find_base_seq(std::string virus_type, std::string lineage, std::string base_seq_regex);
static std::vector<seqdb::SeqdbEntrySeq> collect(std::string virus_type, std::string lineage, const seqdb::SeqdbEntrySeq& base_seq);
static std::pair<size_t, std::vector<seqdb::SeqdbEntrySeq>> pick(const std::vector<seqdb::SeqdbEntrySeq>& sequences, const seqdb::SeqdbEntrySeq& base_seq, size_t number_to_pick, size_t hamming_distance_threshold);
static size_t common_length(entry_seq_iter_t first, entry_seq_iter_t last);
static void destroy_failed_sequence_end(std::string seq_id, std::string& aa, std::string& nucs, std::string base_seq_aa, std::string base_seq_nucs);

// ----------------------------------------------------------------------

constexpr const char* sUsage = " [options] --flu <H1, H3, BV, BY> --base-seq <\"VICTORIA/830/2013 MDCK2\"> --output-nuc <source.fasta.xz>\n";

int main(int argc, char* const argv[])
{
    try {
        using namespace std::string_literals;
        argc_argv args(argc, argv,
                       {
                           {"--db-dir", ""},
                           {"--flu", "", "(mandatory) H1, H3, BV, BY"},
                           {"--recent", 4000, "number of the recent sequences to export (excluding base sequence)"},
                           {"--hamming-distance-threshold", 160, "export only sequences having hamming distance to the base sequence less than threshold"},
                           {"--base-seq", "", "(mandatory) sequence to become root of the tree (regex)"},
                           {"--output-aa", "", "fasta file with amino-acids"},
                           {"--output-nucs", "", "(mandatory) fasta file with nucleotides"},
                           {"--output-names", "", "text file with exported sequence names"},

                           {"-v", false},
                           {"--verbose", false},
                           {"-h", false},
                           {"--help", false},
                       });
        if (args["-h"] || args["--help"] || args.number_of_arguments() || !args["--flu"] || !args["--base-seq"] || !args["--output-nucs"]) {
            throw std::runtime_error("Usage: "s + args.program() + sUsage + args.usage_options());
        }
        const bool verbose = args["-v"] || args["--verbose"];
        const auto [virus_type, lineage] = parse_flu(std::string(args["--flu"]));
        // const size_t recent = args["--recent"];
        // const size_t hamming_distance_threshold = args["--hamming-distance-threshold"];
        seqdb::setup_dbs(std::string(args["--db-dir"]), verbose ? seqdb::report::yes : seqdb::report::no);

        const auto base_seq = find_base_seq(virus_type, lineage, std::string(args["--base-seq"]));
        const auto all_sequences_sorted_by_date = collect(virus_type, lineage, base_seq);
        auto [aa_common_length, to_export] = pick(all_sequences_sorted_by_date, base_seq, args["--recent"], args["--hamming-distance-threshold"]);
        // add base_seq
        to_export.insert(to_export.begin(), base_seq);

        // write
        const std::string base_seq_nucs = base_seq.seq().nucleotides(true, 0, aa_common_length * 3);
        const std::string base_seq_aa = base_seq.seq().amino_acids(true, 0, aa_common_length);
        acmacs::file::ofstream output_nucs(args["--output-nucs"]);
        acmacs::file::ofstream output_aa(args["--output-aa"] ? args["--output-aa"].str() : "/dev/null"s);
        for (const auto& seq : to_export) {
            std::string nucs = seq.seq().nucleotides(true, 0, aa_common_length * 3);
            std::string aa = seq.seq().amino_acids(true, 0, aa_common_length);
            destroy_failed_sequence_end(seq.seq_id(seqdb::SeqdbEntrySeq::encoded_t::yes), aa, nucs, base_seq_aa, base_seq_nucs);
            output_nucs.stream() << '>' << seq.seq_id(seqdb::SeqdbEntrySeq::encoded_t::yes) << '\n' << nucs << '\n';
            output_aa.stream() << '>' << seq.seq_id(seqdb::SeqdbEntrySeq::encoded_t::yes) << '\n' << aa << '\n';
        }
        if (args["--output-names"]) {
            acmacs::file::ofstream output(args["--output-names"]);
            for (const auto& seq : to_export)
                output.stream() << seq.seq_id(seqdb::SeqdbEntrySeq::encoded_t::yes) << '\n';
        }
        return 0;
    }
    catch (std::exception& err) {
        std::cerr << "ERROR: " << err.what() << '\n';
        return 1;
    }
}

// ----------------------------------------------------------------------

std::pair<std::string, std::string> parse_flu(std::string flu)
{
    flu = string::upper(flu);
    if (flu == "H1")
        return {"A(H1N1)", ""};
    else if (flu == "H3")
        return {"A(H3N2)", ""};
    else if (flu[0] == 'B' && flu[1] == 'V')
        return {"B", "VICTORIA"};
    else if (flu[0] == 'B' && flu[1] == 'Y')
        return {"B", "YAMAGATA"};
    else
        throw std::runtime_error("Unrecognized --flu argument: \"" + flu + "\" (H1, H3, BV, BY expected)");

} // parse_flu

// ----------------------------------------------------------------------

seqdb::SeqdbEntrySeq find_base_seq(std::string virus_type, std::string lineage, std::string base_seq_regex)
{
    auto begin = seqdb::get().begin();
    begin.filter_subtype(virus_type).filter_lineage(lineage).filter_aligned(true).filter_gene("HA").filter_name_regex(base_seq_regex);
    std::vector<seqdb::SeqdbEntrySeq> sequences;
    std::copy(begin, seqdb::get().end(), std::back_inserter(sequences));
    if (sequences.empty())
        throw std::runtime_error("No base seqs found");
    if (sequences.size() > 1) {
        for (const auto& es : sequences)
            std::cerr << "ERROR: base_seq: " << es.make_name() << '\n';
        throw std::runtime_error("Too many base seqs found");
    }
    std::cerr << "INFO: base_seq: " << sequences[0].make_name() << '\n';
    return sequences[0];

} // find_base_seq

// ----------------------------------------------------------------------

std::vector<seqdb::SeqdbEntrySeq> collect(std::string virus_type, std::string lineage, const seqdb::SeqdbEntrySeq& base_seq)
{
    auto begin = seqdb::get().begin();
    begin.filter_subtype(virus_type).filter_lineage(lineage).filter_aligned(true).filter_gene("HA");
    std::vector<seqdb::SeqdbEntrySeq> sequences;
    std::copy_if(begin, seqdb::get().end(), std::back_inserter(sequences), [&base_seq](const auto& es) { return es != base_seq; });
    std::cerr << "INFO: " << sequences.size() << " total sequences found\n";

    std::sort(sequences.begin(), sequences.end(), [](const auto& es1, const auto& es2) { return es1.seq_id(seqdb::SeqdbEntrySeq::encoded_t::no) < es2.seq_id(seqdb::SeqdbEntrySeq::encoded_t::no); });
    const auto uniq_end = std::unique(sequences.begin(), sequences.end(), [](const auto& es1, const auto& es2) { return es1.seq_id(seqdb::SeqdbEntrySeq::encoded_t::no) == es2.seq_id(seqdb::SeqdbEntrySeq::encoded_t::no); });
    if (uniq_end != sequences.end()) {
        std::cerr << "INFO: seq_id duplicates removed: " << (sequences.end() - uniq_end) << '\n';
        sequences.erase(uniq_end, sequences.end());
    }

    // most recent first
    std::sort(sequences.begin(), sequences.end(), [](const auto& es1, const auto& es2) { return es1.entry().date() > es2.entry().date(); });
    return sequences;

} // collect

// ----------------------------------------------------------------------

std::pair<size_t, std::vector<seqdb::SeqdbEntrySeq>> pick(const std::vector<seqdb::SeqdbEntrySeq>& sequences, const seqdb::SeqdbEntrySeq& base_seq, size_t number_to_pick, size_t hamming_distance_threshold)
{
    const size_t aa_common_length = common_length(sequences.begin(), sequences.size() > number_to_pick ? sequences.begin() + static_cast<std::vector<seqdb::SeqdbEntrySeq>::difference_type>(number_to_pick) : sequences.end());
      // std::cerr << "DEBUG: common_length: " << aa_common_length << '\n';
    const std::string base_seq_nucs = base_seq.seq().nucleotides(true, 0, aa_common_length * 3);
    const std::string base_seq_nucs_last_21 = base_seq_nucs.substr(aa_common_length * 3 - 21);
    std::vector<seqdb::SeqdbEntrySeq> result;
    // acmacs::Counter<size_t> counter;
    for (auto seqp = sequences.begin(); seqp != sequences.end() && result.size() < number_to_pick; ++seqp) {
        const std::string nucs = seqp->seq().nucleotides(true, 0, aa_common_length * 3);
        const auto hamming_distance_with_base = acmacs::seqdb::hamming_distance(base_seq_nucs, nucs);
        if (hamming_distance_with_base < hamming_distance_threshold) {
            result.push_back(*seqp);
        }

        // const std::string last_21 = nucs.substr(aa_common_length * 3 - 21);
        // const auto hamming_distance_last_21 = string::hamming_distance(base_seq_nucs_last_21, last_21);
        // counter.add(hamming_distance_last_21);
        // if (hamming_distance_last_21 > 6)
        //     std::cerr << "DEBUG: " << seqp->seq_id(seqdb::SeqdbEntrySeq::encoded_t::no) << " HD:" << hamming_distance_last_21 << '\n';
    }
    // std::cerr << "DEBUG: humming_at_last_21: " << counter << '\n';

    return {aa_common_length, result};

} // pick

// ----------------------------------------------------------------------

static size_t common_length(entry_seq_iter_t first, entry_seq_iter_t last)
{
    return acmacs::Counter(first, last, [](const auto& es) { return es.seq().amino_acids(true).size(); }).max().first;

} // common_length

// ----------------------------------------------------------------------

void destroy_failed_sequence_end(std::string seq_id, std::string& aa, std::string& nucs, std::string /*base_seq_aa*/, std::string base_seq_nucs)
{
      // found in B/Vic before SSM TC1 2018-08-14, messy sequence end

    constexpr size_t last_aa_pos_to_check = 13;
    const auto nucs_size = nucs.size();
    const auto aa_size = aa.size();

    std::vector<size_t> dist(last_aa_pos_to_check + 1, 0);
    for (size_t last_pos = last_aa_pos_to_check; last_pos > 0; --last_pos) {
        dist[last_pos] = acmacs::seqdb::hamming_distance(base_seq_nucs.substr(nucs_size - last_pos * 3), nucs.substr(nucs_size - last_pos * 3));
    }
    for (size_t last_pos = last_aa_pos_to_check; last_pos > 1; --last_pos) {
        if (dist[last_pos] > dist[last_pos - 1] && dist[last_pos] >= last_pos) {
              // too much mess found at the end, perhaps sequencing failure, replace end with '-'
            std::cerr << "INFO: " << seq_id << " last " << last_pos * 3 << " nucs replaced with -, sequencing failure at the end suspected\n"
                      << "        base seq end: " << base_seq_nucs.substr(nucs_size - last_aa_pos_to_check * 3) << '\n'
                      << "             seq end: " << nucs.substr(nucs_size - last_aa_pos_to_check * 3) << '\n';
            nucs.resize(nucs_size - last_pos * 3);
            nucs.resize(nucs_size, '-');
            aa.resize(aa_size - last_pos);
            aa.resize(aa_size, 'X');
            std::cerr << "   resulting seq end: " << nucs.substr(nucs_size - last_aa_pos_to_check * 3) << '\n';
            std::cerr << "DEBUG: " << dist << '\n';
            break;
        }
    }

} // destroy_failed_sequence_end

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
