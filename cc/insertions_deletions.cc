#include "insertions_deletions.hh"

using namespace seqdb;

// ----------------------------------------------------------------------

InsertionsDeletionsDetector::InsertionsDeletionsDetector(Seqdb& aSeqdb, std::string aVirusType)
{
    auto iter = aSeqdb.begin();
    iter.filter_subtype(aVirusType);
    for (; iter != aSeqdb.end(); ++iter) {
        try {
            mEntries.emplace_back(*iter);
        }
        catch (SequenceNotAligned&) {
        }
    }
      // std::cerr << mEntries.size() << " seqs for " << aVirusType << std::endl;

} // InsertionsDeletionsDetector::InsertionsDeletionsDetector

// ----------------------------------------------------------------------

void InsertionsDeletionsDetector::detect()
{

} // InsertionsDeletionsDetector::detect

// ----------------------------------------------------------------------

    // std::map<std::string, std::vector<size_t>> by_virus_type;
    // split_by_virus_type(by_virus_type);
    //   // std::cerr << by_virus_type << std::endl;
    // for (const auto& vt_indices: by_virus_type) {
    //     std::vector<std::string> aas;
    //     // std::set<std::string> last_aa;
    //     // std::set<size_t> aa_len;
    //     std::vector<std::map<char,std::vector<std::string>>> aas_per_pos;
    //     std::string common;
    //     for (size_t index: vt_indices.second) {
    //         for (const auto& seq: mEntries[index].seqs()) {
    //             if (seq.aligned()) {
    //                 const std::string aa = seq.amino_acids(true);
    //                 aas.push_back(aa);
    //                 // last_aa.emplace(1, aa.back());
    //                 // aa_len.insert(aa.size());
    //                 if (aa.size() > aas_per_pos.size())
    //                     aas_per_pos.resize(aa.size());
    //                 for (size_t pos = 0; pos < aa.size(); ++pos) {
    //                     if (aa[pos] != 'X' && aa[pos] != '-')
    //                         aas_per_pos[pos].emplace(aa[pos], std::vector<std::string>{}).first->second.push_back(aa);
    //                 }
    //             }
    //         }
    //     }
    //     std::cerr << vt_indices.first << " " << aas.size() << std::endl;
    //     if (vt_indices.first == "B") {
    //         size_t last_common_pos = 0;
    //         for (size_t pos = 0; pos < aas_per_pos.size(); ++pos) {
    //             if (aas_per_pos[pos].size() == 1)
    //                 last_common_pos = pos;
    //         }
    //         std::cerr << "last_common_pos: " << last_common_pos << std::endl;
    //         for (size_t pos = last_common_pos; pos < std::min(last_common_pos + 10, aas_per_pos.size()); ++pos) {
    //             std::cerr << pos;
    //             for (const auto& e: aas_per_pos[pos])
    //                 std::cerr << ' ' << e.first << ':' << e.second.size();
    //             std::cerr << std::endl;
    //             if (pos == 161) {
    //                 for (const auto& seq: aas_per_pos[pos]['D'])
    //                     std::cerr << seq << std::endl;
    //             }
    //         }
    //     }
    //     // std::cerr << "last AA: " << last_aa << std::endl;
    //     // std::cerr << "AA len: " << aa_len << std::endl;
    //     // if (vt_indices.first == "B")
    //     //     std::copy(aas.begin(), aas.end(), std::ostream_iterator<std::string>(std::cerr, "\n"));
    // }


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
