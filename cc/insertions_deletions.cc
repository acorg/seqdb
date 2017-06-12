#include "insertions_deletions.hh"

using namespace seqdb;

// ----------------------------------------------------------------------

class AAsPerPosEntry : public std::string
{
 public:
    inline AAsPerPosEntry() {}
    inline void add(char aAA) { if (aAA != 'X' && aAA != '-' && find(aAA) == npos) append(1, aAA); }
    inline bool common() const { return size() == 1; }
};

class AAsPerPos : public std::vector<AAsPerPosEntry>
{
 public:
    inline AAsPerPos() {}
    inline void adjust_size(size_t new_size) { if (size() < new_size) resize(new_size); }
    inline void update(std::string aa) { adjust_size(aa.size()); for (size_t pos = 0; pos < aa.size(); ++pos) operator[](pos).add(aa[pos]); }
    inline void collect(const InsertionsDeletionsDetector::Entries& entries) { std::for_each(entries.begin(), entries.end(), [this](const auto& entry) { update(entry.amino_acids); }); }
    template <typename Func> inline void collect_if(const InsertionsDeletionsDetector::Entries& entries, Func aFunc) { std::for_each(entries.begin(), entries.end(), [this,&aFunc](const auto& entry) { if (aFunc(entry)) update(entry.amino_acids); }); }
    inline size_t last_common() const { size_t last = static_cast<size_t>(-1); for (size_t pos = 0; pos < size(); ++pos) { if (operator[](pos).common()) last = pos; } return last; }
};

// ----------------------------------------------------------------------

InsertionsDeletionsDetector::InsertionsDeletionsDetector(Seqdb& aSeqdb, std::string aVirusType)
    : mVirusType(aVirusType)
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
    AAsPerPos aas_per_pos;
    aas_per_pos.collect(mEntries);
    const size_t last_common = aas_per_pos.last_common();
    std::cerr << "last_common: " << last_common << std::endl;
    if (last_common < 500) {
        std::vector<size_t> last_common_for_aa;
        for (char aa: aas_per_pos[last_common + 1]) {
            AAsPerPos aas_per_pos_for_aa;
            aas_per_pos_for_aa.collect_if(mEntries, [&](const Entry& entry) -> bool { return entry.amino_acids[last_common + 1] == aa; });
            last_common_for_aa.push_back(aas_per_pos_for_aa.last_common());
        }
        std::cerr << "last_common " << (last_common + 1) << ' ' << last_common_for_aa << std::endl;
    }
    else {
        std::cout << "No insertions/deletions for " << mVirusType << std::endl;
    }

      // for (size_t pos = last_common + 1; pos < std::min(last_common + 10, aas_per_pos.size()); ++pos)
      //     std::cerr << pos << ' ' << aas_per_pos[pos] << std::endl;

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
