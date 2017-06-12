#include <iomanip>

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
    inline std::vector<size_t> common_pos() const { std::vector<size_t> common; for (size_t pos = 0; pos < size(); ++pos) { if (operator[](pos).common()) { common.push_back(pos); } } return common; }
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
            if (mEntries.back().amino_acids.size() > mLongest.size())
                mLongest = mEntries.back().amino_acids;
        }
        catch (SequenceNotAligned&) {
        }
    }
    // std::cerr << "Longest " << mLongest.size() << "                 " << mLongest << std::endl;
      // std::cerr << mEntries.size() << " seqs for " << aVirusType << std::endl;

} // InsertionsDeletionsDetector::InsertionsDeletionsDetector

// ----------------------------------------------------------------------

void InsertionsDeletionsDetector::detect()
{
    align_to_longest();

    AAsPerPos aas_per_pos;
    aas_per_pos.collect(mEntries);
    const auto common_pos = aas_per_pos.common_pos();
    std::cerr << mVirusType << ": last common: " << common_pos.back() << " total: " << common_pos.size() /* <<  ' '  << common_pos */ << std::endl;

    size_t num_with_deletions = 0;
    for (auto& entry: mEntries) {
        if (!entry.pos_number.empty()) {
            entry.apply_pos_number();
            ++num_with_deletions;
            // if (entry.pos_number.front().second > 1)
            //     std::cerr << entry.entry_seq.make_name() << std::endl << entry.pos_number << ' ' << entry.amino_acids << std::endl;
        }
    }
    if (num_with_deletions)
        std::cout << mVirusType << ": " << num_with_deletions << " sequences with deletions detected" << std::endl;

} // InsertionsDeletionsDetector::detect

// ----------------------------------------------------------------------

void InsertionsDeletionsDetector::align_to_longest()
{
    const size_t min_common = static_cast<size_t>(std::floor(mLongest.size() * 0.8));
      // std::cerr << mVirusType << ": min_common: " << min_common << std::endl;
    for (auto& entry: mEntries) {
        entry.align_to(mLongest, min_common);
    }

} // InsertionsDeletionsDetector::align_to_longest

// ----------------------------------------------------------------------

class NoAdjustPos : public std::exception { public: using std::exception::exception; };

void InsertionsDeletionsDetector::Entry::align_to(std::string master, size_t min_common)
{
    while (amino_acids.size() > min_common && number_of_common(master) < min_common) {
        try {
            const size_t adjust_pos = find_adjust_pos(master);
            size_t max_common = 0;
            size_t best_num_insert = 0, num_insert;
            for (num_insert = 1; num_insert <= 5; ++num_insert) {
                amino_acids.insert(adjust_pos, 1, '-');
                const size_t common = number_of_common(master);
                if (common > max_common) {
                    max_common = common;
                    best_num_insert = num_insert;
                }
            }
            amino_acids.erase(adjust_pos, num_insert - best_num_insert - 1);
            pos_number.emplace_back(adjust_pos, best_num_insert);
        }
        catch (NoAdjustPos&) {
              // std::cerr << "NoAdjustPos: NC " << number_of_common(master) << std::endl;
            break;
        }
    }

} // InsertionsDeletionsDetector::Entry::align_to

// ----------------------------------------------------------------------

size_t InsertionsDeletionsDetector::Entry::find_adjust_pos(std::string master) const
{
    constexpr const size_t min_gap = 5;

    const size_t last_pos = std::min(amino_acids.size(), master.size());
    for (size_t pos = 0, previous_common = 0; pos < last_pos; ++pos) {
        if (common(amino_acids[pos], master[pos])) {
            if ((pos - previous_common) > min_gap) {
                return previous_common + 1;
            }
            previous_common = pos;
        }
    }
      // std::cerr << "NoAdjustPos common:" << number_of_common(master) << std::endl << master << std::endl << amino_acids << std::endl;
    throw NoAdjustPos{};

} // InsertionsDeletionsDetector::Entry::find_adjust_pos

// ----------------------------------------------------------------------

size_t InsertionsDeletionsDetector::Entry::number_of_common(std::string master) const
{
    const size_t last_pos = std::min(amino_acids.size(), master.size());
    size_t number_of_common = 0;
    for (size_t pos = 0; pos < last_pos; ++pos)
        if (common(amino_acids[pos], master[pos]))
            ++number_of_common;
    return number_of_common;

} // InsertionsDeletionsDetector::Entry::number_of_common

// ----------------------------------------------------------------------

void InsertionsDeletionsDetector::Entry::apply_pos_number()
{
    for (const auto& pn: pos_number) {
        entry_seq.seq().add_deletions(pn.first, pn.second);
    }

} // InsertionsDeletionsDetector::Entry::apply_pos_number

// ----------------------------------------------------------------------


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
