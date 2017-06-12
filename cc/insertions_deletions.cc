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

class NoAdjustPos : public std::exception { public: using std::exception::exception; };
class NotAlignedTo : public std::exception { public: using std::exception::exception; };

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
    mMaster = mEntries.front().amino_acids;
    // std::cerr << "Master " << mMaster.size() << "                 " << mMaster << std::endl;
      // std::cerr << mEntries.size() << " seqs for " << aVirusType << std::endl;

} // InsertionsDeletionsDetector::InsertionsDeletionsDetector

// ----------------------------------------------------------------------

void InsertionsDeletionsDetector::detect()
{
    align_to_master();

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

void InsertionsDeletionsDetector::align_to_master()
{
    const size_t min_common = static_cast<size_t>(std::floor(mMaster.size() * 0.8));
      // std::cerr << mVirusType << ": min_common: " << min_common << std::endl;
    bool restart = true;
    while (restart) {
        restart = false;
        for (auto& entry: mEntries) {
            try {
                entry.pos_number = entry.align_to(mMaster, entry.amino_acids, min_common);
            }
            catch (NotAlignedTo&) {
                try {
                      // perhaps master has deletions
                    entry.align_to(entry.amino_acids, mMaster, min_common);
                      // yes, change master to current and restart
                    mMaster = entry.amino_acids;
                    // std::cerr << mVirusType << ": master changed to " << mMaster.size() << std::endl;
                    restart = true;
                    revert();
                    break;
                }
                catch (NotAlignedTo&) {
                    std::cerr << mVirusType << ": cannot find deletions in " << entry.entry_seq.make_name() << std::endl;
                }
            }
        }
    }

} // InsertionsDeletionsDetector::align_to_master

// ----------------------------------------------------------------------

std::vector<std::pair<size_t, size_t>> InsertionsDeletionsDetector::Entry::align_to(std::string master, std::string to_align, size_t min_common)
{
    std::vector<std::pair<size_t, size_t>> pos_number;
    if (to_align.size() > min_common) {
        while (true) {
            const size_t current_common = number_of_common(to_align, master);
            if (current_common > min_common)
                break;
            try {
                const size_t adjust_pos = find_adjust_pos(master, to_align);
                size_t max_common = current_common;
                size_t best_num_insert = 0, num_insert;
                for (num_insert = 1; num_insert <= 5; ++num_insert) {
                    to_align.insert(adjust_pos, 1, '-');
                    const size_t common = number_of_common(to_align, master);
                    if (common > max_common) {
                        max_common = common;
                        best_num_insert = num_insert;
                    }
                }
                to_align.erase(adjust_pos, num_insert - best_num_insert - 1);
                if (best_num_insert) {
                    pos_number.emplace_back(adjust_pos, best_num_insert);
                }
                else {
                      // std::cerr << "No best_num_insert: current_common:" << current_common << " num_insert:" << num_insert << std::endl << master << std::endl << to_align << std::endl;
                    throw NotAlignedTo{};
                }
            }
            catch (NoAdjustPos&) {
                  // std::cerr << "NoAdjustPos: NC " << number_of_common(to_align, master) << std::endl;
                throw NotAlignedTo{};
            }
        }
    }
    return pos_number;

} // InsertionsDeletionsDetector::Entry::align_to

// ----------------------------------------------------------------------

size_t InsertionsDeletionsDetector::Entry::find_adjust_pos(std::string master, std::string to_align)
{
    constexpr const size_t min_gap = 5;

    const size_t last_pos = std::min(to_align.size(), master.size());
    for (size_t pos = 0, previous_common = 0; pos < last_pos; ++pos) {
        if (common(to_align[pos], master[pos])) {
            if ((pos - previous_common) > min_gap) {
                return previous_common + 1;
            }
            previous_common = pos;
        }
    }
    std::cerr << "NoAdjustPos common:" << number_of_common(to_align, master) << std::endl << master << std::endl << to_align << std::endl;
    throw NoAdjustPos{};

} // InsertionsDeletionsDetector::Entry::find_adjust_pos

// ----------------------------------------------------------------------

size_t InsertionsDeletionsDetector::Entry::number_of_common(std::string a, std::string b)
{
    const size_t last_pos = std::min(a.size(), b.size());
    size_t number_of_common = 0;
    for (size_t pos = 0; pos < last_pos; ++pos)
        if (common(a[pos], b[pos]))
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
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
