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
    inline void collect(const InsertionsDeletionsDetector::Entries& entries) { std::for_each(entries.begin(), entries.end(), [this](const auto& entry) { this->update(entry.amino_acids); }); }
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
    std::cerr << mVirusType << ": master " << mMaster << std::endl;
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
            // std::cerr << '+' << std::endl;
            try {
                entry.pos_number = entry.align_to(mMaster, entry.amino_acids, min_common);
            }
            catch (NotAlignedTo&) {
                try {
                      // perhaps master has deletions
                    entry.align_to(entry.amino_acids, mMaster, min_common);
                      // yes, change master to current and restart
                    mMaster = entry.amino_acids;
                    std::cerr << mVirusType << ": master changed to " << mMaster << std::endl;
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

class adjust_pos
{
 public:
    static inline adjust_pos begin(std::string to_align, std::string master, size_t pos = 0) { return {to_align, master, pos}; }
    static inline adjust_pos end(std::string to_align, std::string master) { return {to_align, master}; }

    inline bool operator==(const adjust_pos& an) const { return /* mToAlign == an.mToAlign && mMaster == an.mMaster && */ mPos == an.mPos; }
    inline bool operator!=(const adjust_pos& an) const { return ! operator==(an); }

    inline size_t operator*() const { return mPos; }
      // inline const size_t *operator->() const { return &mPos; }
    inline adjust_pos& operator++() { ++mPos; find(); return *this; }

 private:
    inline adjust_pos(std::string to_align, std::string master, size_t pos)
        : mToAlign(to_align), mMaster(master), mLastPos(std::min(to_align.size(), master.size())), mPos(pos) { find(); }
    inline adjust_pos(std::string to_align, std::string master)
        : mToAlign(to_align), mMaster(master), mLastPos(std::min(to_align.size(), master.size())), mPos(mLastPos) { find(); }
    std::string mToAlign, mMaster;
    size_t mLastPos, mPos;

    inline void find() { while (mPos < mLastPos && (InsertionsDeletionsDetector::Entry::common(mToAlign[mPos], mMaster[mPos]) || mToAlign[mPos] == '-' || mMaster[mPos] == '-')) ++mPos; }

}; // class adjust_pos

std::vector<std::pair<size_t, size_t>> InsertionsDeletionsDetector::Entry::align_to(std::string master, std::string& to_align, size_t min_common)
{
    std::vector<std::pair<size_t, size_t>> pos_number;
    if (to_align.size() > min_common) {
        size_t start = 0;
        while (start < to_align.size()) {
            const size_t current_common = number_of_common(to_align, master);
            if (current_common > min_common)
                break;
            try {
                adjust_pos pos = adjust_pos::begin(to_align, master, start), pos_end = adjust_pos::end(to_align, master);
                start = to_align.size();
                for (; pos != pos_end; ++pos) {
                      // std::cerr << *pos << std::endl;
                    const size_t best_num_insert = check_adjust_pos(master, to_align, *pos, current_common);
                    if (best_num_insert) {
                        pos_number.emplace_back(*pos, best_num_insert);
                        start = *pos + 1;
                        if (best_num_insert == 2)
                            std::cerr << *pos << ' '  << to_align << std::endl;
                        break;
                    }
                }

                // for (size_t adjust_pos = 0; adjust_pos < to_align.size(); ++adjust_pos) {
                //     adjust_pos = find_adjust_pos(master, to_align, adjust_pos);
                //     const size_t best_num_insert = check_adjust_pos(master, to_align, adjust_pos, current_common);
                //     if (best_num_insert) {
                //         pos_number.emplace_back(adjust_pos, best_num_insert);
                //         if (best_num_insert == 2)
                //             std::cerr << adjust_pos << ' '  << to_align << std::endl;
                //         break;
                //     }
                // }
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

size_t InsertionsDeletionsDetector::Entry::check_adjust_pos(std::string master, std::string& to_align, size_t adjust_pos, size_t current_common)
{
    constexpr const size_t max_num_deletions = 2;

    size_t max_common = current_common;
    size_t best_num_insert = 0, num_insert;
    for (num_insert = 1; num_insert <= max_num_deletions; ++num_insert) {
        to_align.insert(adjust_pos, 1, '-');
        const size_t common = number_of_common(to_align, master);
        if (common > max_common) {
            std::cerr << "adjust_pos:" << adjust_pos << " num_insert:" << num_insert << " common:" << common << " current_common:" << current_common << std::endl;
            max_common = common;
            best_num_insert = num_insert;
        }
    }
    to_align.erase(adjust_pos, num_insert - best_num_insert - 1);
    return best_num_insert;

} // InsertionsDeletionsDetector::Entry::check_adjust_pos

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
