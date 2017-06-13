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
    if (!mEntries.empty()) {
        mMaster = mEntries.front().amino_acids;
        // std::cerr << mVirusType << ": master  " << mEntries.front().entry_seq.entry().virus_type() << ' ' << mEntries.front().entry_seq.make_name() << " " << mMaster << std::endl;
    }

} // InsertionsDeletionsDetector::InsertionsDeletionsDetector

// ----------------------------------------------------------------------

void InsertionsDeletionsDetector::detect()
{
    if (!mEntries.empty() /* && mVirusType == "B" */) {
        align_to_master();

        AAsPerPos aas_per_pos;
        aas_per_pos.collect(mEntries);
        const auto common_pos = aas_per_pos.common_pos();
          // std::cerr << mVirusType << ": last common: " << common_pos.back() << " total: " << common_pos.size() /* <<  ' '  << common_pos */ << std::endl;

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
            std::cout << mVirusType << ": " << num_with_deletions << " sequences with deletions detected, total sequences: " << mEntries.size() << std::endl;
    }

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

// ----------------------------------------------------------------------

static inline size_t number_of_common(std::string a, size_t start_a, std::string b, size_t start_b)
{
    size_t num = 0;
    for (; start_a < a.size() && start_b < b.size(); ++start_a, ++start_b) {
        if (InsertionsDeletionsDetector::Entry::common(a[start_a], b[start_b]))
            ++num;
    }
    return num;
}

static inline size_t number_of_common(std::string a, std::string b)
{
    return number_of_common(a, 0, b, 0);
}

static inline size_t number_of_common_before(std::string a, std::string b, size_t last)
{
    size_t num = 0;
    last = std::min(std::min(last, a.size()), b.size());
    for (size_t pos = 0; pos < last; ++pos) {
        if (InsertionsDeletionsDetector::Entry::common(a[pos], b[pos]))
            ++num;
    }
    return num;
}

// ----------------------------------------------------------------------

struct DeletionPos
{
    inline DeletionPos(size_t aPos, size_t aNumDeletions, size_t aNumCommon) : pos(aPos), num_deletions(aNumDeletions), num_common(aNumCommon) {}
    inline bool operator<(const DeletionPos& aNother) const { return num_common == aNother.num_common ? pos < aNother.pos : num_common > aNother.num_common; }
    size_t pos, num_deletions, num_common;
    friend inline std::ostream& operator<<(std::ostream& out, const DeletionPos& aPos) { return out << "pos:" << aPos.pos << " num_deletions:" << aPos.num_deletions << " num_common:" << aPos.num_common; }
};

using DeletionPosSet = std::vector<DeletionPos>;

static inline void update(DeletionPosSet& pos_set, std::string master, std::string to_align, size_t pos, size_t left_common)
{
    constexpr const size_t max_num_deletions = 2;
    const size_t last_pos = std::min(master.size(), to_align.size());
    if ((pos + max_num_deletions) < last_pos) {
        for (size_t num_insert = 1; num_insert <= max_num_deletions; ++num_insert) {
            if (InsertionsDeletionsDetector::Entry::common(master[pos + num_insert], to_align[pos])) {
                pos_set.emplace_back(pos, num_insert, left_common + number_of_common(master, pos + num_insert, to_align, pos));
            }
        }
    }
}

std::vector<std::pair<size_t, size_t>> InsertionsDeletionsDetector::Entry::align_to(std::string master, std::string& to_align, size_t min_common)
{
    std::vector<std::pair<size_t, size_t>> pos_number;
    if (to_align.size() > min_common) {
        size_t start = 0;
        while (start < to_align.size()) {
            const size_t current_common = number_of_common(to_align, master);
            try {
                DeletionPosSet pos_set;
                adjust_pos pos = adjust_pos::begin(to_align, master, start), pos_end = adjust_pos::end(to_align, master);
                start = to_align.size();
                for (; pos != pos_end; ++pos) {
                    update(pos_set, master, to_align, *pos, number_of_common_before(master, to_align, *pos));
                }
                if (!pos_set.empty()) {
                    const auto& del_pos = *std::min_element(pos_set.begin(), pos_set.end());
                    if (del_pos.num_common > current_common) {
                        // std::cerr << "del_pos: " << del_pos << std::endl;
                        to_align.insert(del_pos.pos, del_pos.num_deletions, '-');
                        start = del_pos.pos + del_pos.num_deletions + 1;
                        pos_number.emplace_back(del_pos.pos, del_pos.num_deletions);
                    }
                }
            }
            catch (NoAdjustPos&) {
                  // std::cerr << "NoAdjustPos: NC " << number_of_common(to_align, master) << std::endl;
                throw NotAlignedTo{};
            }
        }
    }
    // if (!pos_number.empty()) // && pos_number.front().second > 1)
    //     std::cerr << pos_number << ' '  << to_align << std::endl;
    return pos_number;

} // InsertionsDeletionsDetector::Entry::align_to

// ----------------------------------------------------------------------

void InsertionsDeletionsDetector::Entry::apply_pos_number()
{
    for (const auto& pn: pos_number) {
        entry_seq.seq().add_deletions(pn.first, pn.second);
    }

} // InsertionsDeletionsDetector::Entry::apply_pos_number

// ----------------------------------------------------------------------

void BLineageDetector::detect()
{
    auto iter = mSeqdb.begin();
    iter.filter_subtype("B");
    for (; iter != mSeqdb.end(); ++iter) {
        auto entry_seq = *iter;
        const auto& seq = entry_seq.seq();
        try {
            const char a162 = seq.amino_acid_at(162), a163 = seq.amino_acid_at(163);
            const std::string stored_lineage = entry_seq.entry().lineage();
            const std::string detected_lineage = (a162 != '-' && a163 == '-') ? "YAMAGATA" : "VICTORIA"; // if deletion in both 162 and 163, it's Vic 2016-2017 outlier
            if (stored_lineage.empty())
                entry_seq.entry().lineage(detected_lineage);
            else if (stored_lineage != detected_lineage)
                std::cerr << "WARNING: lineage conflict: " << entry_seq.make_name() << "  stored: " << stored_lineage << " detected by sequence: " << detected_lineage << std::endl << "  " << seq.amino_acids(true) << std::endl;
        }
        catch (SequenceNotAligned&) {
        }
    }

} // BLineageDetector::detect

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
