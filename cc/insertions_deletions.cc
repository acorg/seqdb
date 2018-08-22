#include <iomanip>
#include <array>

#include "acmacs-base/debug.hh"
#include "acmacs-base/counter.hh"
#include "seqdb/insertions_deletions.hh"
#include "seqdb/amino-acids.hh"

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

// class NoAdjustPos : public std::exception { public: using std::exception::exception; };
class SwitchMaster : public std::exception { public: using std::exception::exception; };

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
    choose_master();

} // InsertionsDeletionsDetector::InsertionsDeletionsDetector

// ----------------------------------------------------------------------

void InsertionsDeletionsDetector::choose_master()
{
    if (!mEntries.empty()) {
        size_t master_number_aa = 0;
        if (mVirusType == "A(H1N1)")
            master_number_aa = NORMAL_SEQUENCE_AA_LENGTH_H1;
        else if (mVirusType == "A(H3N2)")
            master_number_aa = NORMAL_SEQUENCE_AA_LENGTH_H3;
        else if (mVirusType == "B")
            master_number_aa = NORMAL_SEQUENCE_AA_LENGTH_B;
        else if (mEntries.size() > 10)
            std::cerr << "WARNING: unknown normal sequence size for " << mVirusType << '\n';

        std::string master_name;
        if (master_number_aa > 0) {
            for (auto& entry : mEntries) {
                if (entry.amino_acids.size() == master_number_aa) {
                    mMaster = entry.amino_acids;
                    master_name = entry.entry_seq.make_name();
                    break;
                }
            }
        }
        if (mMaster.empty() || master_name.empty()) {
            mMaster = mEntries.front().amino_acids;
            master_name = mEntries.front().entry_seq.make_name();
            master_switching_allowed_ = true;
        }
        else {
            master_switching_allowed_ = false;
        }
        std::cout << "INFO: " << mVirusType << ": master: " << master_name << "\n                 " << mMaster << std::endl;
    }
    // else
    //     std::cerr << "ERROR: InsertionsDeletionsDetector::choose_master: no entries\n";

} // InsertionsDeletionsDetector::choose_master

// ----------------------------------------------------------------------

void InsertionsDeletionsDetector::detect()
{
    if (!mEntries.empty() /* && mVirusType == "B" */) {
        std::cout << "INFO: detecting insertions/deletions in " << mVirusType << '\n';

        // acmacs::Counter counter(mEntries.begin(), mEntries.end(), [](const auto& entry) { return entry.amino_acids.size(); });
        // std::cerr << "DEBUG: seq-lengths: " << counter << '\n';

        align_to_master();

        AAsPerPos aas_per_pos;
        aas_per_pos.collect(mEntries);
        const auto common_pos = aas_per_pos.common_pos();
          // std::cerr << mVirusType << ": last common: " << common_pos.back() << " total: " << common_pos.size() /* <<  ' '  << common_pos */ << std::endl;

        size_t num_with_deletions = 0;
        for (auto& entry: mEntries) {
            if (!entry.pos_number.empty()) {
                try {
                    entry.apply_pos_number();
                    ++num_with_deletions;
                    // if (entry.pos_number.size() > 1)
                    //     std::cerr << entry.entry_seq.make_name() << ' ' << entry.pos_number << ' ' << entry.amino_acids << '\n';
                }
                catch (InvalidShift&) {
                }
            }
        }
        if (num_with_deletions)
            std::cout << "INFO: " << mVirusType << ": " << num_with_deletions << " sequences with deletions detected, total sequences: " << mEntries.size() << std::endl;
    }

} // InsertionsDeletionsDetector::detect

// ----------------------------------------------------------------------

void InsertionsDeletionsDetector::align_to_master()
{
    // const size_t min_common = static_cast<size_t>(std::floor(mMaster.size() * 0.8));
    // std::cerr << mVirusType << ": min_common: " << min_common << std::endl;
    bool restart = true;
    while (restart) {
        restart = false;
        for (auto& entry : mEntries) {
            try {
                entry.pos_number = entry.align_to(mMaster, entry.amino_acids, entry.entry_seq);
            }
            catch (SwitchMaster&) {
                if (master_switching_allowed_ && entry.amino_acids.size() > std::lround(mMaster.size() * 0.9)) { // do not switch master, if new sequence is too short
                    try {
                        // perhaps master has deletions
                        entry.align_to(entry.amino_acids, mMaster, entry.entry_seq);
                        // yes, change master to current and restart
                        mMaster = entry.amino_acids;
                        std::cout << "INFO: " << mVirusType << ": master changed to " << entry.entry_seq.make_name() << '\n' << mMaster << '\n';
                        restart = true;
                        revert();
                        break;
                    }
                    catch (SwitchMaster&) {
                        // std::cerr << mVirusType << ": cannot find deletions in " << entry.entry_seq.make_name() << std::endl;
                    }
                }
                // else {
                //     std::cout << "INFO: " << mVirusType << ": master NOT changed to " << entry.entry_seq.make_name() << ": too short: " << entry.amino_acids.size() << " amino-acids\n";
                // }
            }
        }
    }

} // InsertionsDeletionsDetector::align_to_master

// ----------------------------------------------------------------------

static inline size_t number_of_common(const std::string a, size_t start_a, const std::string b, size_t start_b)
{
    size_t num = 0;
    for (; start_a < a.size() && start_b < b.size(); ++start_a, ++start_b) {
        if (InsertionsDeletionsDetector::Entry::common(a[start_a], b[start_b]))
            ++num;
    }
    return num;
}

static inline size_t number_of_common(const std::string a, const std::string b)
{
    return number_of_common(a, 0, b, 0);
}

static inline size_t number_of_common_before(const std::string a, const std::string b, size_t last)
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

class adjust_pos
{
 public:
    static inline adjust_pos begin(std::string to_align, std::string master, size_t pos = 0) { return {to_align, master, pos}; }
    static inline adjust_pos end(std::string to_align, std::string master) { return {to_align, master}; }

    bool operator==(const adjust_pos& an) const { return /* mToAlign == an.mToAlign && mMaster == an.mMaster && */ mPos == an.mPos; }
    bool operator!=(const adjust_pos& an) const { return ! operator==(an); }

    size_t operator*() const { return mPos; }
      // const size_t *operator->() const { return &mPos; }
    adjust_pos& operator++() { ++mPos; find(); return *this; }

 private:
    adjust_pos(std::string to_align, std::string master, size_t pos)
        : mToAlign(to_align), mMaster(master), mLastPos(std::min(to_align.size(), master.size())), mPos(pos)
        {
            // acmacs::debug dbg(true);
            // dbg << "adjust_pos1 pos:" << mPos << " last_pos:" << mLastPos << '\n';
            find();
            // dbg << "adjust_pos2 pos:" << mPos << " last_pos:" << mLastPos << '\n';
        }
    adjust_pos(std::string to_align, std::string master)
        : mToAlign(to_align), mMaster(master), mLastPos(std::min(to_align.size(), master.size())), mPos(mLastPos) {}
    std::string mToAlign, mMaster;
    size_t mLastPos, mPos;

    void find()
    {
        while (mPos < mLastPos && (InsertionsDeletionsDetector::Entry::common(mToAlign[mPos], mMaster[mPos]) || mToAlign[mPos] == '-' || mMaster[mPos] == '-'))
            ++mPos;
    }

}; // class adjust_pos

// ----------------------------------------------------------------------

struct DeletionPos
{
    DeletionPos(size_t aPos, size_t aNumDeletions, size_t aNumCommon) : pos(aPos), num_deletions(aNumDeletions), num_common(aNumCommon) {}
    bool operator<(const DeletionPos& aNother) const { return num_common == aNother.num_common ? pos < aNother.pos : num_common > aNother.num_common; }
    size_t pos, num_deletions, num_common;
    void fix(size_t aPos, size_t aNumDeletions, size_t aNumCommon) { pos = aPos; num_deletions = aNumDeletions; num_common = aNumCommon; }
    friend inline std::ostream& operator<<(std::ostream& out, const DeletionPos& aPos) { return out << "pos:" << aPos.pos << " num_deletions:" << aPos.num_deletions << " num_common:" << aPos.num_common; }
};

using DeletionPosSet = std::vector<DeletionPos>;

static inline void update(DeletionPosSet& pos_set, std::string master, std::string to_align, size_t pos, size_t common_before)
{
    constexpr const size_t max_num_deletions = 5;
    const size_t last_pos = std::min(master.size(), to_align.size());
    // acmacs::debug dbg(true);
    // dbg << "update1 " << pos << " common_before:" << common_before << '\n';
    if ((pos + max_num_deletions) < last_pos) {
        for (size_t num_insert = 1; num_insert <= max_num_deletions; ++num_insert) {
            if (InsertionsDeletionsDetector::Entry::common(master[pos + num_insert], to_align[pos])) {
                pos_set.emplace_back(pos, num_insert, common_before + number_of_common(master, pos + num_insert, to_align, pos));
                // dbg << "update2 " << pos << " deletions:" << num_insert << '\n';
            }
        }
    }
}

std::vector<std::pair<size_t, size_t>> InsertionsDeletionsDetector::Entry::align_to(const std::string master, std::string& to_align, const SeqdbEntrySeq& entry_seq)
{
    // acmacs::debug dbg(false);

    constexpr const bool yamagata_163_hack = true;
    constexpr const bool victoria_tripledel2017_hack = true;
    // dbg << '\n' << entry_seq.make_name() << '\n' << "align: " << to_align << '\n' << "maste: " << master << '\n';
    std::vector<std::pair<size_t, size_t>> pos_number;
    size_t start = 0;
    size_t best_common = number_of_common(master, to_align);
    const std::string to_align_orig = to_align;
    while (start < to_align.size()) {
        const size_t current_common = number_of_common(master, to_align);
        // try {
        DeletionPosSet pos_set;
        adjust_pos pos = adjust_pos::begin(to_align, master, start), pos_end = adjust_pos::end(to_align, master);
        start = to_align.size();
        for (; pos != pos_end; ++pos) {
            update(pos_set, master, to_align, *pos, number_of_common_before(master, to_align, *pos));
        }
        // dbg << "pos_set: " << pos_set << '\n';
        if (!pos_set.empty()) {
            auto& del_pos = *std::min_element(pos_set.begin(), pos_set.end());
            // dbg << "del_pos: " << del_pos << '\n';
            if (del_pos.num_common > current_common) {
                if (yamagata_163_hack && entry_seq.entry().virus_type() == "B" && del_pos.num_deletions == 1 && del_pos.pos > (163 - 1) && del_pos.pos <= (166 - 1)) {
                      // std::cout << "INFO: yamagata_163_hack applied for " << entry_seq.make_name() << '\n';
                    // yamagata deletion must be at 163
                    // David Burke 2017-08-17: deletions ( and insertions) of amino acids usually occur in regions of the protein structure where it changes direction ( loops ).
                    // In the case of HA, this is after VPK and before NKTAT/YKNAT.
                    to_align.insert(163 - 1, 1, '-');
                    del_pos.fix(163 - 1, 1, number_of_common(master, to_align)); // -1 because we count from zero here
                }
                else if (victoria_tripledel2017_hack && entry_seq.entry().virus_type() == "B" && del_pos.num_deletions == 3 && del_pos.pos == (164 - 1)) {
                    // The triple deletion is 162, 163 and 164 (pos 1 based). this is the convention that has been chosen (Sarah 2018-08-16 08:31)
                    to_align.insert(162 - 1, del_pos.num_deletions, '-');
                    del_pos.fix(162 - 1, del_pos.num_deletions, number_of_common(master, to_align)); // -1 because we count from zero here
                }
                else {
                    to_align.insert(del_pos.pos, del_pos.num_deletions, '-');
                }
                // dbg << "del_pos: " << del_pos << '\n';
                start = del_pos.pos + del_pos.num_deletions + 1;
                pos_number.emplace_back(del_pos.pos, del_pos.num_deletions);
                if (best_common < del_pos.num_common)
                    best_common = del_pos.num_common;
            }
        }
        // }
        // catch (NoAdjustPos&) {
        //       // std::cerr << "NoAdjustPos: NC " << number_of_common(to_align, master) << '\n';
        //     throw SwitchMaster{};
        // }
    }
    // if (!pos_number.empty()) // && pos_number.front().second > 1)
    //     dbg << pos_number << " best-common:" << best_common << '\n' << "align: " << to_align << '\n' << "maste: " << master << '\n';
    if (best_common < static_cast<size_t>(master.size() * 0.7)) {
        // dbg << "Too bad matching (common:" << best_common << " threshold:" << (master.size() * 0.7) << "), should we switch master?" << '\n';
        to_align = to_align_orig;
        throw SwitchMaster{};
    }
    // dbg << to_align << '\n';
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
            // if (std::initializer_list<size_t> pos_to_check = {160, 161, 164, 165, 166, 167, 68, 169, 170}; std::any_of(std::begin(pos_to_check), std::end(pos_to_check), [&seq](size_t pos) { return seq.amino_acid_at(pos) == '-'; })) {
            //       // throw std::runtime_error("Invalid deletion in B sequence for " + iter.make_name());
            //     std::cerr << "ERROR: Invalid deletion in B sequence for " << iter.make_name() << std::endl;
            // }
            const std::string stored_lineage = entry_seq.entry().lineage();
            const std::string detected_lineage = (seq.amino_acid_at(162) != '-' && (seq.amino_acid_at(163) == '-' || seq.amino_acid_at(164) == '-' || seq.amino_acid_at(165) == '-' || seq.amino_acid_at(166) == '-')) ? "YAMAGATA" : "VICTORIA"; // if deletion in both 162 and 163, it's Vic 2016-2017 outlier
            // std::cerr << detected_lineage << ' ' << entry_seq.make_name() << '\n' << seq.amino_acids(true) << '\n';
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
