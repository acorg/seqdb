#pragma once

#include <stdexcept>
#include <string>
#include <limits>

// ----------------------------------------------------------------------

class InvalidShift : public std::runtime_error
{
 public:
    inline InvalidShift() : std::runtime_error("invalid shift") {}
};

class Shift
{
 public:
    typedef int ShiftT;

    static constexpr ShiftT NotAligned = std::numeric_limits<ShiftT>::max();
    static constexpr ShiftT AlignmentFailed = std::numeric_limits<ShiftT>::max() - 1;
    // static constexpr ShiftT SequenceTooShort = std::numeric_limits<ShiftT>::max() - 2;

    inline Shift(ShiftT aShift = NotAligned) : mShift(aShift) {}
    inline Shift& operator=(ShiftT aShift) { mShift = aShift; return *this; }
    inline Shift& operator=(std::string::difference_type aShift) { mShift = static_cast<ShiftT>(aShift); return *this; }
    inline Shift& operator-=(std::string::difference_type a) { mShift -= static_cast<ShiftT>(a); return *this; }
    inline Shift& operator-=(size_t a) { mShift -= static_cast<ShiftT>(a); return *this; }
    inline bool operator==(Shift aShift) const { try { return mShift == aShift.mShift; } catch (InvalidShift&) { return false; } }
    inline bool operator!=(Shift aShift) const { return !operator==(aShift); }

    inline bool aligned() const { return mShift != NotAligned && mShift != AlignmentFailed /* && mShift != SequenceTooShort */; }
    inline bool alignment_failed() const { return mShift == AlignmentFailed; }

    inline operator ShiftT() const
        {
            switch (mShift) {
              case NotAligned:
              case AlignmentFailed:
              // case SequenceTooShort:
                  throw InvalidShift();
              default:
                  break;
            }
            return static_cast<int>(mShift);
        }

    inline operator std::string() const
        {
            switch (mShift) {
              case NotAligned:
                  return "NotAligned";
              case AlignmentFailed:
                  return "AlignmentFailed";
              default:
                  return std::string("Aligned: ") + std::to_string(mShift);
            }
        }

    inline void reset() { mShift = NotAligned; }

    // inline ShiftT to_json() const
    //     {
    //         try {
    //             return *this;
    //         }
    //         catch (InvalidShift&) {
    //             throw json::no_value();
    //         }
    //     }

    // inline void from_json(ShiftT& source)
    //     {
    //         mShift = source;
    //     }

 private:
    ShiftT mShift;

}; // class Shift

inline std::ostream& operator << (std::ostream& out, Shift aShift) { return out << static_cast<std::string>(aShift); }

// ----------------------------------------------------------------------
