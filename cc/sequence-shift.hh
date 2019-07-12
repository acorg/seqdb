#pragma once

#include <stdexcept>
#include <string>
#include <limits>

// ----------------------------------------------------------------------

namespace seqdb
{
    class InvalidShift : public std::runtime_error
    {
     public:
        inline InvalidShift() : std::runtime_error("invalid shift") {}
    };

    class Shift
    {
     public:
        using ShiftT = int;

         // all entries have shift (i.e. aligned) starting 2019-07-10

        static constexpr ShiftT NotAligned = std::numeric_limits<ShiftT>::max();
        static constexpr ShiftT AlignmentFailed = std::numeric_limits<ShiftT>::max() - 1;
          // static constexpr ShiftT SequenceTooShort = std::numeric_limits<ShiftT>::max() - 2;

        Shift() : mShift(0) {} // all entries have shift (i.e. aligned) by default starting 2019-07-10
        Shift(ShiftT aShift) : mShift(aShift) {}
        Shift& operator=(ShiftT aShift) { mShift = aShift; return *this; }
        Shift& operator=(std::string::difference_type aShift) { mShift = static_cast<ShiftT>(aShift); return *this; }
        Shift& operator-=(std::string::difference_type a) { mShift -= static_cast<ShiftT>(a); return *this; }
        Shift& operator-=(size_t a) { mShift -= static_cast<ShiftT>(a); return *this; }
        bool operator==(Shift aShift) const { try { return mShift == aShift.mShift; } catch (InvalidShift&) { return false; } }
        bool operator!=(Shift aShift) const { return !operator==(aShift); }

        bool aligned() const { return mShift != NotAligned && mShift != AlignmentFailed /* && mShift != SequenceTooShort */; }
        bool alignment_failed() const { return mShift == AlignmentFailed; }

        operator ShiftT() const
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

        ShiftT raw() const { return mShift; }
        ShiftT& raw() { return mShift; }

        operator std::string() const
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

        void reset() { mShift = NotAligned; }

          // ShiftT to_json() const
          //     {
          //         try {
          //             return *this;
          //         }
          //         catch (InvalidShift&) {
          //             throw json::no_value();
          //         }
          //     }

          // void from_json(ShiftT& source)
          //     {
          //         mShift = source;
          //     }

     private:
        ShiftT mShift;

    }; // class Shift

} // namespace seqdb

// ----------------------------------------------------------------------

inline std::ostream& operator << (std::ostream& out, seqdb::Shift aShift)
{
    return out << static_cast<std::string>(aShift);
}

// ----------------------------------------------------------------------
