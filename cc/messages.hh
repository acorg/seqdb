#pragma once

#include <sstream>
#include "acmacs-base/string.hh"

// ----------------------------------------------------------------------

namespace seqdb
{
    class Messages
    {
     public:
        inline Messages() = default;

        inline std::ostream& warning() { return mWarnings; }

        inline operator std::string() const { return string::strip(mWarnings.str()); }
        inline operator bool() const { return !mWarnings.str().empty(); }

        inline void add(const Messages& aSource)
            {
                mWarnings << static_cast<std::string>(aSource);
            }

     private:
        std::stringstream mWarnings;

    }; // class Messages

} // namespace seqdb

inline std::ostream& operator << (std::ostream& out, const seqdb::Messages& m)
{
    return out << static_cast<std::string>(m);
}

// ----------------------------------------------------------------------
