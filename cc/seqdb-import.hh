#pragma once

#include <string>

// ----------------------------------------------------------------------

namespace seqdb
{
    class Seqdb;
    void seqdb_import(std::string aFilename, Seqdb& aSeqdb);
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
