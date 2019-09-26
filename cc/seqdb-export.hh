#pragma once

#include <string>

// ----------------------------------------------------------------------

namespace seqdb
{
    class Seqdb;
    void seqdb_export(std::string_view aFilename, const Seqdb& aSeqdb, size_t aIndent);
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
