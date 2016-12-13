#pragma once

#include <string>

#include "seqdb/seqdb.hh"

// ----------------------------------------------------------------------

namespace seqdb
{
    void seqdb_export(std::string aFilename, const Seqdb& aSeqdb, size_t aIndent);
    void seqdb_import(std::string aFilename, Seqdb& aSeqdb);
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
