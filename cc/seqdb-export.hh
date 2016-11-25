#pragma once

#include <string>

#include "seqdb.hh"

// ----------------------------------------------------------------------

void seqdb_export(std::string aFilename, const Seqdb& aSeqdb, size_t aIndent);
// void seqdb_export_pretty(std::string aFilename, const Seqdb& aSeqdb);
void seqdb_import(std::string aFilename, Seqdb& aSeqdb);

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
