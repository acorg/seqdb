#include <stack>

#include "rapidjson/reader.h"

#include "seqdb-export.hh"
#include "acmacs-base/read-file.hh"

// ----------------------------------------------------------------------

void seqdb_export(std::string aFilename, const Seqdb& aSeqdb)
{
    JsonWriter writer("seqdb");
    writer.StartObject();
    writer.Key("  version");
    writer.String("sequence-database-v2");
    // writer << JsonKey::Antigens << aSeqdb.antigens() << JsonKey::Sera << aSeqdb.sera() << JsonKey::Tables << aSeqdb.charts() << EndObject;
    acmacs_base::write_file(aFilename, writer);
}

// ----------------------------------------------------------------------

void seqdb_export_pretty(std::string aFilename, const Seqdb& aSeqdb)
{
    JsonPrettyWriter writer("seqdb");

    writer.StartObject();
    writer.Key("  version");
    writer.String("sequence-database-v2");
    // writer << JsonKey::Antigens << aSeqdb.antigens() << JsonKey::Sera << aSeqdb.sera() << JsonKey::Tables << aSeqdb.charts() << EndObject;

    acmacs_base::write_file(aFilename, writer);

} // seqdb_export_pretty

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------

void seqdb_import(std::string buffer, Seqdb& aSeqdb)
{
    if (buffer == "-")
        buffer = acmacs_base::read_stdin();
    else
        buffer = acmacs_base::read_file(buffer);
    if (buffer[0] == '{') { // && buffer.find("\"  version\": \"sequence-database-v2\"") != std::string::npos) {
        // SeqdbReaderEventHandler handler{aSeqdb};
        // rapidjson::Reader reader;
        // rapidjson::StringStream ss(buffer.c_str());
        // reader.Parse(ss, handler);
        // if (reader.HasParseError())
        //     throw Error("cannot import seqdb: data parsing failed at state " + std::to_string(static_cast<unsigned>(handler.state.top())) + " at pos " + std::to_string(reader.GetErrorOffset()) + ": " +  GetParseError_En(reader.GetParseErrorCode()) + "\n" + buffer.substr(reader.GetErrorOffset(), 50));
        // if (!handler.in_init_state())
        //     throw Error("internal: not in init state on parsing completion");
    }
    else
        throw std::runtime_error("cannot import seqdb: unrecognized source format");

} // seqdb_import

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
