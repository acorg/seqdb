#include <stack>

#include "rapidjson/reader.h"

#include "seqdb-export.hh"
#include "json-writer.hh"
#include "acmacs-base/read-file.hh"

// ----------------------------------------------------------------------

template <typename RW> inline JsonWriterT<RW>& operator <<(JsonWriterT<RW>& writer, const Seqdb& seqdb)
{
    writer.StartObject();
    writer.Key("  version");
    writer.String("sequence-database-v2");
    writer.Key("data");
    return writer << StartArray
                  << EndArray
                  << EndObject;

    // << StartObject
    //               << JsonKey::TableId << chart.table_id()
    //               << if_not_empty(JsonKey::Virus, chart.chart_info().virus())
    //               << if_not_empty(JsonKey::VirusType, chart.chart_info().virus_type())
    //               << if_not_empty(JsonKey::Assay, chart.chart_info().assay())
    //               << if_not_empty(JsonKey::Date, chart.chart_info().date())
    //               << if_not_empty(JsonKey::Lab, chart.chart_info().lab())
    //               << if_not_empty(JsonKey::Rbc, chart.chart_info().rbc())
    //               << if_not_empty(JsonKey::Name, chart.chart_info().name())
    //           // conflicts with JsonKey::Sera << if_not_empty(JsonKey::VirusSubset, chart.chart_info().subset())
    //               << JsonKey::Antigens << chart.antigens()
    //               << JsonKey::Sera << chart.sera()
    //               << JsonKey::Titers << chart.titers()
    //               << EndObject;
}

// ----------------------------------------------------------------------

void seqdb_export(std::string aFilename, const Seqdb& aSeqdb, size_t aIndent)
{
    acmacs_base::write_file(aFilename, json(aSeqdb, "seqdb", aIndent));
}

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
