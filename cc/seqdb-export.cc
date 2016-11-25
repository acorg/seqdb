#include <stack>

#include "rapidjson/reader.h"

#include "seqdb-export.hh"
#include "json-writer.hh"

// ----------------------------------------------------------------------

static constexpr const char* SEQDB_JSON_DUMP_VERSION = "sequence-database-v2";

// ----------------------------------------------------------------------

class if_aligned
{
 public:
    inline if_aligned(SeqdbJsonKey key, Shift shift) : mKey(key), mShift(shift) {}

    template <typename RW> friend inline JsonWriterT<RW>& operator <<(JsonWriterT<RW>& writer, const if_aligned& data)
        {
            if (data.mShift.aligned())
                writer << data.mKey << data.mShift.raw();
            return writer;
        }

 private:
    SeqdbJsonKey mKey;
    Shift mShift;
};

// ----------------------------------------------------------------------

template <typename RW> inline JsonWriterT<RW>& operator <<(JsonWriterT<RW>& writer, const SeqdbSeq& seq)
{
    return writer << StartObject
                  << if_not_empty(SeqdbJsonKey::Passages, seq.passages())
                  << if_not_empty(SeqdbJsonKey::Nucleotides, seq.nucleotides(false))
                  << if_not_empty(SeqdbJsonKey::AminoAcids, seq.amino_acids(false))
                  << if_aligned(SeqdbJsonKey::NucleotideShift, seq.nucleotides_shift())
                  << if_aligned(SeqdbJsonKey::AminoAcidShift, seq.amino_acids_shift())
                  << if_not_empty(SeqdbJsonKey::LabIds, seq.lab_ids())
                  << if_not_empty(SeqdbJsonKey::Gene, seq.gene())
                  << if_not_empty(SeqdbJsonKey::HiNames, seq.hi_names())
                  << if_not_empty(SeqdbJsonKey::Reassortant, seq.reassortant())
                  << if_not_empty(SeqdbJsonKey::Clades, seq.clades())
                  << EndObject;
}

// ----------------------------------------------------------------------

template <typename RW> inline JsonWriterT<RW>& operator <<(JsonWriterT<RW>& writer, const SeqdbEntry& entry)
{
    return writer << StartObject
                  << if_not_empty(SeqdbJsonKey::Name, entry.name())
                  << if_not_empty(SeqdbJsonKey::Continent, entry.continent())
                  << if_not_empty(SeqdbJsonKey::Country, entry.country())
                  << if_not_empty(SeqdbJsonKey::Dates, entry.dates())
                  << if_not_empty(SeqdbJsonKey::Lineage, entry.lineage())
                  << if_not_empty(SeqdbJsonKey::VirusType, entry.virus_type())
                  << SeqdbJsonKey::SequenceSet << entry.seqs()
                  << EndObject;
}

// ----------------------------------------------------------------------

template <typename RW> inline JsonWriterT<RW>& operator <<(JsonWriterT<RW>& writer, const Seqdb& seqdb)
{
    return writer << StartObject
                  << JsonObjectKey("  version") << SEQDB_JSON_DUMP_VERSION
                  << JsonObjectKey("data") << seqdb.entries()
                  << EndObject;
}

// ----------------------------------------------------------------------

void seqdb_export(std::string aFilename, const Seqdb& aSeqdb, size_t aIndent)
{
    export_to_json(aSeqdb, "seqdb", aFilename, aIndent);
}

// ----------------------------------------------------------------------

void seqdb_import(std::string buffer, Seqdb& aSeqdb)
{
    if (buffer == "-")
        buffer = acmacs_base::read_stdin();
    else
        buffer = acmacs_base::read_file(buffer);
    if (buffer[0] == '{') { // && buffer.find("\"  version\": " + SEQDB_JSON_DUMP_VERSION) != std::string::npos) {
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
