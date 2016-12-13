#include "seqdb-export.hh"
#include "json-keys.hh"
#include "acmacs-base/json-writer.hh"
#include "acmacs-base/json-reader.hh"

using namespace seqdb;

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
                  << if_not_empty(SeqdbJsonKey::LabIds, seq.lab_ids_raw())
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

void seqdb::seqdb_export(std::string aFilename, const Seqdb& aSeqdb, size_t aIndent)
{
    export_to_json(aSeqdb, SEQDB_JSON_DUMP_VERSION, aFilename, aIndent);
}

// ----------------------------------------------------------------------

using HandlerBase = json_reader::HandlerBase<Seqdb>;
using StringListHandler = json_reader::StringListHandler<Seqdb>;
using MapListHandler = json_reader::MapListHandler<Seqdb>;

// ----------------------------------------------------------------------

class SeqHandler : public HandlerBase
{
 public:
    inline SeqHandler(Seqdb& aSeqdb, std::vector<SeqdbSeq>& aSeqs)
        : HandlerBase(aSeqdb), mSeqs(aSeqs), mStarted(false), mKey(SeqdbJsonKey::Unknown) {}

    inline virtual HandlerBase* StartArray()
        {
            if (mStarted)
                throw json_reader::Failure();
            mStarted = true;
            return nullptr;
        }

    inline virtual HandlerBase* StartObject()
        {
            if (!mStarted)
                throw json_reader::Failure();
            mSeqs.emplace_back();
            return nullptr;
        }

    inline virtual HandlerBase* EndObject()
        {
            return nullptr;
        }

    inline virtual HandlerBase* Key(const char* str, rapidjson::SizeType length)
        {
            if (!mStarted)
                throw json_reader::Failure();
            if (length != 1)
                throw json_reader::Failure();
            mKey = static_cast<SeqdbJsonKey>(*str);
            HandlerBase* result = nullptr;
#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wswitch-enum"
#endif
            switch (mKey) {
              case SeqdbJsonKey::AminoAcids:
              case SeqdbJsonKey::Nucleotides:
              case SeqdbJsonKey::Gene:
              case SeqdbJsonKey::AminoAcidShift:
              case SeqdbJsonKey::NucleotideShift:
                  break;
              case SeqdbJsonKey::Clades:
                  result = new StringListHandler(mTarget, mSeqs.back().clades());
                  break;
              case SeqdbJsonKey::HiNames:
                  result = new StringListHandler(mTarget, mSeqs.back().hi_names());
                  break;
              case SeqdbJsonKey::LabIds:
                  result = new MapListHandler(mTarget, mSeqs.back().lab_ids_raw());
                  break;
              case SeqdbJsonKey::Passages:
                  result = new StringListHandler(mTarget, mSeqs.back().passages());
                  break;
              case SeqdbJsonKey::Reassortant:
                  result = new StringListHandler(mTarget, mSeqs.back().reassortant());
                  break;
              default:
                  break;
            }
#pragma GCC diagnostic pop
            return result;
        }

    inline virtual HandlerBase* Int(int i)
        {
#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wswitch-enum"
#endif
            switch (mKey) {
              case SeqdbJsonKey::AminoAcidShift:
                  mSeqs.back().amino_acids_shift_raw() = i;
                  break;
              case SeqdbJsonKey::NucleotideShift:
                  mSeqs.back().nucleotides_shift_raw() = i;
                  break;
              default:
                  throw json_reader::Failure();
            }
#pragma GCC diagnostic pop
            return nullptr;
        }

    inline virtual HandlerBase* Uint(unsigned u)
        {
            return Int(static_cast<int>(u));
        }

    inline virtual HandlerBase* String(const char* str, rapidjson::SizeType length)
        {
#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wswitch-enum"
#endif
            switch (mKey) {
              case SeqdbJsonKey::AminoAcids:
                  mSeqs.back().amino_acids().assign(str, length);
                  break;
              case SeqdbJsonKey::Nucleotides:
                  mSeqs.back().nucleotides().assign(str, length);
                  break;
              case SeqdbJsonKey::Gene:
                  mSeqs.back().gene().assign(str, length);
                  break;
              default:
                  throw json_reader::Failure();
            }
#pragma GCC diagnostic pop
            return nullptr;
        }

 private:
    std::vector<SeqdbSeq>& mSeqs;
    bool mStarted;
    SeqdbJsonKey mKey;

}; // class SeqHandler

// ----------------------------------------------------------------------

class DataHandler : public HandlerBase
{
 public:
    inline DataHandler(Seqdb& aSeqdb, std::vector<SeqdbEntry>& aData)
        : HandlerBase(aSeqdb), mData(aData), mStarted(false), mKey(SeqdbJsonKey::Unknown) {}

    inline virtual HandlerBase* StartArray()
        {
            if (mStarted) {
                throw json_reader::Failure();
            }
            else {
                mStarted = true;
            }
            return nullptr;
        }

    inline virtual HandlerBase* StartObject()
        {
            if (!mStarted)
                throw json_reader::Failure();
            mData.emplace_back();
            return nullptr;
        }

    inline virtual HandlerBase* EndObject()
        {
            return nullptr;
        }

    inline virtual HandlerBase* Key(const char* str, rapidjson::SizeType length)
        {
            if (!mStarted)
                throw json_reader::Failure();
            if (length != 1)
                throw json_reader::Failure();
            mKey = static_cast<SeqdbJsonKey>(*str);
            HandlerBase* result = nullptr;
#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wswitch-enum"
#endif
            switch (mKey) {
              case SeqdbJsonKey::Dates:
                  result = new StringListHandler(mTarget, mData.back().dates());
                  break;
              case SeqdbJsonKey::SequenceSet:
                  result = new SeqHandler(mTarget, mData.back().seqs());
                  break;
              default:
                  break;
            }
#pragma GCC diagnostic pop
            return result;
        }

    inline virtual HandlerBase* String(const char* str, rapidjson::SizeType length)
        {
#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wswitch-enum"
#endif
            switch (mKey) {
              case SeqdbJsonKey::Name:
                  mData.back().name().assign(str, length);
                  break;
              case SeqdbJsonKey::Continent:
                  mData.back().continent().assign(str, length);
                  break;
              case SeqdbJsonKey::Country:
                  mData.back().country().assign(str, length);
                  break;
              case SeqdbJsonKey::Lineage:
                  mData.back().lineage().assign(str, length);
                  break;
              case SeqdbJsonKey::VirusType:
                  mData.back().virus_type().assign(str, length);
                  break;
              default:
                  throw json_reader::Failure();
            }
#pragma GCC diagnostic pop
            return nullptr;
        }

 private:
    std::vector<SeqdbEntry>& mData;
    bool mStarted;
    SeqdbJsonKey mKey;

}; // class DataHandler

// ----------------------------------------------------------------------

class SeqdbRootHandler : public HandlerBase
{
 private:
    enum class Keys { Unknown, Version, Data };

 public:
    inline SeqdbRootHandler(Seqdb& aSeqdb) : HandlerBase(aSeqdb), mKey(Keys::Unknown) {}

    inline virtual HandlerBase* Key(const char* str, rapidjson::SizeType length)
        {
            HandlerBase* result = nullptr;
            const std::string found_key(str, length);
            if (found_key == "  version")
                mKey = Keys::Version;
            else if (found_key == "data")
                result = new DataHandler(mTarget, mTarget.entries());
            else
                result = HandlerBase::Key(str, length);
            return result;
        }

    inline virtual HandlerBase* String(const char* str, rapidjson::SizeType length)
        {
            HandlerBase* result = nullptr;
            switch (mKey) {
              case Keys::Version:
                  if (strncmp(str, SEQDB_JSON_DUMP_VERSION, std::min(length, static_cast<rapidjson::SizeType>(strlen(SEQDB_JSON_DUMP_VERSION))))) {
                      std::cerr << "Unsupported version: \"" << std::string(str, length) << '"' << std::endl;
                      throw json_reader::Failure();
                  }
                  break;
              case Keys::Data:
              case Keys::Unknown:
                  result = HandlerBase::String(str, length);
                  break;
            }
            return result;
        }

 private:
    Keys mKey;

};

// ----------------------------------------------------------------------

void seqdb::seqdb_import(std::string aFilename, Seqdb& aSeqdb)
{
    json_reader::read_from_file<Seqdb, SeqdbRootHandler>(aFilename, aSeqdb);

} // seqdb_import

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
