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

class Error : public std::runtime_error { public: using std::runtime_error::runtime_error; };
class Failure : public std::exception { public: using std::exception::exception; };
class Pop : public std::exception { public: using std::exception::exception; };

class HandlerBase
{
 public:
    inline HandlerBase(Seqdb& aSeqdb) : mSeqdb(aSeqdb), mIgnore(false) {}
    virtual ~HandlerBase();

    inline virtual HandlerBase* StartObject() { std::cerr << "HandlerBase StartObject " << typeid(*this).name() << std::endl; throw Failure(); }
    inline virtual HandlerBase* EndObject() { throw Pop(); }
    inline virtual HandlerBase* StartArray() { throw Failure(); }
    inline virtual HandlerBase* EndArray() { throw Pop(); }
    inline virtual HandlerBase* Double(double d) { std::cerr << "Double: " << d << std::endl; throw Failure(); }
    inline virtual HandlerBase* Int(int i) { std::cerr << "Int: " << i << std::endl; throw Failure(); }
    inline virtual HandlerBase* Uint(unsigned u) { std::cerr << "Uint: " << u << std::endl; throw Failure(); }

    inline virtual HandlerBase* Key(const char* str, rapidjson::SizeType length)
        {
            if ((length == 1 && *str == '_') || (length > 0 && *str == '?')) {
                mIgnore = true;
            }
            else {
                std::cerr << "Key: \"" << std::string(str, length) << '"' << std::endl;
                throw Failure();
            }
            return nullptr;
        }

    inline virtual HandlerBase* String(const char* str, rapidjson::SizeType length)
        {
            if (mIgnore) {
                mIgnore = false;
            }
            else {
                std::cerr << "String: \"" << std::string(str, length) << '"' << std::endl;
                throw Failure();
            }
            return nullptr;
        }

 protected:
    Seqdb& mSeqdb;
    bool mIgnore;
};

HandlerBase::~HandlerBase()
{
}

// ----------------------------------------------------------------------

class DataHandler : public HandlerBase
{
 public:
    inline DataHandler(Seqdb& aSeqdb, std::vector<SeqdbEntry>& aData)
        : HandlerBase(aSeqdb), mData(aData), mStarted(false), mKey(SeqdbJsonKey::Unknown) {}

    inline virtual HandlerBase* StartArray()
        {
            if (mStarted) {
                throw Failure();
            }
            else {
                mStarted = true;
            }
            return nullptr;
        }

    inline virtual HandlerBase* StartObject()
        {
            if (!mStarted)
                throw Failure();
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
                throw Failure();
            if (length != 1)
                throw Failure();
            mKey = static_cast<SeqdbJsonKey>(*str);
            return nullptr;
        }

    inline virtual HandlerBase* Int(int i)
        {
            switch (mKey) {
              case SeqdbJsonKey::AminoAcidShift:
                  break;
              case SeqdbJsonKey::NucleotideShift:
                  break;
              default:
                  throw Failure();
            }
            return nullptr;
        }

    inline virtual HandlerBase* String(const char* str, rapidjson::SizeType length)
        {
            switch (mKey) {
              case SeqdbJsonKey::Name:
                  mData.back().name().assign(str, length);
                  break;
              default:
                  throw Failure();
            }
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
            Keys new_key = Keys::Unknown;
            const std::string found_key(str, length);
            if (found_key == "  version")
                mKey = Keys::Version;
            else if (found_key == "data")
                result = new DataHandler(mSeqdb, mSeqdb.entries());
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
                      throw Failure();
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

class DocRootHandler : public HandlerBase
{
 public:
    inline DocRootHandler(Seqdb& aSeqdb) : HandlerBase(aSeqdb) {}

    inline virtual HandlerBase* StartObject() { return new SeqdbRootHandler(mSeqdb); }
};

// ----------------------------------------------------------------------

class SeqdbReaderEventHandler : public rapidjson::BaseReaderHandler<rapidjson::UTF8<>, SeqdbReaderEventHandler>
{
 private:
    template <typename... Args> inline bool handler(HandlerBase* (HandlerBase::*aHandler)(Args... args), Args... args)
        {
            try {
                auto new_handler = ((*mHandler.top()).*aHandler)(args...);
                if (new_handler)
                    mHandler.emplace(new_handler);
            }
            catch (Pop&) {
                if (mHandler.empty())
                    return false;
                mHandler.pop();
            }
            catch (Failure&) {
                return false;
            }
            return true;
        }

 public:
    inline SeqdbReaderEventHandler(Seqdb& aSeqdb)
        : mSeqdb(aSeqdb)
        {
            mHandler.emplace(new DocRootHandler(mSeqdb));
        }

    inline bool StartObject() { return handler(&HandlerBase::StartObject); }
    inline bool EndObject(rapidjson::SizeType /*memberCount*/) { return handler(&HandlerBase::EndObject); }
    inline bool StartArray() { return handler(&HandlerBase::StartArray); }
    inline bool EndArray(rapidjson::SizeType /*elementCount*/) { return handler(&HandlerBase::EndArray); }
    inline bool Key(const char* str, rapidjson::SizeType length, bool /*copy*/) { return handler(&HandlerBase::Key, str, length); }
    inline bool String(const Ch* str, rapidjson::SizeType length, bool /*copy*/) { return handler(&HandlerBase::String, str, length); }
    inline bool Int(int i) { return handler(&HandlerBase::Int, i); }
    inline bool Uint(unsigned u) { return handler(&HandlerBase::Uint, u); }
    inline bool Double(double d) { return handler(&HandlerBase::Double, d); }

      // inline bool Bool(bool /*b*/) { return false; }
      // inline bool Null() { std::cout << "Null()" << std::endl; return false; }
      // inline bool Int64(int64_t i) { std::cout << "Int64(" << i << ")" << std::endl; return false; }
      // inline bool Uint64(uint64_t u) { std::cout << "Uint64(" << u << ")" << std::endl; return false; }

 private:
    Seqdb& mSeqdb;
    std::stack<std::unique_ptr<HandlerBase>> mHandler;

      // ----------------------------------------------------------------------

    // friend void seqdb_import(std::string, Seqdb&);
};

// ----------------------------------------------------------------------

void seqdb_import(std::string buffer, Seqdb& aSeqdb)
{
    if (buffer == "-")
        buffer = acmacs_base::read_stdin();
    else
        buffer = acmacs_base::read_file(buffer);
    if (buffer[0] == '{') { // && buffer.find("\"  version\": " + SEQDB_JSON_DUMP_VERSION) != std::string::npos) {
        SeqdbReaderEventHandler handler{aSeqdb};
        rapidjson::Reader reader;
        rapidjson::StringStream ss(buffer.c_str());
        reader.Parse(ss, handler);
        if (reader.HasParseError())
            throw Error("cannot import seqdb: data parsing failed at pos " + std::to_string(reader.GetErrorOffset()) + ": " +  GetParseError_En(reader.GetParseErrorCode()) + "\n" + buffer.substr(reader.GetErrorOffset(), 50));
    }
    else
        throw std::runtime_error("cannot import seqdb: unrecognized source format");

} // seqdb_import

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
