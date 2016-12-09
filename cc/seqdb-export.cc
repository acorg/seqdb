#include <stack>

#include "seqdb-export.hh"
#include "json-writer.hh"
#include "acmacs-base/json-reader.hh"

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

void seqdb_export(std::string aFilename, const Seqdb& aSeqdb, size_t aIndent)
{
    export_to_json(aSeqdb, "seqdb", aFilename, aIndent);
}

// ----------------------------------------------------------------------

typedef json_reader::HandlerBase<Seqdb> HandlerBase;
typedef json_reader::StringListHandler<Seqdb> StringListHandler;
typedef json_reader::MapListHandler<Seqdb> MapListHandler;

// ----------------------------------------------------------------------

// class Error : public std::runtime_error { public: using std::runtime_error::runtime_error; };
// class Failure : public std::exception { public: using std::exception::exception; };
// class Pop : public std::exception { public: using std::exception::exception; };

// class HandlerBase
// {
//  public:
//     inline HandlerBase(Seqdb& aSeqdb) : mSeqdb(aSeqdb), mIgnore(false) {}
//     virtual ~HandlerBase();

//     inline virtual HandlerBase* StartObject() { std::cerr << "HandlerBase StartObject " << typeid(*this).name() << std::endl; throw Failure(); }
//     inline virtual HandlerBase* EndObject() { throw Pop(); }
//     inline virtual HandlerBase* StartArray() { throw Failure(); }
//     inline virtual HandlerBase* EndArray() { throw Pop(); }
//     inline virtual HandlerBase* Double(double d) { std::cerr << "Double: " << d << std::endl; throw Failure(); }
//     inline virtual HandlerBase* Int(int i) { std::cerr << "Int: " << i << std::endl; throw Failure(); }
//     inline virtual HandlerBase* Uint(unsigned u) { std::cerr << "Uint: " << u << std::endl; throw Failure(); }

//     inline virtual HandlerBase* Key(const char* str, rapidjson::SizeType length)
//         {
//             if ((length == 1 && *str == '_') || (length > 0 && *str == '?')) {
//                 mIgnore = true;
//             }
//             else {
//                 std::cerr << "Key: \"" << std::string(str, length) << '"' << std::endl;
//                 throw Failure();
//             }
//             return nullptr;
//         }

//     inline virtual HandlerBase* String(const char* str, rapidjson::SizeType length)
//         {
//             if (mIgnore) {
//                 mIgnore = false;
//             }
//             else {
//                 std::cerr << "String: \"" << std::string(str, length) << '"' << std::endl;
//                 throw Failure();
//             }
//             return nullptr;
//         }

//  protected:
//     Seqdb& mSeqdb;
//     bool mIgnore;
// };

// HandlerBase::~HandlerBase()
// {
// }


// ----------------------------------------------------------------------

// class StringListHandler : public HandlerBase
// {
//  public:
//     inline StringListHandler(Seqdb& aSeqdb, std::vector<std::string>& aList)
//         : HandlerBase(aSeqdb), mList(aList), mStarted(false) {}

//     inline virtual HandlerBase* StartArray()
//         {
//             if (mStarted)
//                 throw Failure();
//             mStarted = true;
//             return nullptr;
//         }

//     inline virtual HandlerBase* EndObject() { throw Failure(); }

//     inline virtual HandlerBase* String(const char* str, rapidjson::SizeType length)
//         {
//             mList.emplace_back(str, length);
//             return nullptr;
//         }

//  private:
//     std::vector<std::string>& mList;
//     bool mStarted;

// }; // class StringListHandler

// // ----------------------------------------------------------------------

// class MapListHandler : public HandlerBase
// {
//  public:
//     inline MapListHandler(Seqdb& aSeqdb, std::map<std::string, std::vector<std::string>>& aMap)
//         : HandlerBase(aSeqdb), mMap(aMap), mStarted(false) {}

//     inline virtual HandlerBase* StartObject()
//         {
//             if (mStarted)
//                 throw Failure();
//             mStarted = true;
//             return nullptr;
//         }

//     inline virtual HandlerBase* EndArray() { throw Failure(); }

//     inline virtual HandlerBase* Key(const char* str, rapidjson::SizeType length)
//         {
//             return new StringListHandler(mSeqdb, mMap[{str, length}]);
//         }

//  private:
//     std::map<std::string, std::vector<std::string>>& mMap;
//     bool mStarted;

// }; // class MapListHandler

// // ----------------------------------------------------------------------

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

// class DocRootHandler : public HandlerBase
// {
//  public:
//     inline DocRootHandler(Seqdb& aSeqdb) : HandlerBase(aSeqdb) {}

//     inline virtual HandlerBase* StartObject() { return new SeqdbRootHandler(mSeqdb); }
// };

// // ----------------------------------------------------------------------

// class SeqdbReaderEventHandler : public rapidjson::BaseReaderHandler<rapidjson::UTF8<>, SeqdbReaderEventHandler>
// {
//  private:
//     template <typename... Args> inline bool handler(HandlerBase* (HandlerBase::*aHandler)(Args... args), Args... args)
//         {
//             try {
//                 auto new_handler = ((*mHandler.top()).*aHandler)(args...);
//                 if (new_handler)
//                     mHandler.emplace(new_handler);
//             }
//             catch (Pop&) {
//                 if (mHandler.empty())
//                     return false;
//                 mHandler.pop();
//             }
//             catch (Failure&) {
//                 return false;
//             }
//             return true;
//         }

//  public:
//     inline SeqdbReaderEventHandler(Seqdb& aSeqdb)
//         : mSeqdb(aSeqdb)
//         {
//             mHandler.emplace(new DocRootHandler(mSeqdb));
//         }

//     inline bool StartObject() { return handler(&HandlerBase::StartObject); }
//     inline bool EndObject(rapidjson::SizeType /*memberCount*/) { return handler(&HandlerBase::EndObject); }
//     inline bool StartArray() { return handler(&HandlerBase::StartArray); }
//     inline bool EndArray(rapidjson::SizeType /*elementCount*/) { return handler(&HandlerBase::EndArray); }
//     inline bool Key(const char* str, rapidjson::SizeType length, bool /*copy*/) { return handler(&HandlerBase::Key, str, length); }
//     inline bool String(const Ch* str, rapidjson::SizeType length, bool /*copy*/) { return handler(&HandlerBase::String, str, length); }
//     inline bool Int(int i) { return handler(&HandlerBase::Int, i); }
//     inline bool Uint(unsigned u) { return handler(&HandlerBase::Uint, u); }
//     inline bool Double(double d) { return handler(&HandlerBase::Double, d); }

//       // inline bool Bool(bool /*b*/) { return false; }
//       // inline bool Null() { std::cout << "Null()" << std::endl; return false; }
//       // inline bool Int64(int64_t i) { std::cout << "Int64(" << i << ")" << std::endl; return false; }
//       // inline bool Uint64(uint64_t u) { std::cout << "Uint64(" << u << ")" << std::endl; return false; }

//  private:
//     Seqdb& mSeqdb;
//     std::stack<std::unique_ptr<HandlerBase>> mHandler;

//       // ----------------------------------------------------------------------

//     // friend void seqdb_import(std::string, Seqdb&);
// };

typedef json_reader::ReaderEventHandler<Seqdb, SeqdbRootHandler> SeqdbReaderEventHandler;


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
            throw json_reader::Error("cannot import seqdb: data parsing failed at pos " + std::to_string(reader.GetErrorOffset()) + ": " +  GetParseError_En(reader.GetParseErrorCode()) + "\n" + buffer.substr(reader.GetErrorOffset(), 50));
    }
    else
        throw json_reader::Error("cannot import seqdb: unrecognized source format");

} // seqdb_import

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
