#include "seqdb-import.hh"
#include "seqdb/seqdb.hh"
#include "json-keys.hh"

#include "acmacs-base/read-file.hh"
#include "acmacs-base/json-importer.hh"
namespace jsi = json_importer;

// ----------------------------------------------------------------------

namespace seqdb
{
    static constexpr const char* SEQDB_JSON_DUMP_VERSION = "sequence-database-v2";

      // ----------------------------------------------------------------------

    class SeqdbDataFile
    {
     public:
        inline SeqdbDataFile(Seqdb& aSeqdb) : mSeqdb(aSeqdb) {}

        inline void indentation(const char* /*str*/, size_t /*length*/)
            {
                  // mIndentation.assign(str, length);
                  // std::cerr << "Indentation: " << mIndentation << std::endl;
            }

        inline void version(const char* str, size_t length)
            {
                const std::string version{str, length};
                if (version != SEQDB_JSON_DUMP_VERSION)
                    throw std::runtime_error("Unsupported seqdb version: \"" + version + "\"");
            }

        inline std::vector<SeqdbEntry>& seqdb() { return mSeqdb.entries(); }

     private:
        Seqdb& mSeqdb;
          // std::string mIndentation;
    };

      // ----------------------------------------------------------------------

    using LabIds = std::map<std::string, std::vector<std::string>>;

    class LabIdStorer : public jsi::StorerBase
    {
     public:
        using Base = jsi::StorerBase;

        inline LabIdStorer(LabIds& aTarget) : mTarget(aTarget), mStarted(false) {}

    inline virtual Base* StartObject()
        {
            if (mStarted)
                return jsi::storers::_i::failure(typeid(*this).name() + std::string(": unexpected StartObject event"));
            mTarget.clear();
            mStarted = true;
            return nullptr;
        }

    inline virtual Base* EndObject()
        {
            return jsi::storers::_i::pop();
        }

    inline virtual Base* Key(const char* str, rapidjson::SizeType length)
        {
            mLab.assign(str, length);
            return nullptr;
        }

    inline virtual Base* StartArray()
        {
            if (mLab.empty())
                return jsi::storers::_i::failure(typeid(*this).name() + std::string(": unexpected StartArray event"));
            mTarget.emplace(mLab, std::vector<std::string>{});
            return nullptr;
        }

    inline virtual Base* EndArray()
        {
            mLab.clear();
            return nullptr;
        }

    inline virtual Base* String(const char* str, rapidjson::SizeType length)
        {
            if (mLab.empty())
                return jsi::storers::_i::failure(typeid(*this).name() + std::string(": unexpected String event"));
            mTarget[mLab].emplace_back(str, length);
            return nullptr;
        }

     private:
        LabIds& mTarget;
        bool mStarted;
        std::string mLab;
    };

    void seqdb_import(std::string aFilename, Seqdb& aSeqdb)
    {
        jsi::data<SeqdbSeq> seq_data = {
            {"a", jsi::field(&SeqdbSeq::amino_acids)},
            {"c", jsi::field(&SeqdbSeq::clades)},
            {"g", jsi::field(&SeqdbSeq::gene)},
            {"h", jsi::field(&SeqdbSeq::hi_names)},
            {"l", jsi::field<LabIdStorer, SeqdbSeq, LabIds>(&SeqdbSeq::lab_ids_raw)}, // {"lab": ["lab_id"]},
            {"n", jsi::field(&SeqdbSeq::nucleotides)},
            {"p", jsi::field(&SeqdbSeq::passages)},
            {"r", jsi::field(&SeqdbSeq::reassortant)},
            {"s", jsi::field(&SeqdbSeq::amino_acids_shift_raw)},
            {"t", jsi::field(&SeqdbSeq::nucleotides_shift_raw)},
        };

        using ESS = void (SeqdbEntry::*)(const char*, size_t);
        jsi::data<SeqdbEntry> entry_data = {
            {"N", jsi::field(static_cast<ESS>(&SeqdbEntry::name))},
            {"d", jsi::field(&SeqdbEntry::dates)},
            {"C", jsi::field(static_cast<ESS>(&SeqdbEntry::continent))},
            {"c", jsi::field(static_cast<ESS>(&SeqdbEntry::country))},
            {"l", jsi::field(static_cast<ESS>(&SeqdbEntry::lineage))},
            {"s", jsi::field(&SeqdbEntry::seqs, seq_data)},
            {"v", jsi::field(static_cast<ESS>(&SeqdbEntry::virus_type))},
        };

        jsi::data<SeqdbDataFile> seqdb_data = {
            {"_", jsi::field(&SeqdbDataFile::indentation)},
            {"  version", jsi::field(&SeqdbDataFile::version)},
            {"data", jsi::field(&SeqdbDataFile::seqdb, entry_data)},
        };

        SeqdbDataFile data{aSeqdb};
        jsi::import(acmacs::file::read(aFilename), data, seqdb_data);

    } // seqdb_import
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
