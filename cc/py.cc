#include "acmacs-base/pybind11.hh"
#include "acmacs-chart-2/chart.hh"
#include "seqdb/seqdb.hh"

using namespace seqdb;

// ----------------------------------------------------------------------

struct PySeqdbEntrySeqIterator
{
    PySeqdbEntrySeqIterator(Seqdb& aSeqdb, py::object aRef)
        : mRef(aRef), mCurrent(aSeqdb.begin()), mEnd(aSeqdb.end())
        {
        }

    SeqdbEntrySeq next()
        {
            if (mCurrent == mEnd)
                throw py::stop_iteration();
            auto r = *mCurrent;
            ++mCurrent;
            return r;
        }

    PySeqdbEntrySeqIterator& filter_lab(std::string aLab) { mCurrent.filter_lab(aLab); return *this; }
    PySeqdbEntrySeqIterator& filter_labid(std::string aLab, std::string aId) { mCurrent.filter_labid(aLab, aId); return *this; }
    PySeqdbEntrySeqIterator& filter_subtype(std::string aSubtype) { mCurrent.filter_subtype(aSubtype); return *this; }
    PySeqdbEntrySeqIterator& filter_lineage(std::string aLineage) { mCurrent.filter_lineage(aLineage); return *this; }
    PySeqdbEntrySeqIterator& filter_continent(std::string aContinent) { mCurrent.filter_continent(aContinent); return *this; }
    PySeqdbEntrySeqIterator& filter_country(std::string aCountry) { mCurrent.filter_country(aCountry); return *this; }
    PySeqdbEntrySeqIterator& filter_aligned(bool aAligned) { mCurrent.filter_aligned(aAligned); return *this; }
    PySeqdbEntrySeqIterator& filter_gene(std::string aGene) { mCurrent.filter_gene(aGene); return *this; }
    PySeqdbEntrySeqIterator& filter_date_range(std::string aBegin, std::string aEnd) { mCurrent.filter_date_range(aBegin, aEnd); return *this; }
    PySeqdbEntrySeqIterator& filter_hi_name(bool aHasHiName) { mCurrent.filter_hi_name(aHasHiName); return *this; }
    PySeqdbEntrySeqIterator& filter_name_regex(std::string aNameRegex) { mCurrent.filter_name_regex(aNameRegex); return *this; }
    PySeqdbEntrySeqIterator& filter_clade(std::string aClade) { mCurrent.filter_clade(aClade); return *this; }

    py::object mRef; // keep a reference
    SeqdbIterator mCurrent;
    SeqdbIterator mEnd;

}; // struct PySeqdbSeqIterator

// ----------------------------------------------------------------------

struct PySeqdbEntryIterator
{
    PySeqdbEntryIterator(Seqdb& aSeqdb, py::object aRef)
        : mRef(aRef), mCurrent(aSeqdb.begin_entry()), mEnd(aSeqdb.end_entry())
        {
        }

    auto& next()
        {
            if (mCurrent == mEnd)
                throw py::stop_iteration();
            auto& r = *mCurrent;
            ++mCurrent;
            return r;
        }

    py::object mRef; // keep a reference
    std::vector<SeqdbEntry>::iterator mCurrent;
    std::vector<SeqdbEntry>::iterator mEnd;

}; // struct PySeqdbEntryIterator

// ----------------------------------------------------------------------

struct PySeqdbSeqIterator
{
    PySeqdbSeqIterator(SeqdbEntry& aSeqdbEntry)
        : mCurrent(aSeqdbEntry.begin_seq()), mEnd(aSeqdbEntry.end_seq())
        {
        }

    auto& next()
        {
            if (mCurrent == mEnd)
                throw py::stop_iteration();
            auto& r = *mCurrent;
            ++mCurrent;
            return r;
        }

    std::vector<SeqdbSeq>::iterator mCurrent;
    std::vector<SeqdbSeq>::iterator mEnd;

}; // struct PySeqdbSeqIterator

// ----------------------------------------------------------------------

PYBIND11_MODULE(seqdb_backend, m)
{
    m.doc() = "Seqdb access plugin";

      // ----------------------------------------------------------------------

    py::class_<SeqdbSeq>(m, "SeqdbSeq")
            .def("has_lab", &SeqdbSeq::has_lab, py::arg("lab"))
            .def("cdcids", &SeqdbSeq::cdcids)
            .def("update_clades", &SeqdbSeq::update_clades, py::arg("virus_type"), py::arg("lineage"), py::arg("name") = "*name not available*")
            .def_property_readonly("passages", static_cast<const std::vector<std::string>& (SeqdbSeq::*)() const>(&SeqdbSeq::passages))
            .def_property_readonly("reassortant", static_cast<const std::vector<std::string>& (SeqdbSeq::*)() const>(&SeqdbSeq::reassortant))
            .def_property_readonly("hi_names", static_cast<const std::vector<std::string>& (SeqdbSeq::*)() const>(&SeqdbSeq::hi_names))
            .def("add_hi_name", &SeqdbSeq::add_hi_name, py::arg("hi_name"))
            .def("amino_acids", static_cast<std::string (SeqdbSeq::*)(bool, size_t, size_t) const>(&SeqdbSeq::amino_acids), py::arg("aligned"), py::arg("left_part_size") = int(0), py::arg("resize") = int(0), py::doc("if aligned and left_part_size > 0 - include signal peptide and other stuff to the left from the beginning of the aligned sequence."))
            .def("nucleotides", static_cast<std::string (SeqdbSeq::*)(bool, size_t, size_t) const>(&SeqdbSeq::nucleotides), py::arg("aligned"), py::arg("left_part_size") = int(0), py::arg("resize") = int(0), py::doc("if aligned and left_part_size > 0 - include signal peptide and other stuff to the left from the beginning of the aligned sequence."))
            .def("amino_acids_shift", &SeqdbSeq::amino_acids_shift)
            .def("nucleotides_shift", &SeqdbSeq::nucleotides_shift)
            .def("amino_acid_at", &SeqdbSeq::amino_acid_at, py::arg("pos"), py::arg("ignore_errors") = false, py::doc("pos starts from 1"))
            .def("lab", &SeqdbSeq::lab)
            .def("lab_id", &SeqdbSeq::lab_id)
            .def("lab_ids", &SeqdbSeq::lab_ids_for_lab, py::arg("lab"))
            .def("lab_ids", &SeqdbSeq::lab_ids)
            .def("passage", &SeqdbSeq::passage)
            .def("gene", static_cast<std::string (SeqdbSeq::*)() const>(&SeqdbSeq::gene))
            .def("clades", py::overload_cast<>(&SeqdbSeq::clades))
            ;

    py::class_<SeqdbEntry>(m, "SeqdbEntry")
            .def_property_readonly("name", static_cast<std::string_view (SeqdbEntry::*)() const>(&SeqdbEntry::name))
            .def_property_readonly("continent", static_cast<std::string_view (SeqdbEntry::*)() const>(&SeqdbEntry::continent))
            .def_property_readonly("country", static_cast<std::string_view (SeqdbEntry::*)() const>(&SeqdbEntry::country))
            .def_property_readonly("virus_type", static_cast<std::string_view (SeqdbEntry::*)() const>(&SeqdbEntry::virus_type))
              // .def("set_virus_type", static_cast<void (SeqdbEntry::*)(std::string)>(&SeqdbEntry::virus_type), py::arg("virus_type"))
            .def("add_date", &SeqdbEntry::add_date, py::arg("date"))
              // .def("add_or_update_sequence", &SeqdbEntry::add_or_update_sequence, py::arg("sequence"), py::arg("passage"), py::arg("reassortant"), py::arg("lab"), py::arg("lab_id"), py::arg("gene"))
            .def("cdcids", &SeqdbEntry::cdcids)
            .def("date", &SeqdbEntry::date)
            .def_property_readonly("lineage", static_cast<std::string_view (SeqdbEntry::*)() const>(&SeqdbEntry::lineage))
            .def("update_lineage", [](SeqdbEntry& entry, std::string lineage) -> std::string { Messages msg; entry.update_lineage(lineage, msg); return msg; }, py::arg("lineage"))
            .def("__iter__", [](py::object entry) { return PySeqdbSeqIterator(entry.cast<SeqdbEntry&>()); })
            ;

    py::class_<SeqdbEntrySeq>(m, "SeqdbEntrySeq")
            .def_property_readonly("entry", static_cast<const SeqdbEntry& (SeqdbEntrySeq::*)() const>(&SeqdbEntrySeq::entry), py::return_value_policy::reference)
            .def_property_readonly("seq", static_cast<const SeqdbSeq& (SeqdbEntrySeq::*)() const>(&SeqdbEntrySeq::seq), py::return_value_policy::reference)
            .def("make_name", &SeqdbEntrySeq::make_name, py::arg("passage_separator") = std::string(" "))
            .def("seq_id", [](const SeqdbEntrySeq& e, bool encoded) { return e.seq_id(encoded ? seqdb::SeqdbEntrySeq::encoded_t::yes : seqdb::SeqdbEntrySeq::encoded_t::no); }, py::arg("encoded") = false)
            .def("__bool__", &seqdb::SeqdbEntrySeq::operator bool)
            ;

    py::class_<PySeqdbEntrySeqIterator>(m, "PySeqdbEntrySeqIterator")
            .def("__iter__", [](PySeqdbEntrySeqIterator& it) { return it; })
            .def("__next__", &PySeqdbEntrySeqIterator::next)
            .def("filter_lab", &PySeqdbEntrySeqIterator::filter_lab)
            .def("filter_labid", &PySeqdbEntrySeqIterator::filter_labid, py::arg("lab"), py::arg("id"))
            .def("filter_subtype", &PySeqdbEntrySeqIterator::filter_subtype)
            .def("filter_lineage", &PySeqdbEntrySeqIterator::filter_lineage)
            .def("filter_continent", &PySeqdbEntrySeqIterator::filter_continent)
            .def("filter_country", &PySeqdbEntrySeqIterator::filter_country)
            .def("filter_aligned", &PySeqdbEntrySeqIterator::filter_aligned)
            .def("filter_gene", &PySeqdbEntrySeqIterator::filter_gene)
            .def("filter_clade", &PySeqdbEntrySeqIterator::filter_clade)
            .def("filter_date_range", &PySeqdbEntrySeqIterator::filter_date_range)
            .def("filter_hi_name", &PySeqdbEntrySeqIterator::filter_hi_name)
            .def("filter_name_regex", &PySeqdbEntrySeqIterator::filter_name_regex)
            ;

    py::class_<PySeqdbEntryIterator>(m, "PySeqdbEntryIterator")
            .def("__iter__", [](PySeqdbEntryIterator& it) { return it; })
            .def("__next__", &PySeqdbEntryIterator::next, py::return_value_policy::reference)
            ;

    py::class_<PySeqdbSeqIterator>(m, "PySeqdbSeqIterator")
            .def("__iter__", [](PySeqdbSeqIterator& it) { return it; })
            .def("__next__", &PySeqdbSeqIterator::next, py::return_value_policy::reference)
            ;

    py::class_<Seqdb>(m, "Seqdb")
            .def(py::init<>())
              // .def("from_json", &Seqdb::from_json, py::doc("reads seqdb from json"))
              //.def("load", &Seqdb::load, py::arg("filename") = std::string(), py::doc("reads seqdb from file containing json"))
              //.def("json", &Seqdb::to_json, py::arg("indent") = size_t(0))
            .def("save", &Seqdb::save, py::arg("filename") = std::string(), py::arg("indent") = size_t(0), py::doc("writes seqdb into file in json format"))
            .def("number_of_entries", &Seqdb::number_of_entries)
            .def("number_of_seqs", &Seqdb::number_of_seqs)
            .def("find_by_name", static_cast<SeqdbEntry* (Seqdb::*)(std::string_view)>(&Seqdb::find_by_name), py::arg("name"), py::return_value_policy::reference, py::doc("returns entry found by name or None"))
            .def("find_by_seq_id", [](seqdb::Seqdb& aSeqdb, std::string seq_id, bool ignore) { return aSeqdb.find_by_seq_id(seq_id, ignore ? seqdb::Seqdb::ignore_not_found::yes : seqdb::Seqdb::ignore_not_found::no); }, py::arg("seq_id"), py::arg("ignore_not_found") = false, py::return_value_policy::reference, py::doc("returns entry_seq found by seq_id"))
              // .def("new_entry", &Seqdb::new_entry, py::arg("name"), py::return_value_policy::reference, py::doc("creates and inserts into the database new entry with the passed name, returns that name, throws if database already has entry with that name."))
            .def("add_sequence", &Seqdb::add_sequence, py::arg("name"), py::arg("virus_type"), py::arg("lineage"), py::arg("lab"), py::arg("date"), py::arg("lab_id"), py::arg("passage"), py::arg("reassortant"), py::arg("sequence"), py::arg("gene"), py::doc("adds sequence to the database, inserts new entry if necessary."))
            .def("report_not_aligned_after_adding", &Seqdb::report_not_aligned_after_adding, py::doc("writes to stderr list of not aligned sequences recently added via add_sequence."))
            .def("cleanup", &Seqdb::cleanup, py::arg("remove_short_sequences") = true)
            .def("detect_insertions_deletions", &Seqdb::detect_insertions_deletions)
            .def("detect_b_lineage", &Seqdb::detect_b_lineage)
            .def("update_clades", &Seqdb::update_clades)
            .def("report", &Seqdb::report)
            .def("report_identical", &Seqdb::report_identical)
            .def("report_not_aligned", &Seqdb::report_not_aligned, py::arg("prefix_size"), py::doc("returns report with AA prefixes of not aligned sequences."))
            .def("iter_seq", [](py::object seqdb) { return PySeqdbEntrySeqIterator(seqdb.cast<Seqdb&>(), seqdb); })
            .def("iter_entry", [](py::object seqdb) { return PySeqdbEntryIterator(seqdb.cast<Seqdb&>(), seqdb); })
            .def("all_hi_names", &Seqdb::all_hi_names, py::doc("returns list of all hi_names (\"h\") found in seqdb."))
            .def("all_passages", &Seqdb::all_passages, py::doc("returns list of all passages found in seqdb."))
            .def("remove_hi_names", &Seqdb::remove_hi_names, py::doc("removes all hi_names (\"h\") found in seqdb (e.g. before matching again)."))
            .def("match_hidb", [](seqdb::Seqdb& aSeqdb, bool verbose, bool greedy) { aSeqdb.match_hidb(verbose ? seqdb::report::yes : seqdb::report::no, greedy); }, py::arg("verbose") = false, py::arg("greedy") = true, py::doc("match all names against hidb, returns list of not found locations"))
            .def("build_hi_name_index", &Seqdb::build_hi_name_index)
            .def("find_hi_name", &Seqdb::find_hi_name, py::arg("name"), py::return_value_policy::reference, py::doc("returns entry_seq found by hi name or None"))
            .def("aa_at_positions_for_antigens", [](const seqdb::Seqdb& aSeqdb, const acmacs::chart::Antigens& aAntigens, const std::vector<size_t>& aPositions, bool aVerbose) {
                                                     std::map<std::string, std::vector<size_t>> r; aSeqdb.aa_at_positions_for_antigens(aAntigens, aPositions, r, aVerbose ? seqdb::report::yes : seqdb::report::no); return r; }, py::arg("antigens"), py::arg("positions"), py::arg("verbose"))
            .def("match_antigens", [](const seqdb::Seqdb& aSeqdb, const acmacs::chart::Antigens& aAntigens, std::string aChartVirusType, bool aVerbose) {
                                       std::vector<seqdb::SeqdbEntrySeq> r; aSeqdb.match(aAntigens, r, aChartVirusType, aVerbose ? seqdb::report::yes : seqdb::report::no); return r; }, py::arg("antigens"), py::arg("virus_type") = std::string{}, py::arg("verbose") = true)
            ;

    m.def("setup_dbs", [](std::string db_dir, bool aVerbose) { seqdb::setup_dbs(db_dir, aVerbose ? seqdb::report::yes : seqdb::report::no); }, py::arg("db_dir"), py::arg("verbose") = false);
    m.def("seqdb_setup", [](std::string filename, bool aVerbose) { seqdb::setup(filename, aVerbose ? seqdb::report::yes : seqdb::report::no); }, py::arg("filename"), py::arg("verbose") = false);
    m.def("get_seqdb", [](bool aTimer) { return seqdb::get(seqdb::ignore_errors::no, do_report_time(aTimer)); }, py::arg("timer") = false, py::return_value_policy::reference);


}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
