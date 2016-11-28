#include "acmacs-base/pybind11.hh"
#include "seqdb.hh"

// ----------------------------------------------------------------------

struct PySeqdbEntrySeqIterator
{
    inline PySeqdbEntrySeqIterator(Seqdb& aSeqdb, py::object aRef)
        : mRef(aRef), mCurrent(aSeqdb.begin()), mEnd(aSeqdb.end())
        {
        }

    inline SeqdbEntrySeq next()
        {
            if (mCurrent == mEnd)
                throw py::stop_iteration();
            auto r = *mCurrent;
            ++mCurrent;
            return r;
        }

    inline PySeqdbEntrySeqIterator& filter_lab(std::string aLab) { mCurrent.filter_lab(aLab); return *this; }
    inline PySeqdbEntrySeqIterator& filter_labid(std::string aLab, std::string aId) { mCurrent.filter_labid(aLab, aId); return *this; }
    inline PySeqdbEntrySeqIterator& filter_subtype(std::string aSubtype) { mCurrent.filter_subtype(aSubtype); return *this; }
    inline PySeqdbEntrySeqIterator& filter_lineage(std::string aLineage) { mCurrent.filter_lineage(aLineage); return *this; }
    inline PySeqdbEntrySeqIterator& filter_aligned(bool aAligned) { mCurrent.filter_aligned(aAligned); return *this; }
    inline PySeqdbEntrySeqIterator& filter_gene(std::string aGene) { mCurrent.filter_gene(aGene); return *this; }
    inline PySeqdbEntrySeqIterator& filter_date_range(std::string aBegin, std::string aEnd) { mCurrent.filter_date_range(aBegin, aEnd); return *this; }
    inline PySeqdbEntrySeqIterator& filter_hi_name(bool aHasHiName) { mCurrent.filter_hi_name(aHasHiName); return *this; }
    inline PySeqdbEntrySeqIterator& filter_name_regex(std::string aNameRegex) { mCurrent.filter_name_regex(aNameRegex); return *this; }

    py::object mRef; // keep a reference
    SeqdbIterator mCurrent;
    SeqdbIterator mEnd;

}; // struct PySeqdbSeqIterator

// ----------------------------------------------------------------------

struct PySeqdbEntryIterator
{
    inline PySeqdbEntryIterator(Seqdb& aSeqdb, py::object aRef)
        : mRef(aRef), mCurrent(aSeqdb.begin_entry()), mEnd(aSeqdb.end_entry())
        {
        }

    inline auto& next()
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
    inline PySeqdbSeqIterator(SeqdbEntry& aSeqdbEntry)
        : mCurrent(aSeqdbEntry.begin_seq()), mEnd(aSeqdbEntry.end_seq())
        {
        }

    inline auto& next()
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

PYBIND11_PLUGIN(seqdb_backend)
{
    py::module m("seqdb_backend", "Seqdb access plugin");

      // ----------------------------------------------------------------------

    py::class_<SeqdbSeq>(m, "SeqdbSeq")
            .def("has_lab", &SeqdbSeq::has_lab, py::arg("lab"))
            .def("cdcids", &SeqdbSeq::cdcids)
            .def("update_clades", &SeqdbSeq::update_clades, py::arg("virus_type"), py::arg("lineage"))
            .def_property_readonly("passages", static_cast<const std::vector<std::string>& (SeqdbSeq::*)() const>(&SeqdbSeq::passages))
            .def_property_readonly("reassortant", static_cast<const std::vector<std::string>& (SeqdbSeq::*)() const>(&SeqdbSeq::reassortant))
            .def_property_readonly("hi_names", static_cast<const std::vector<std::string>& (SeqdbSeq::*)() const>(&SeqdbSeq::hi_names))
            .def("add_hi_name", &SeqdbSeq::add_hi_name, py::arg("hi_name"))
            .def("amino_acids", static_cast<std::string (SeqdbSeq::*)(bool, size_t) const>(&SeqdbSeq::amino_acids), py::arg("aligned"), py::arg("left_part_size") = int(0), py::doc("if aligned and left_part_size > 0 - include signal peptide and other stuff to the left from the beginning of the aligned sequence."))
            .def("nucleotides", static_cast<std::string (SeqdbSeq::*)(bool, size_t) const>(&SeqdbSeq::nucleotides), py::arg("aligned"), py::arg("left_part_size") = int(0), py::doc("if aligned and left_part_size > 0 - include signal peptide and other stuff to the left from the beginning of the aligned sequence."))
            .def("amino_acids_shift", &SeqdbSeq::amino_acids_shift)
            .def("nucleotides_shift", &SeqdbSeq::nucleotides_shift)
            .def("lab", &SeqdbSeq::lab)
            .def("lab_id", &SeqdbSeq::lab_id)
            .def("lab_ids", &SeqdbSeq::lab_ids_for_lab, py::arg("lab"))
            .def("lab_ids", &SeqdbSeq::lab_ids)
            .def("passage", &SeqdbSeq::passage)
            .def("gene", static_cast<std::string (SeqdbSeq::*)() const>(&SeqdbSeq::gene))
            .def("clades", static_cast<const std::vector<std::string>& (SeqdbSeq::*)() const>(&SeqdbSeq::clades))
            ;

    py::class_<SeqdbEntry>(m, "SeqdbEntry")
            .def_property_readonly("name", static_cast<std::string (SeqdbEntry::*)() const>(&SeqdbEntry::name))
            .def_property("continent", static_cast<std::string (SeqdbEntry::*)() const>(&SeqdbEntry::continent), static_cast<void (SeqdbEntry::*)(std::string)>(&SeqdbEntry::continent))
            .def_property("country", static_cast<std::string (SeqdbEntry::*)() const>(&SeqdbEntry::country), static_cast<void (SeqdbEntry::*)(std::string)>(&SeqdbEntry::country))
            .def_property("virus_type", static_cast<std::string (SeqdbEntry::*)() const>(&SeqdbEntry::virus_type), static_cast<void (SeqdbEntry::*)(std::string)>(&SeqdbEntry::virus_type))
              // .def("set_virus_type", static_cast<void (SeqdbEntry::*)(std::string)>(&SeqdbEntry::virus_type), py::arg("virus_type"))
            .def("add_date", &SeqdbEntry::add_date, py::arg("date"))
            .def("add_or_update_sequence", &SeqdbEntry::add_or_update_sequence, py::arg("sequence"), py::arg("passage"), py::arg("reassortant"), py::arg("lab"), py::arg("lab_id"), py::arg("gene"))
            .def("cdcids", &SeqdbEntry::cdcids)
            .def("date", &SeqdbEntry::date)
            .def_property_readonly("lineage", static_cast<std::string (SeqdbEntry::*)() const>(&SeqdbEntry::lineage))
            .def("update_lineage", [](SeqdbEntry& entry, std::string lineage) -> std::string { Messages msg; entry.update_lineage(lineage, msg); return msg; }, py::arg("lineage"))
            .def("__iter__", [](py::object entry) { return PySeqdbSeqIterator(entry.cast<SeqdbEntry&>()); })
            ;

    py::class_<SeqdbEntrySeq>(m, "SeqdbEntrySeq")
            .def_property_readonly("entry", static_cast<const SeqdbEntry& (SeqdbEntrySeq::*)() const>(&SeqdbEntrySeq::entry), py::return_value_policy::reference)
            .def_property_readonly("seq", static_cast<const SeqdbSeq& (SeqdbEntrySeq::*)() const>(&SeqdbEntrySeq::seq), py::return_value_policy::reference)
            .def("make_name", &SeqdbEntrySeq::make_name, py::arg("passage_separator") = std::string(" "))
            .def("seq_id", &SeqdbEntrySeq::seq_id)
            ;

    py::class_<PySeqdbEntrySeqIterator>(m, "PySeqdbEntrySeqIterator")
            .def("__iter__", [](PySeqdbEntrySeqIterator& it) { return it; })
            .def("__next__", &PySeqdbEntrySeqIterator::next)
            .def("filter_lab", &PySeqdbEntrySeqIterator::filter_lab)
            .def("filter_labid", &PySeqdbEntrySeqIterator::filter_labid, py::arg("lab"), py::arg("id"))
            .def("filter_subtype", &PySeqdbEntrySeqIterator::filter_subtype)
            .def("filter_lineage", &PySeqdbEntrySeqIterator::filter_lineage)
            .def("filter_aligned", &PySeqdbEntrySeqIterator::filter_aligned)
            .def("filter_gene", &PySeqdbEntrySeqIterator::filter_gene)
            .def("filter_date_range", &PySeqdbEntrySeqIterator::filter_date_range)
            .def("filter_hi_name", &PySeqdbEntrySeqIterator::filter_hi_name)
            .def("filter_name_regex", &PySeqdbEntrySeqIterator::filter_name_regex)
            ;

    py::class_<PySeqdbEntryIterator>(m, "PySeqdbEntryIterator")
            .def("__iter__", [](PySeqdbEntryIterator& it) { return it; })
            .def("__next__", &PySeqdbEntryIterator::next, py::return_value_policy::reference);
            ;

    py::class_<PySeqdbSeqIterator>(m, "PySeqdbSeqIterator")
            .def("__iter__", [](PySeqdbSeqIterator& it) { return it; })
            .def("__next__", &PySeqdbSeqIterator::next, py::return_value_policy::reference);
            ;

    py::class_<Seqdb>(m, "Seqdb")
            .def(py::init<>())
              // .def("from_json", &Seqdb::from_json, py::doc("reads seqdb from json"))
            .def("load", &Seqdb::load, py::arg("filename") = std::string(), py::doc("reads seqdb from file containing json"))
              //.def("json", &Seqdb::to_json, py::arg("indent") = size_t(0))
            .def("save", &Seqdb::save, py::arg("filename") = std::string(), py::arg("indent") = size_t(0), py::doc("writes seqdb into file in json format"))
            .def("number_of_entries", &Seqdb::number_of_entries)
            .def("find_by_name", static_cast<SeqdbEntry* (Seqdb::*)(std::string)>(&Seqdb::find_by_name), py::arg("name"), py::return_value_policy::reference, py::doc("returns entry found by name or None"))
            .def("new_entry", &Seqdb::new_entry, py::arg("name"), py::return_value_policy::reference, py::doc("creates and inserts into the database new entry with the passed name, returns that name, throws if database already has entry with that name."))
            .def("cleanup", &Seqdb::cleanup, py::arg("remove_short_sequences") = true)
            .def("report", &Seqdb::report)
            .def("report_identical", &Seqdb::report_identical)
            .def("report_not_aligned", &Seqdb::report_not_aligned, py::arg("prefix_size"), py::doc("returns report with AA prefixes of not aligned sequences."))
            .def("iter_seq", [](py::object seqdb) { return PySeqdbEntrySeqIterator(seqdb.cast<Seqdb&>(), seqdb); })
            .def("iter_entry", [](py::object seqdb) { return PySeqdbEntryIterator(seqdb.cast<Seqdb&>(), seqdb); })
            .def("all_hi_names", &Seqdb::all_hi_names, py::doc("returns list of all hi_names (\"h\") found in seqdb."))
            .def("remove_hi_names", &Seqdb::remove_hi_names, py::doc("removes all hi_names (\"h\") found in seqdb (e.g. before matching again)."))
            .def("match_hidb", &Seqdb::match_hidb, py::arg("hidb_dir"), py::doc("match all names against hidb"))
            ;

      // ----------------------------------------------------------------------

    return m.ptr();
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
