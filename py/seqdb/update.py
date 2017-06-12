# -*- Python -*-
# license
# license.
# ----------------------------------------------------------------------

import os, re, pprint
import logging; module_logger = logging.getLogger(__name__)
from pathlib import Path
from acmacs_base.files import backup_file
from . import Seqdb, fasta as fasta_m

# ----------------------------------------------------------------------

def create(hidb_dir, seqdb_filename, fasta_files, match_hidb, add_clades, save, report_all_passages, report_identical, report_not_aligned_prefixes, verbose):
    db = Seqdb(str(hidb_dir))
    db_updater = SeqdbUpdater(db, filename=seqdb_filename, load=False)
    for filename in fasta_files:
        data = fasta_m.read_fasta_with_name_parsing(fasta_file=filename, lab="", virus_type="")
        # module_logger.info('{} entries to update seqdb with'.format(len(data)))
        # pprint.pprint(data)
        db_updater.add(data)
    # module_logger.info('Sequences: {} Entries: {}'.format(db.number_of_seqs(), db.number_of_entries()))
    db_updater.detect_insertions_deletions()
    if report_all_passages:
        passages = db.all_passages()
        module_logger.info('Passages: {}\n  {}'.format(len(passages), "\n  ".join(passages)))
    if match_hidb:
        db_updater.match_hidb(verbose=verbose)
        db.build_hi_name_index()          # to report duplicates
    if add_clades:
        db_updater.add_clades()               # clades must be updated after matching with hidb, because matching provides info about B lineage
    # print(db.report())
    if report_identical:
        print(db.report_identical())
    if report_not_aligned_prefixes:
        print(db.report_not_aligned(report_not_aligned_prefixes))
    print(db.report())
    if save:
        db_updater.save(indent=1)

# ----------------------------------------------------------------------

class SeqdbUpdater:

    def __init__(self, seqdb, filename, load=False):
        self.seqdb = seqdb
        self.filename = Path(filename)
        # self.normalize_names = normalize_names
        if load:
            self.load()

    def load(self):
        if self.filename.is_file():
            self.seqdb.load(filename=str(self.filename))

    def save(self, indent):
        if self.filename.is_file():
            backup_file(self.filename)
        self.seqdb.save(filename=str(self.filename), indent=indent)

    # def set_hidb(self, hidb):
    #     self.hidb = hidb

    def add(self, data):
        """data is list of dicts {"date":, "lab":, "name":, "passage":, "reassortant":, "sequence":, "virus_type":, "gene":}"""
        num_added = 0
        for entry in data:
            if entry.get("name"):
                entry["name"] = self.fix_name(entry["name"])
                from .normalize import passage, reassortant
                if entry.get("passage"):
                    entry["passage"] = passage(entry["passage"])
                if entry.get("reassortant"):
                    entry["reassortant"] = reassortant(entry["reassortant"])
                # annotatitions?
                # module_logger.debug('before add_sequence {}'.format({k:v for k,v in entry.items() if k != "sequence"}))
                message = self.seqdb.add_sequence(name=entry["name"], virus_type=entry.get("virus_type", ""), lineage=entry.get("lineage", ""),
                                            lab=entry.get("lab", ""), date=entry.get("date", ""), lab_id=entry.get("lab_id", ""),
                                            passage=entry.get("passage", ""), reassortant=entry.get("reassortant", ""),
                                            sequence=entry.get("sequence", ""), gene=entry.get("gene", ""))
                if message:
                    module_logger.warning("{}: {}".format(entry["name"], message.replace("\n", " ")))
            else:
                module_logger.warning('Cannot add entry without name: {}'.format(entry["lab_id"]))
        messages = self.seqdb.cleanup(remove_short_sequences=True)
        if messages:
            module_logger.warning(messages)

    def add_clades(self):
        if self.seqdb.number_of_seqs():
            for entry_seq in self.seqdb.iter_seq():
                entry_seq.seq.update_clades(virus_type=entry_seq.entry.virus_type, lineage=entry_seq.entry.lineage)

    def match_hidb(self, verbose=False):
        self.seqdb.remove_hi_names()
        self.seqdb.match_hidb(verbose=verbose)

    def detect_insertions_deletions(self):
        self.seqdb.detect_insertions_deletions()

    # ----------------------------------------------------------------------
    # convert

    sNameReplacements = [
        [re.compile(r"CAPETOWN", re.I), "CAPE TOWN"],
        [re.compile(r"(/IRAN/\d+/2015)\s*\(\d.*", re.I), r"\1"],
        [re.compile(r"/MALI\s*0*(\d+)\s*CI/", re.I), r"/MALI/\1 CI/"],
        [re.compile(r"/ST[\.\-\s]*PETERSBURG/", re.I), r"/SAINT PETERSBURG/"],
        [re.compile(r"/MANAGUA/(\d+)_(\d+)/", re.I), r"/MANAGUA/\1.\2/"],
        [re.compile(r"/NAKHONRATCHAISMA/", re.I), r"/NAKHON RATCHASIMA/"],
        [re.compile(r"\s+\([\d\-]+\)$", re.I), r""],   # "A/EGYPT/1323/2016 (96)", "B/ICELAND/33/2015 (15-02274)"
        [re.compile(r"\s+\(VS\d+\)$", re.I), r""],   # "B/BELGIUM/G668/2014 (VS0246)"
        [re.compile(r"/(\d+)_LIKE", re.I), r"/\1"],   # "A/WUHAN/395/95_LIKE"
        [re.compile(r"/SICHUAN[/\-]QINGYANG/?(\d+)/", re.I), r"/SICHUAN QINGYANG/\1/"],
        [re.compile(r"/GENOVA/", re.I), r"/GENOA/"],
        [re.compile(r"/JILILN-NANGUAN/", re.I), r"/JILILN NANGUAN/"],
        [re.compile(r"/NEW MEXCIO/", re.I), r"/NEW MEXICO/"],
        [re.compile(r"/YUNNAN-MENGZI/", re.I), r"/YUNNAN MENGZI/"],
        [re.compile(r"/PTO MONTT/?(\d+)/", re.I), r"/PUERTO MONTT/\1/"],
        [re.compile(r"/PERTH/16/2009 V0152-14\d", re.I), r"/PERTH/16/2009"], # reassortant?
        [re.compile(r"/JILIL?N[\-\s]*NANGUAN/", re.I), r"/JILIN NANGUAN/"],
        [re.compile(r"B/PHUKET/3073/2013 BVR-1B", re.I), r"B/PHUKET/3073/2013"], # reassortant?
        [re.compile(r"/NEIDERSACHSEN/", re.I), r"/NIEDERSACHSEN/"],
        # [re.compile(r"", re.I), r""],
        # [re.compile(r"", re.I), r""],

        # CDC-LV is reassortant?
        [re.compile(r"[\s\-]+CDC-LV\d+[A-Z]?$", re.I), ""],   # [A(H1N1)/SOUTH AFRICA/3626/2013 CDC-LV14A]
        ]

    def fix_name(self, name):
        for repl in self.sNameReplacements:
            name = repl[0].sub(repl[1], name)
        return name

    # ----------------------------------------------------------------------
    # ----------------------------------------------------------------------

    # # sReName = re.compile(r"^(A\(H\d+N\d+\)|B)/")

    # def add(self, data):
    #     """data is list of dicts {"date":, "lab":, "name":, "passage":, "sequence":, "virus_type":, "gene":}"""
    #     num_added = 0
    #     for entry in data:
    #         self._add_sequence(entry)
    #     messages = self.seqdb.cleanup(remove_short_sequences=True)
    #     if messages:
    #         module_logger.warning(messages)

    # def _add_sequence(self, data):
    #     self._normalize(data)
    #     # module_logger.debug('{}'.format({k:v for k,v in data.items() if k not in ["sequence"]}))
    #     name = data.get("name")
    #     if name:
    #         if name[1] in ["/", "("] and data.get("virus_type") and name[0] != data["virus_type"][0]:
    #             module_logger.warning('Virus type ({}) and name ({}) mismatch'.format(data["virus_type"], name))
    #         entry = self.seqdb.find_by_name(name)
    #         if entry is None:
    #             # if self.normalize_names and not self.sReName.match(name):
    #             #     module_logger.warning('Suspicious name {!r}'.format(name))
    #             entry = self.seqdb.new_entry(name)
    #             if data.get("virus_type"):
    #                 entry.virus_type = data["virus_type"] # entry = {"N": name, "s": [], "v": data["virus_type"]}
    #         self._update_db_entry(entry, data)
    #     else:
    #         module_logger.warning('Cannot add entry without name: {}'.format(data["lab_id"]))

    # def _update_db_entry(self, entry, data):
    #     if data.get("virus_type") and entry.virus_type != data["virus_type"]:
    #         raise RuntimeError("Cannot add {!r} to {!r} for {!r}".format(data["virus_type"], entry.virus_type, data.get("name")))
    #     if data.get("location", {}).get("country"):
    #         entry.country = data["location"]["country"]
    #     if data.get("location", {}).get("continent"):
    #         entry.continent = data["location"]["continent"]
    #     if data.get("date"):
    #         entry.add_date(data["date"])
    #     if data.get("annotatitions"):
    #         module_logger.warning('Sequence {} has annotatitions {}'.format(data["name"], data["annotatitions"]))
    #     message = entry.add_or_update_sequence(sequence=data["sequence"], passage=data.get("passage", ""), reassortant=data.get("reassortant", ""), lab=data.get("lab", ""), lab_id=data.get("lab_id", ""), gene=data.get("gene", ""))
    #     if message:
    #         module_logger.warning("{}: {}".format(data["name"], message.replace("\n", " ")))
    #     # if self.hidb:
    #     #     self._match_hidb(entry, data)

    # # ----------------------------------------------------------------------

    # def _normalize(self, data):
    #     from .normalize import passage, reassortant
    #     if data.get("passage"):
    #         data["passage"] = passage(data["passage"])
    #     if data.get("reassortant"):
    #         data["reassortant"] = reassortant(data["reassortant"])
    #     if data.get("name") and data.get("virus_type"):
    #         self._add_virus_type_to_name(data)
    #         data["name"] = data["name"].replace("/LYON/CHU/", "/LYON CHU/").replace(" /", "/").replace("/ ", "/")

    # # ----------------------------------------------------------------------

    # def _add_virus_type_to_name(self, data):
    #     # name can be in different subtypes, e.g. there is A/DAKAR/03/2014 in H1pdm and H3 (both NIMR)
    #     if data["name"][:2] == "A/":
    #         data["name"] = data["virus_type"] + data["name"][1:]

    # ----------------------------------------------------------------------
    # ----------------------------------------------------------------------

    # def _match_hidb(self, seqdb_entry):
    #     hi_entry = self.hidb.find_antigen_by_name(seqdb_entry.name, virus_type=seqdb_entry.virus_type)
    #     # if "3000384319" in seqdb_entry.cdcids():
    #     #     module_logger.debug('M1 {} {} -> {}'.format(seqdb_entry.name, seqdb_entry.cdcids(), hi_entry))
    #     if not hi_entry: # and seqdb_member.seq.has_lab("CDC"):
    #         # module_logger.debug('find by cdcid {} {}'.format(seqdb_entry.name, seqdb_entry.cdcids()))
    #         for cdcid in seqdb_entry.cdcids():
    #             hi_entry = self.hidb.find_antigen_by_cdcid(cdcid, virus_type=seqdb_entry.virus_type)
    #             if hi_entry:
    #                 break
    #         if hi_entry and seqdb_entry.name == "A(H3N2)/TEXAS/88/2016":
    #             module_logger.debug('by cdcid {} {} -> {}'.format(seqdb_entry.name, seqdb_entry.cdcids(), hi_entry["N"]))
    #     if hi_entry:
    #         try:
    #             matches = self.hidb.match(name=seqdb_entry.name,
    #                                       seq_passages_reassortant=[{"p": seq.passages or [""], "r": seq.reassortant or [""]} for seq in seqdb_entry],
    #                                       hi_entry=hi_entry)
    #         except:
    #             module_logger.error('hi_entry {}'.format(hi_entry))
    #             #module_logger.error('e2l {}'.format(e2l))
    #             raise
    #         # if "B/INDIANA/25/2015" in seqdb_entry.name: # "RV2366" in name:
    #         #     module_logger.debug('{}\n{}\nseq_passages_reassortant:{}\nhi_entry:{}\n'.format(seqdb_entry.name, pprint.pformat(matches), pprint.pformat([{"p": seq.passages or [""], "r": seq.reassortant or [""]} for seq in seqdb_entry]), pprint.pformat(hi_entry, width=200)))
    #         if matches:
    #             self._apply_matches(matches, seqdb_entry.name, seqdb_entry, hi_entry)

    # def _apply_matches(self, matches, name, seqdb_entry, hi_entry):
    #     # matches is list of dicts {"s": seq_passage, "r": seq_reassortant, "h": hi_variant}
    #     for m in matches:
    #         for hi_name in [self._make_hi_name(name, m["h"]), self._make_hi_name(hi_entry["N"], m["h"])] + [self._make_hi_name(other_name, m["h"]) for other_name in hi_entry.get("o", [])]:
    #             # if "TEXAS/88/2016" in name:
    #             #     module_logger.debug('{} hi_name: {} hi-entry: {}\n{}'.format(name, hi_name, hi_entry["N"], m))
    #             if hi_name not in self.hidb_already_matched:
    #                 for seq in seqdb_entry:
    #                     if m["p"] in seq.passages and set(m["r"]) & set(seq.reassortant):
    #                         self.hidb_already_matched.add(hi_name)
    #                         seq.add_hi_name(hi_name)
    #                         msg = seqdb_entry.update_lineage(lineage=hi_entry.get("l", ""))
    #                         if msg:
    #                             module_logger.warning('Lineage mismatch for {} {}: seq:{} hi:{}'.format(name, hi_name, seqdb_entry.lineage, hi_entry.get("l", "")))
    #                         for date in hi_entry.get("d", []):
    #                             seqdb_entry.add_date(date)
    #                         break

    # # sVariantFieldOrder = [hidb_m.sReassortantKey, hidb_m.sSerumIdKey, hidb_m.sPassageKey, hidb_m.sAnnotationsKey, hidb_m.sExtraKey, hidb_m.sSerumSpeciesKey]

    # def _make_hi_name(self, name, hi_variant):
    #     return name + " " + " ".join(hi_variant[f] for f in self.sVariantFieldOrder if hi_variant.get(f))

    # ----------------------------------------------------------------------

    # def _normalize(self, data):
    #     if self.normalize_names:
    #         normalized_by_acmacs = self._normalize_by_acmacs(data) # {k: self._normalize_x(k, data) for k in ["name", "passage", "date"]}
    #         # module_logger.info('normalized_by_acmacs\n{}'.format(pprint.pformat(normalized_by_acmacs, width=200)))
    #         for entry in data:
    #             if entry.get("name"):
    #                 norm = normalized_by_acmacs["name"][entry["name"]]
    #                 for err in norm.get("errors", []):
    #                     if "[ag_name] extra:" in str(err):
    #                         module_logger.info("Name parsing for {}: {}".format(entry["name"], err))
    #                     else:
    #                         module_logger.warning("Name parsing for {}: {}".format(entry["name"], err))
    #                 if norm["name"][:2] == "A/" and entry["virus_type"][0] == "A":
    #                     norm["name"] = "{}{}".format(entry["virus_type"], norm["name"][1:])
    #                 elif norm["name"][1] != "/" and norm["name"].count("/") == 2:
    #                     norm["name"] = "{}/{}".format(entry["virus_type"], norm["name"])
    #                 if norm["name"] != entry["name"]:
    #                     entry["raw_name"] = entry["name"]
    #                 entry.update({f: v for f, v in norm.items() if f not in {"errors"}})
    #             for key in ["passage", "date"]:
    #                 if entry.get(key):
    #                     norm = normalized_by_acmacs[key][entry[key]]
    #                     if norm != entry[key]:
    #                         entry["raw_" + key] = entry[key]
    #                         entry[key] = norm

    # def _normalize_by_acmacs(self, data):
    #     return {
    #         "name": acmacs.normalize_names({e["name"] for e in data if e.get("name")}),
    #         "passage": acmacs.normalize_passages({e["passage"] for e in data if e.get("passage")}),
    #         "date": acmacs.normalize_dates({e["date"] for e in data if e.get("date")})
    #         }


# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
