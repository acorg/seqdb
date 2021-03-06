# -*- Python -*-
# license
# license.

"""
Functions for reading and generating fasta files.
"""

import os, re, collections, operator, random as random_m
import logging; module_logger = logging.getLogger(__name__)
from acmacs_base.files import read_text, write_binary
from acmacs_base.encode_name import encode
from . import normalize

# ======================================================================

class FastaReaderError (Exception):
    pass

# ======================================================================

def export_from_seqdb(seqdb, filename, output_format, amino_acids, lab, virus_type, lineage, continent, country, gene, clade, start_date, end_date, random, recent, base_seq, include_seq, name_format, aligned, truncate_left, encode_name, wrap, truncate_to_most_common_length, hamming_distance_threshold, hamming_distance_report, sort_by, with_hi_name, name_match):

    def make_entry(e):
        r = {
            "e": e,
            "n": name_format.format(hi_name_or_seq_name_with_passage=e.make_name(), name=e.entry.name,
                                        date=e.entry.date(), lab_id=e.seq.lab_id(), passage=e.seq.passage(),
                                        lab=e.seq.lab(), gene=e.seq.gene(), seq_id=e.seq_id()),
            "d": e.entry.date(),
            }
        return r

    def get_sequence(e, left_part_size):
        try:
            e["s"] = e["e"].seq.amino_acids(aligned=aligned, left_part_size=left_part_size) if amino_acids else e["e"].seq.nucleotides(aligned=aligned, left_part_size=left_part_size)
        except Exception as err:
            module_logger.warning("Cannot get sequence for {}: {}".format(e["n"], err))
            e["s"] = None
        return e

    def left_part(e):
        return - (e["e"].seq.amino_acids_shift() if amino_acids else e["e"].seq.nucleotides_shift())

    def exclude_by_hamming_distance(e1, e2, threshold):
        hd = hamming_distance(e1["s"], e2["s"], e1["n"], e2["n"])
        if hd >= threshold:
            module_logger.info('{!r} excluded because hamming distance to {!r} is {} (threshold: {})'.format(e2["n"], e1["n"], hd, threshold))
            r = False
        else:
            r = True
        return r

    def seq_present(seqs, e1):
        present = [e2_no for e2_no, e2 in enumerate(seqs) if e1["n"] == e2["n"]]
        return present[0] if present else None

    # ----------------------------------------------------------------------

    if start_date is not None:
        module_logger.info('Requested start date: ' + start_date)
    if end_date is not None:
        module_logger.info('Requested end date:   ' + end_date)

    virus_type, lineage = normalize.virus_type_lineage(virus_type, lineage)
    iter = (seqdb.iter_seq()
            .filter_lab(normalize.lab(lab) or "")
            .filter_subtype(virus_type or "")
            .filter_lineage(lineage or "")
            .filter_continent(continent or "")
            .filter_country(country or "")
            .filter_aligned(aligned)
            .filter_gene(gene)
            .filter_clade(clade or "")
            .filter_date_range(normalize.date(start_date), normalize.date(end_date))
            .filter_hi_name(with_hi_name)
            )
    if name_match is not None:
        iter = iter.filter_name_regex(name_match)
    sequences = [make_entry(e) for e in iter]
    left_part_size = 0 if truncate_left else max(left_part(seq) for seq in sequences)
    if left_part_size:
        module_logger.info('Left part size (signal peptide): {}'.format(left_part_size))
    sequences = [get_sequence(seq, left_part_size) for seq in sequences]

    # report empty sequences
    empty = [seq for seq in sequences if not seq["s"]]
    if empty:
        module_logger.warning('The following sequences are empty and are not exported:\n  {}'.format("\n  ".join(seq["n"] for seq in empty)))

    # remove empty sequences
    sequences = [seq for seq in sequences if seq["s"]]
    if not sequences:
        raise ValueError("No sequences found for exporting")

    # avoid repeated names
    sequences.sort(key=operator.itemgetter("n"))
    prev_name = None
    repeat_no = 0
    for seq in sequences:
        if seq["n"] == prev_name:
            repeat_no += 1
            seq["n"] += "__" + str(repeat_no)
        else:
            prev_name = seq["n"]
            repeat_no = 0

    sorted_by = "name"
    if recent is not None and len(sequences) > recent:
        sequences.sort(key=operator.itemgetter("d"))
        sorted_by = "date"
        sequences = sequences[len(sequences) - recent:]

    if random is not None and len(sequences) > random:
        module_logger.info("choosing {} from {} at random".format(random, len(sequences)))
        random_m.seed()
        chosen = set()
        while len(chosen) < random:
            chosen.add(random_m.randrange(len(sequences)))
        sequences = [sequences[num] for num in chosen]

    if sort_by and sort_by != sorted_by:
        if sort_by == "date":
            sequences.sort(key=operator.itemgetter("d"))
        elif args.sort_by == "name":
            sequences.sort(key=operator.itemgetter("n"))
        else:
            raise ValueError("Unrecognized sort_by argument")

    module_logger.info('{} sequences to export'.format(len(sequences)))
    module_logger.info('Oldest: {}  Newest: {}'.format(sequences[0]["d"], sequences[-1]["d"]))

    # base seq is always the first one in the file, regardless of sorting, to ease specifying the outgroup for GARLI
    if base_seq:
        base_seqs = [get_sequence(make_entry(e), left_part_size) for e in seqdb.iter_seq().filter_name_regex(base_seq)]
        if len(base_seqs) != 1:
            raise ValueError("{} base sequences selected:\n{}".format(len(base_seqs), "\n".join(repr(s["n"]) for s in base_seqs)))
        module_logger.info('base_seq: {}'.format(base_seqs[0]["n"]))
        base_seq_present = seq_present(seqs=sequences, e1=base_seqs[0])
        if base_seq_present:
            # base seq is already there, remove it
            del sequences[base_seq_present]
        sequences = base_seqs + sequences

    if include_seq:
        include_seqs = [get_sequence(make_entry(e1), left_part_size) for seq in include_seq for e1 in seqdb.iter_seq().filter_name_regex(seq)]
        module_logger.info('include_seqs:\n  {}'.format("\n  ".join(ss["n"] for ss in include_seqs)))
        include_seqs = list(filter(lambda e1: not seq_present(seqs=sequences, e1=e1), include_seqs)) # remove already present sequences
        sequences.extend(include_seqs)

    # convert to most common length BEFORE excluding by hamming distance threshold
    if truncate_to_most_common_length:
        truncate_to_most_common(sequences, fill="X" if amino_acids else "-")

    if hamming_distance_threshold:
        sequences = list(filter(lambda s: exclude_by_hamming_distance(sequences[0], s, hamming_distance_threshold), sequences))
    if len(sequences) < 1:
        raise ValueError("Too few ({}) sequences found for exporting".format(len(sequences)))

    if hamming_distance_report:
        hamming_distances = sorted(([e["n"], hamming_distance(sequences[0]["s"], e["s"], sequences[0]["n"], e["n"])] for e in sequences[1:]), key=operator.itemgetter(1), reverse=True)
    else:
        hamming_distances = None

    exp = exporter(output=filename, output_format=output_format, encode_name=encode_name, wrap=wrap)
    module_logger.info('Writing {} sequences'.format(len(sequences)))
    for ss in sequences:
        exp.write(name=ss["n"], sequence=ss["s"])
    return {"base_seq": encode(sequences[0]["n"]) if base_seq else None, "filename": filename, "hamming_distances": hamming_distances, "number_of_sequences": len(sequences)}

# ----------------------------------------------------------------------

def most_common_length(sequences):
    len_stat = collections.Counter(len(e["s"]) for e in sequences)
    return len_stat.most_common(1)[0][0]

# ----------------------------------------------------------------------

def truncate_to_most_common(sequences, fill):
    mclen = most_common_length(sequences)
    module_logger.info('Truncating/extending sequences to the most common length: {}'.format(mclen))
    for entry in sequences:
        ls = len(entry["s"])
        if ls > mclen:
            entry["s"] = entry["s"][:mclen]
        elif ls < mclen:
            entry["s"] += fill * (mclen - ls)
    return mclen

# ----------------------------------------------------------------------

def hamming_distance(s1, s2, n1, n2):
    l = min(len(s1), len(s2))
    hd = sum(1 for pos in range(l) if s1[pos] != s2[pos])
    # module_logger.debug('HD: {} {!r} {!r}'.format(hd, n1, n2))
    return hd

# ----------------------------------------------------------------------

def read_from_file(filename):
    """Yields tuple (name, sequence) for each entry in the file"""
    yield from read_from_string(read_text(filename), filename)

# ----------------------------------------------------------------------

def read_from_string(source, filename):
    """Yields tuple (name, sequence) for each entry in the string"""
    sequence = []
    name = None

    for line_no, line in enumerate(source.splitlines(), start=1):
        if not line or line[0] == ';':               # empty or comment line
            pass
        elif line[0] == '>':
            if name or sequence:
                yield (name, _check_sequence("".join(sequence), name, filename, line_no))
            sequence = []
            name = line[1:].strip()
        else:
            if not name:
                raise FastaReaderError('{filename}:{line_no}: sequence without name'.format(filename=filename, line_no=line_no))
            sequence.append(line.replace("/", "-").upper()) # / found in H1pdm sequences
    if name:
        yield (name, _check_sequence("".join(sequence), name, filename, line_no))

# ----------------------------------------------------------------------

sReSequence = re.compile(r"^[A-Za-z\-~:\*\.]+$")

def _check_sequence(sequence, name, filename, line_no):
    if not sequence:
        raise FastaReaderError('{filename}:{line_no}: {name!r} without sequence'.format(name=name, filename=filename, line_no=line_no))
    if not sReSequence.match(sequence):
        raise FastaReaderError('{filename}:{line_no}: invalid sequence read: {sequence}'.format(sequence=sequence, filename=filename, line_no=line_no))
    return sequence

# ----------------------------------------------------------------------

sReReassortant = re.compile(r"^(.+/\d+)\s*((?:X|BX|NYMC(?:[\s\-]X)|VI|NIB(?:SC)?|RESVIR|IVR|RG(?:\s*ID)?|(?:IBCDC|UW)RG|PR8-IDCDC-RG)[-\s]?[0-9A-Z]+)\s*$")

def detect_reassortant(entry):
    # module_logger.debug('detect_reassortant {!r} -- {!r}'.format(entry.get("name"), entry.get("reassortant")))
    if entry.get("name") and not entry.get("reassortant"):
        m = sReReassortant.match(entry["name"])
        if m and m.group(2):
            entry["reassortant"] = m.group(2).strip()
            entry["name"] = m.group(1).strip()
    # module_logger.debug('detect_reassortant2 {!r} -- {!r}'.format(entry.get("name"), entry.get("reassortant")))

# ----------------------------------------------------------------------

def read_fasta(fasta_file):
    """Returns list of dict {"name":, "sequence":}"""

    def make_entry(raw_name, sequence):
        entry = {"sequence": sequence, "name": raw_name.upper()}
        detect_reassortant(entry)
        return entry

    r = [make_entry(raw_name, sequence) for raw_name, sequence in read_from_string(read_text(fasta_file), fasta_file)]
    module_logger.debug('{} sequences imported from {}'.format(len(r), fasta_file))
    return r

# ----------------------------------------------------------------------

def read_fasta_with_name_parsing(fasta_file, lab, virus_type, **_):
    """Returns list of dict {"name":, "sequence":, "date":, "lab":}"""
    np = name_parser()

    def make_entry(raw_name, sequence):
        n_entry = np.parse(raw_name, lab=lab)
        if not n_entry:
            raise RuntimeError("Cannot parse name: {!r}".format(raw_name))
        entry = {"sequence": sequence, "lab": lab, "virus_type": virus_type}
        entry.update(n_entry)
        for f in ["name", "passage", "lab_id", "virus_type", "lineage"]:
            if entry.get(f):
                entry[f] = entry[f].upper()
        detect_reassortant(entry)
        return entry

    # module_logger.info(f'reading {fasta_file}')
    r = [make_entry(raw_name, sequence) for raw_name, sequence in read_from_string(read_text(fasta_file), fasta_file)]
    module_logger.info('{} sequences imported from {}'.format(len(r), fasta_file))
    return r

# ----------------------------------------------------------------------

sNameParser = None

def name_parser():
    global sNameParser
    if not sNameParser:
        sNameParser = NameParser()
    return sNameParser

class NameParser:

    def __init__(self):
        self.parsers = (
            (re.compile(r'^(?P<name>[^|]+)\s+\|\s+(?:(?P<year1>\d+)-(?P<month1>\d+)-(?P<day1>\d+)|(?P<year2>\d+)-(?P<month2>\d+)\s+\(day unknown\)|(?P<year3>\d+)\s+\(month and day unknown\))\s+\|\s+(?P<passage>[^\|]*)\s+\|\s+(?P<lab_id>[^\|]*)?\s+\|\s+(?P<lab>[A-Za-z ]+)\s+\|\s+(?P<virus_type>[AB]\s*/\s*H\d+N\d+)\s+\|\s*(?P<lineage>[A-Za-z0-9]+)?\s*$', re.I), self.gisaid), # name | date | passage | lab_id | lab | virus_type | lineage
            (re.compile(r'^(?P<name>[^|]+)\s+\|\s+(?:(?P<year1>\d+)-(?P<month1>\d+)-(?P<day1>\d+)|(?P<year2>\d+)-(?P<month2>\d+)\s+\(day unknown\)|(?P<year3>\d+)\s+\(month and day unknown\))\s+\|\s+(?P<passage>[^\|]*)\s+\|\s+(?P<lab_id>[^\|]*)?\s+\|\s+(?P<lab>[A-Za-z ]+)\s*$', re.I), self.gisaid), # name | date | passage | lab_id | lab
            (re.compile(r'^(?P<name>[^|]+)\s+\|\s+(?:(?P<year1>\d+)-(?P<month1>\d+)-(?P<day1>\d+)|(?P<year2>\d+)-(?P<month2>\d+)\s+\(day unknown\)|(?P<year3>\d+)\s+\(month and day unknown\))\s+\|\s+(?P<passage>[^\|]+)\s+\|\s+(?P<lab_id>[^\|]+)?\s+\|.*$', re.I), self.gisaid), # name | date | passage | lab_id | something-else
            (re.compile(r'^(?P<name>[^|]+)\s+\|\s+(?:(?P<year1>\d+)-(?P<month1>\d+)-(?P<day1>\d+)|(?P<year2>\d+)-(?P<month2>\d+)\s+\(day unknown\)|(?P<year3>\d+)\s+\(month and day unknown\))\s+\|\s+(?P<passage>[^\|]+)?\s*\|\s*(?P<lab_id>.+)?\s*$', re.I), self.gisaid), # name | date | passage? | lab_id
            (re.compile(r'^(?P<name1>EPI\d+)\s+\|\s+(?P<gene>HA|NA)\s+\|\s+(?P<designation>[^\|]+)\s+\|\s+(?P<name>EPI_[A-Z_0-9]+)\s+\|\s*(?P<passage>[^\s]+)?\s*\|\s*(?P<flu_type>.+)?\s*$', re.I), self.gisaid_melb), # name1 | gene | designation | name | passage | flu_type
            #(re.compile(r'^(?P<type>B|H3|H1)(?P<location>[A-Z]+)(?P<isolation_number>\d+)(?P<year>\d\d)[A-Z]*\s', re.I), self.nimr_glued),
            (re.compile(r'^(?P<name>[^\s]+)\s+PileUp\sof', re.I), self.nimr_20090914),
            (re.compile(r'^(?P<designation>[^|]+)\s+\|\s+(?P<passage>[^\s]+)\s*\|\s+(?P<name>(?:EPI|201)\d+)\s*$', re.I), self.cdc_20100913), # name | passage | fasta_id
            (re.compile(r'^(?P<name>[^|]+)\s+\|\s+(?P<passage>[^\s]+)\s+\|\s*(?P<lab_id>.+)?\s*$', re.I), self.gisaid_without_date), # name | passage | lab_id
            (re.compile(r'^(?P<name>[^\s]+)\s\s+(?P<date>[\d/]+)\s\s+(?P<passage>[^\s]+)\s*$', re.I), self.melb_20100823), # name  date  passage
            (re.compile(r'^(?P<name>\d+S\d+)\s+\"Contig\s+\d+\"\s+\(\d+,\d+\)$', re.I), self.melb_20110921),
            (re.compile(r'^(?P<name>[^_]+/\d\d\d\d)[\s_]+[^/]*(?:[\s_]?:(?P<passage>[EC]\d*))[^/]*$', re.I), self.name_passage), # CNIC
            (re.compile(r'^(?P<name>[^_]+/\d\d\d\d)[\s_]+(?:Jan|Feb|Mar|Apr|may|Jun|Jul|Aug|Sep|Oct|Nov|Dec)?[^/]*$', re.I), self.name_only), # CNIC
            (re.compile(r'^(?P<name>[^\s]+_4)$', re.I), self.name_only),   # CDC
            (re.compile(r'^(?P<name>[^_]+)[\s_]+[^/]*$', re.I), self.name_only), # CNIC
            (re.compile(r'^(?P<name>.+/\d\d\d\d)[\s\-]*(?P<passage>[\.\dA-Z]+)?$', re.I), self.name_passage), # NIID
            (re.compile(r'.'), self.simple),
            )

    def parse(self, raw_name, lab=None):
        for rex, func_name in self.parsers:
            m = rex.match(raw_name)
            if m:
                # module_logger.debug('NameParser {}'.format(func_name))
                return {k: v for k, v in func_name(raw_name, m, lab=lab).items() if v}
        return None

    def simple(self, raw_name, m, **_):
        return {'name': raw_name}

    # def nimr_glued(self, raw_name, m=None):
    #     def convert_type(t):
    #         if re.match(r'^H\d\d?$', t):
    #             t = 'A'
    #         return t
    #     return {'name': '/'.join((convert_type(m.group('type')), m.group('location'), m.group('isolation_number'), m.group('year')))}

    def nimr_20090914(self, raw_name, m, **_):
        return {'name': m.group('name').upper()}

    mReCdcId = re.compile(r'^\s*(?P<cdcid>\d{8,10})(?:_\d{6}_v\d(?:_\d)?|_\d|\s+.*)?$')

    def gisaid(self, raw_name, m, with_date=True, lab=None, **_):
        groups = m.groupdict()
        # module_logger.debug('gisaid with_date:{} {!r} --> {}'.format(with_date, raw_name, groups))
        year = (with_date and (groups.get('year1') or groups.get('year2') or groups.get('year3'))) or None
        try:
            lab = self._fix_gisaid_lab(m.group('lab'))   # do NOT use groups here!
        except IndexError:
            pass
        lab_id = groups.get('lab_id')
        # module_logger.debug('{} lab_id {} {}'.format(lab, lab_id, self.mReCdcId.match(lab_id)))
        if lab == "CDC" and lab_id is not None:
            m_cdcid = self.mReCdcId.match(lab_id)
            if m_cdcid:
                lab_id = m_cdcid.group("cdcid")
            else:
                if lab_id:
                    module_logger.warning('Not a cdcid: {}'.format(lab_id))
                lab_id = None
        return {
            'name': self._fix_gisaid_name(groups.get('name')),
            'date': year and '-'.join((year, groups.get('month1') or groups.get('month2') or '01', groups.get('day1') or '01')),
            'passage': groups.get('passage'),
            'lab_id': lab_id,
            'lab': lab,
            'virus_type': self._fix_gisaid_virus_type(groups.get("virus_type")),
            'lineage': self._fix_lineage(groups.get("lineage")),
            }

    def _fix_gisaid_lab(self, lab):
        return (lab
                .replace("Centers for Disease Control and Prevention", "CDC")
                .replace("Crick Worldwide Influenza Centre", "NIMR")
                .replace("National Institute for Medical Research", "NIMR")
                .replace("WHO Collaborating Centre for Reference and Research on Influenza", "MELB")
                .replace("National Institute of Infectious Diseases (NIID)", "NIID")
                .replace("National Institute of Infectious Diseases", "NIID")
                .replace("Erasmus Medical Center", "EMC")
                .replace("WHO Chinese National Influenza Center", "CNIC")
                )

    def _fix_gisaid_virus_type(self, virus_type):
        if virus_type:
            m = re.match(r"^\s*([AB])\s*/\s*(H\d+N\d+)\s*$", virus_type)
            if m:
                if m.group(1) == "B":
                    virus_type = "B"
                else:
                    virus_type = "A(" + m.group(2) + ")"
                pass
            else:
                raise ValueError("Unrecognized gisaid flu virus_type: " + virus_type)
        return virus_type

    sReReassortantName = re.compile(r"^((?:NYMC|[IB]VR|B?X-|NIB|RG)[^\(]*)\s*\((?:HY\s+)?([^\)]+)\)?$")

    def _fix_gisaid_name(self, name):
        # module_logger.debug('gisaid name: {!r}'.format(name))
        fixed_name = (name
                 .upper()
                 .replace("ITALY-FVG/", "ITALY/FVG-")   # NIMR uses A/ITALY/FVG-1/2015 in HI
                 .replace("/LYON/CHU/", "/LYON CHU/")
                 .replace("/CAMEROON1", "/CAMEROON/1")
                 .replace("//", "/")
                 )
        m = self.sReReassortantName.match(fixed_name)
        if m:
            fixed_name = m.group(2) + " " + m.group(1)
            # module_logger.debug('fixed reassortant {!r} -> {!r}'.format(name, fixed_name))
        return fixed_name

    def _fix_lineage(self, lineage):
        if lineage:
            lineage = lineage.upper()
            if lineage == "PDM09":
                lineage = "2009PDM"
        return lineage

    def gisaid_without_date(self, raw_name, m, **kw):
        return self.gisaid(raw_name=raw_name, m=m, with_date=False, **kw)

    def gisaid_melb(self, raw_name, m, **_):
        return {'name': m.group('name'), 'gene': m.group('gene')}

    def melb_20100823(self, raw_name, m, **_):
        return {'name': m.group('name').upper(), 'date': m.group('date'), 'passage': m.group('passage')}

    def melb_20110921(self, raw_name, m, **_):
        return {'name': m.group('name').upper()}

    def name_only(self, raw_name, m, **_):
        return {'name': m.group('name').upper()}

    def cdc_20100913(self, raw_name, m, **_):
        return {'name': m.group('name').upper(), 'passage': m.group('passage')}

    def name_passage(self, raw_name, m, **_):
        return {'name': m.group('name').upper(), 'passage': m.group('passage')}

# ----------------------------------------------------------------------

# def generate_one(name, sequence, encode_name, split=True):
#     return ">{}\n{}\n".format((encode(name) if encode_name else name).strip(), (sequence_split(sequence) if split else sequence).strip())

# ----------------------------------------------------------------------

# def sequence_split(sequence, chunk_len=75, separator="\n"):
#     if chunk_len and chunk_len > 0:
#         r = separator.join(sequence[i:i+chunk_len] for i in range(0, len(sequence), chunk_len))
#     else:
#         r = sequence
#     return r

# ----------------------------------------------------------------------

def exporter(output, output_format, encode_name, wrap):
    if output_format == "fasta":
        r = FastaExporter(output=output, encode_name=encode_name, wrap=wrap)
    elif output_format == "phylip":
        r = PhylipExporter(output=output, encode_name=encode_name, wrap=wrap)
    elif output_format == "bioseq": # https://github.com/mnbram/bioseq-mode
        r = BioseqExporter(output=output, encode_name=encode_name, wrap=wrap)
    else:
        raise ValueError("Unrecognized output_format: {}".format(output_format))
    return r

# ----------------------------------------------------------------------

class ExporterBase:

    def __init__(self, output, encode_name, wrap):
        self.output_file = output
        self.output = ""
        module_logger.info('Writing {}'.format(self.output_file))
        self.encode_name = encode_name
        self.wrap = wrap

    def __del__(self):
        if self.output_file == "-":
            print(self.output)
        else:
            write_binary(self.output_file, self.output.encode("utf-8"))

    def sequence_split(self, sequence, chunk_len=60, separator="\n"):
        if chunk_len and chunk_len > 0:
            r = separator.join(sequence[i:i+chunk_len] for i in range(0, len(sequence), chunk_len))
        else:
            r = sequence
        return r

    def make_name(self, name):
        if self.encode_name:
            name = encode(name)
        return name.strip()

class FastaExporter (ExporterBase):

    def write(self, name, sequence):
        self.output += ">{}\n{}\n".format(self.make_name(name), (self.sequence_split(sequence) if self.wrap else sequence).strip())

class PhylipExporter (ExporterBase):

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.names = []
        self.sequences = []

    def __del__(self):
        self.do_write()
        super().__del__()

    def write(self, name, sequence):
        self.names.append(self.make_name(name))
        self.sequences.append(sequence)

    def do_write(self):
        max_s_len = max(len(s) for s in self.sequences)
        max_n_len = max(len(n) for n in self.names)
        self.output += "{} {}\n".format(len(self.sequences), max_s_len)
        if self.wrap:
            raise NotImplementedError("phylip with wrapping")   # http://www.molecularevolution.org/resources/fileformats/phylip_dna
        else:
            for no, (n, s) in enumerate(zip(self.names, self.sequences)):
                self.output += "{:<{}s}  {}{}\n".format(self.names[no], max_n_len, s, "-" * (max_s_len - len(s)))

class BioseqExporter (PhylipExporter):

    def do_write(self):
        self.output = "\n".join(["-*- bioseq -*-\n",
                                 "\n".join(f"{no+1:03d}>{name}" for no, name in enumerate(self.names)),
                                 "\n".join(f"{no+1:03d}#{seq}" for no, seq in enumerate(self.sequences))
                                 ])

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
