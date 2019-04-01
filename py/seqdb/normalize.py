# -*- Python -*-
# license
# license.
# ----------------------------------------------------------------------

# import datetime
import re
import logging; module_logger = logging.getLogger(__name__)

# ----------------------------------------------------------------------

def virus_type_lineage(vt, lin):
    vtn = virus_type(vt)
    if vtn == "B" and lin is None:
        lin = vt[1:]
        if lin and lin[0] == "/":
            lin = vt[2:]
    linn = lineage(lin)
    if vtn == "B" and linn is None:
        raise ValueError("No lineage for " + vt)
    return vtn, linn

# ----------------------------------------------------------------------

sVirusTypes = {
    "B": "B", "BV": "B", "BY": "B", "BVIC": "B", "BYAM": "B", "B/VIC": "B", "B/YAM": "B",
    "H1PDM": "A(H1N1)", "H1SEAS": "A(H1N1)", "H1": "A(H1N1)", "A(H1N1)": "A(H1N1)",
    "H3": "A(H3N2)", "A(H3N2)": "A(H3N2)",
    "ALL": "ALL",
    }

def virus_type(vt):
    if vt is not None:
        try:
            vt = sVirusTypes[vt.upper()]
        except KeyError:
            raise ValueError("Unrecognized virus type {}".format(vt))
    return vt

# ----------------------------------------------------------------------

sLabs = {
    "CDC": "CDC",
    "CNIC": "CNIC",
    "MELB": "MELB", "VIDRL": "MELB",
    "NIID": "NIID",
    "NIMR": "NIMR", "CRICK": "NIMR",
    }

def lab(lab):
    if lab is not None:
        try:
            lab = sLabs[lab.strip().upper()]
        except KeyError:
            raise ValueError("Unrecognized lab: {!r}".format(lab))
    return lab

# ----------------------------------------------------------------------

sLineages = {
    "VICTORIA": "VICTORIA", "VIC": "VICTORIA", "V": "VICTORIA",
    "YAMAGATA": "YAMAGATA", "YAM": "YAMAGATA", "Y": "YAMAGATA",
    }

def lineage(lineage):
    if lineage is not None:
        try:
            lineage = sLineages[lineage.strip().upper()]
        except KeyError:
            raise ValueError("Unrecognized lineage: {!r}".format(lineage))
    return lineage

# ----------------------------------------------------------------------

sReDate = re.compile(r"(?P<year>\d\d\d\d)-?(?:(?P<month>\d\d)-?(?:(?P<day>\d\d))?)?$")

def date(date):
    if not date:
        date = ""
    else:
        m = sReDate.match(date)
        if m:
            year = m.group("year")
            month = m.group("month") or "01"
            day = m.group("day") or "01"
            date = "-".join((year, month, day))
        else:
            raise ValueError("Invalid date: {!r}, YYYY[MM[DD]] expected")
    return date

# ----------------------------------------------------------------------

sRePassages = [
    [re.compile(r"^\s*(CS|ORIGINAL)(?:\([^\)]*\))?\s*$"), r"OR"],
    [re.compile(r"(?<![A-Z_])C(X|\d+)"), r"MDCK\1"],
    [re.compile(r"(?<![A-Z_])S(X|\d+)"), r"SIAT\1"],
    [re.compile(r"(?<![A-Z_])M(X|\d+)"), r"MK\1"],
    [re.compile(r"(?<![A-Z_])(E|SIAT|MK|MDCK)X"), r"\1?"],
    [re.compile(r"\s*,\s*"), r"/"],
    ]

def passage(passage):
    for re_sub in sRePassages:
        passage = re_sub[0].sub(re_sub[1], passage)
    return passage

# ----------------------------------------------------------------------

sReReassortants = [
    [re.compile(r"^(?:NYMC)?\s*B?X\s*-?\s*(\d+[A-C]?)$"), r"NYMC-\1"],
    ]
def reassortant(reassortant):
    for re_sub in sReReassortants:
        reassortant = re_sub[0].sub(re_sub[1], reassortant)
    return reassortant

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
