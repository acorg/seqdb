# -*- Python -*-
# license
# license.

from seqdb_backend import *
from .fasta import export_from_seqdb
from .update import create, SeqdbUpdater

# ----------------------------------------------------------------------

def open(path_to_seqdb):
    seqdb = Seqdb()
    seqdb.load(filename=str(path_to_seqdb))
    return seqdb

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
