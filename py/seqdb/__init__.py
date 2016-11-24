from seqdb_backend import *
from .fasta import export_from_seqdb
from .update import SeqdbUpdater

# ----------------------------------------------------------------------

def open(path_to_seqdb):
    seqdb = Seqdb()
    seqdb.load(filename=str(path_to_seqdb))
    return seqdb

# ----------------------------------------------------------------------
