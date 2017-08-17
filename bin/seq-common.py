#! /usr/bin/env python3
# -*- Python -*-

"""
Prints number of common AAs in two AA sequences
"""

import sys, os, traceback, pprint
if sys.version_info.major != 3: raise RuntimeError("Run script with python3")
from pathlib import Path
sys.path[:0] = [str(Path(os.environ["ACMACSD_ROOT"]).resolve().joinpath("py"))]
import logging; module_logger = logging.getLogger(__name__)

from acmacs_base import timeit
#import seqdb

# ----------------------------------------------------------------------

def main(args):
    s1, s2 = args.input
    num_common = 0
    for pos in range(min(len(s1), len(s2))):
        if s1[pos] == s2[pos]:
            num_common += 1
    print("common:", num_common)

# ----------------------------------------------------------------------

with timeit(sys.argv[0]):
    try:
        import argparse
        # print(sys.argv)
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('input', nargs=2, help='Sequences to compare.')
        parser.add_argument('-d', '--debug', action='store_const', dest='loglevel', const=logging.DEBUG, default=logging.INFO, help='Enable debugging output.')
        args = parser.parse_args()
        logging.basicConfig(level=args.loglevel, format="%(levelname)s %(asctime)s: %(message)s")
        exit_code = main(args)
    except Exception as err:
        logging.error('{}\n{}'.format(err, traceback.format_exc()))
        exit_code = 1
exit(exit_code)

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
