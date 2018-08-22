#! /usr/bin/env python3
# -*- Python -*-

"""
Generates regex to match all the strings.
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
    all_len = min(len(arg)for arg in args.input)
    aas = [set() for pos in range(all_len)]
    for seq in args.input:
        for pos in range(all_len):
            aas[pos].add(seq[pos])
    result = ""
    for pos in range(all_len):
        if len(aas[pos]) == 1:
            result += "".join(list(aas[pos]))
        else:
            result += "[" + "".join(list(aas[pos])) + "]"
    print(aas)
    print(result)

# ----------------------------------------------------------------------

with timeit(sys.argv[0]):
    try:
        import argparse
        # print(sys.argv)
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('input', nargs="*", help='Sequences to compare.')
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
