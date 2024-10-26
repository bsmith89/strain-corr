#!/usr/bin/env python3
"""Normalize a bibliography

Consumes all entries in all input files
and sorts them by citation key.

Should normalize entries via the biblib package:
<https://github.com/aclements/biblib>.

Prints the sorted entries to stdout.

USAGE:
    sort_bib.py <input_file1.bib> [<input_file2.bib> [...]] > <output.bib>

"""


import sys
import itertools
from biblib import bib

def entries(path):
    with open(path) as f:
        db = bib.Parser().parse(f, log_fp=sys.stderr).get_entries()
    for entry in db.values():
        yield entry

def get_key(entry):
    return entry.key

if __name__ == '__main__':
    entries_chain = itertools.chain(*[entries(path) for path in sys.argv[1:]])
    sorted_entries = sorted(entries_chain, key=lambda e: get_key(e))
    for entry in sorted_entries:
        print(entry.to_bib())
