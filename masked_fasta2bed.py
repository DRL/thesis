#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: masked_fasta2bed.py           -f FASTA [-h|--help]

    Options:
        -h --help                       show this
        -f, --fasta FILE                Masekd FASTA file
"""

from __future__ import division
from __future__ import print_function
import sys
from docopt import docopt
import os

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    sys.exit()

def readFasta(infile):
    if not infile or not os.path.exists(infile):
        eprint("[X] - File '%s' does not exist." % (infile))
    with open(infile) as fh:
        header, seqs = '', []
        for l in fh:
            if l[0] == '>':
                if header:
                    yield header, ''.join(seqs)
                header, seqs = l[1:-1].split()[0], []  # Header is split at first whitespace
            else:
                seqs.append(l[:-1])
        yield header, ''.join(seqs)

class Data():
    def __init__(self, infile):
        self.infile = infile
        self.create_bed()

    def create_bed(self):
        region_by_char = {'a': 'agct', 'g': 'agct', 'c': 'agct', 't': 'agct', 'n': 'n',
                          'A': 'AGCT', 'G': 'AGCT', 'C': 'AGCT', 'T': 'AGCT',
                          'R': 'IUPAC', 'Y': 'IUPAC', 'S': 'IUPAC', 'W': 'IUPAC', 'K': 'IUPAC', 'M': 'IUPAC', 'B': 'IUPAC', 'D': 'IUPAC', 'H': 'IUPAC', 'V' : 'IUPAC',
                          'N': 'N'}
        output = []
        count = 0
        for header, seq in readFasta(self.infile):
            count += 1
            start, end, current_region, previous_region = 0, 0, None, None
            for idx, char in enumerate(seq):
                if not previous_region:
                    previous_region = region_by_char[char]
                current_region = region_by_char[char]
                if current_region != previous_region:
                    end = idx
                    output.append("\t".join([header, str(start), str(end), previous_region, str(end-start)]))
                    start = idx
                else:
                    end += 1
                previous_region = current_region
            end = idx+1
            output.append("\t".join([header, str(start), str(end), current_region, str(end-start)]))
        print("\n".join(output))


if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    infile = args['--fasta']
    data = Data(infile)
