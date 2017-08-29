#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage:
    generate_upsetr_expression.py               -f <FILE> [-s <INT>] [-h|--help]

    Options:
        -h, --help                              show this
        -f, --file <FILE>                       binary file with read counts (first row is sample names)
        -s, --sets <INT>                        number of sets to write [default: 100]
"""

from __future__ import division
from docopt import docopt
import os
import sys

def read_file(infile):
    if not infile or not os.path.exists(infile):
        sys.exit("[X] - File '%s' does not exist." % (infile))
    #print "[+] Parsing file %s ..." % (infile)
    with open(infile) as fh:
        for line in fh:
            yield line

class Main():
    def __init__(self, args):
        self.file = args['--file']
        self.set_count = int(args['--sets'])
        self.translation = {'ERR114517.Lindley': 'Lindley',
                            'ERR123952.Bedale': 'Bedale',
                            'ERR123953.Luffness': 'Luffness',
                            'ERR123954.Newton': 'Newton',
                            'ERR123955.Pa1': 'Pa1',
                            'ERR123956.P5A': 'P5A',
                            'ERR123957.P4A': 'P4A',
                            'Gp12-17.WGA.Nextera': 'Pa1_A.WGA_N',
                            'Gp12-17.WGA.Nextera_XT': 'Pa1_A.WGA_NX',
                            'Gp12-17.WGS.Nextera_XT': 'Pa1_A.WGS_NX',
                            'Gp24.WGA.Nextera.vs': 'Pa1_B.WGA_N',
                            'Gp24.WGA.Nextera_XT.vs': 'Pa1_B.WGA_NX',
                            'Gp24.WGS.Nextera_XT.vs': 'Pa1_B.WGS_NX',
                            'Gp4-8.WGA.Nextera': 'Luffness.WGA_N',
                            'Gp4-8.WGA.Nextera_XT': 'Luffness.WGA_NX',
                            'Gp4-8.WGS.Nextera_XT': 'Luffness.WGS_NX'}
        self.counter = {}
        self.parse_file()
        self.generate_output()

    def generate_output(self):
        output = ['expressionInput <- c(']
        output_sets = []
        processed_sets = 0
        for set_name, count in sorted(self.counter.items(), key=lambda k: k[1], reverse=True):
            if processed_sets < self.set_count:
                list_name = list(set_name)
                if len(list_name) == 1:
                    output_sets.append("%s = %s" % (list_name[0], count))
                else:
                    output_sets.append("`%s` = %s" % ("&".join(list_name), count))
            processed_sets += 1
        output.append(", ".join(output_sets))
        output.append(")")
        expression = "".join(output)
        print expression

    def parse_file(self):
        header = []
        for line in read_file(self.file):
            col = line.split()[3:]
            if not header:
                header = [self.translation[x] for x in col]
            else:
                set_name = frozenset([header[idx] for idx, value in enumerate(col) if value != 'NA'])
                if not len(set_name) == 0:
                    if set_name not in self.counter:
                        self.counter[set_name] = 0
                    self.counter[set_name] += 1


if __name__ == "__main__":
    '''
    expressionInput <- c(one = 2, two = 1, three = 2, `one&two` = 1, `one&three` = 4,
    `two&three` = 1, `one&two&three` = 2)
    '''
    __version__ = 0.1
    args = docopt(__doc__)
    main = Main(args)

