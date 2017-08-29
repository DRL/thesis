#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: fasta2counts.py           -f FASTA [-b BED] [-w INT] [-h|--help]

    Options:
        -h --help                       show this
        -f, --fasta FILE                FASTA file
        -b, --bed FILE                  BED File
        -w, --window INT                Window size
"""

from __future__ import division
from __future__ import print_function
import sys
from docopt import docopt
from collections import Counter, OrderedDict
import os
from interval import interval, inf, imath

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    sys.exit()

def read_file(infile):
    if not infile or not os.path.exists(infile):
        eprint("[X] - File '%s' does not exist." % (infile))
    with open(infile) as fh:
        for l in fh:
            yield l

def readFasta(infile):
    header, seqs = '', []
    for l in read_file(infile):
        if l[0] == '>':
            if header:
                yield header, ''.join(seqs)
            header, seqs = l[1:-1].split()[0], []  # Header is split at first whitespace
        else:
            seqs.append(l[:-1])
    yield header, ''.join(seqs)

def calculate_metrics(seq):
    seq = seq.upper()
    length = len(seq)
    base_counter = Counter(seq)
    n_count = base_counter['N']
    gc_count = base_counter['G'] + base_counter['C']
    gc_unnorm = gc_count/length
    agct_count = base_counter['G'] + base_counter['C'] + base_counter['A'] + base_counter['T']
    if agct_count:
        n_frac = n_count/agct_count
        gc_norm = gc_count/agct_count
    else:
        n_frac = 'NA'
        gc_norm = 'NA'
    non_agctn_count = length - base_counter['N'] - agct_count
    return str(length), str(n_count), str(n_frac), str(agct_count), str(gc_unnorm), str(gc_norm), str(non_agctn_count)

class Data():
    def __init__(self, args):
        self.infile = args['--fasta']
        if args['--window']:
            self.window = int(args['--window'])
        else:
            self.window = None
        self.gene_count_by_seq_id = {}
        self.gene_intervals_by_seq_id = {}
        if args['--bed']:
            self.bed_f = args['--bed']
            self.get_gene_counts()
        else:
            self.bed_f = None
        self.write_counts()

    def get_gene_counts(self):
        gene_count_by_seq_id = {}
        gene_intervals_by_seq_id = {}
        for l in read_file(self.bed_f):
            col = l.strip("\n").split("\t")
            seq_id = col[0]
            start = int(col[1])
            end = int(col[2])
            inval = interval([start, end])
            if not seq_id in gene_count_by_seq_id:
                gene_count_by_seq_id[seq_id] = 0
                gene_intervals_by_seq_id[seq_id] = []
            gene_count_by_seq_id[seq_id] += 1
            gene_intervals_by_seq_id[seq_id].append(inval)
        self.gene_count_by_seq_id = gene_count_by_seq_id
        self.gene_intervals_by_seq_id = gene_intervals_by_seq_id

    def write_counts(self):
        basename = ('.').join(self.infile.split('.')[:-1])
        metrics_f = "%s.metrics.txt" % basename
        metrics_window_f = "%s.metrics_window_%s.txt" % (basename, self.window)
        output_metrics = []
        output_metrics_window = []
        header_metrics = "\t".join(['header', 'length', 'N', 'N_frac', 'AGCT', 'GC_unnorm', 'GC_norm', 'non_AGCTN', 'gene_count'])
        header_metrics_window = "\t".join(['header', 'length', 'N', 'N_frac', 'AGCT', 'GC_unnorm', 'GC_norm', 'non_AGCTN', 'gene_count'])
        output_metrics.append(header_metrics)
        output_metrics_window.append(header_metrics_window)
        seq_by_header = OrderedDict()
        for header, seq in readFasta(self.infile):
            seq_by_header[header] = seq
        for header, seq in seq_by_header.items():
            seq = seq.upper()
            length, n_count, n_frac, agct_count, gc_unnorm, gc_norm, non_agctn_count = calculate_metrics(seq)
            if self.bed_f:
                gene_count = str(self.gene_count_by_seq_id.get(header, 0))
            else:
                gene_count = "NA"
            output_metrics.append("\t".join([header, length, n_count, n_frac, agct_count, gc_unnorm, gc_norm, non_agctn_count, gene_count]))
            if self.window:
                for i in xrange(0, len(seq), self.window):
                    window_interval = interval([i, i + self.window])
                    window_id = "%s.%s-%s" % (header, i, i + self.window)
                    window_gene_count = 0
                    if header in self.gene_intervals_by_seq_id:
                        for gene_interval in self.gene_intervals_by_seq_id[header]:
                            if gene_interval in window_interval:
                                # print("%s:%s in %s:%s" % (header, gene_interval, header, window_interval))
                                window_gene_count += 1
                    length, n_count, n_frac, agct_count, gc_unnorm, gc_norm, non_agctn_count = calculate_metrics(seq[i:i + self.window])
                    output_metrics_window.append("\t".join([window_id, length, n_count, n_frac, agct_count, gc_unnorm, gc_norm, non_agctn_count, str(window_gene_count)]))
        with open(metrics_f, 'w') as metrics_fh:
            metrics_fh.write("\n".join(output_metrics) + "\n")
        if self.window:
            with open(metrics_window_f, 'w') as metrics_window_fh:
                metrics_window_fh.write("\n".join(output_metrics_window) + "\n")

if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    data = Data(args)
