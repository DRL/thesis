#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage:
    genome_distributions.py                     -c <FILE> [-n <INT>] [-o <STR>] [-f <STR>]
                                                [-h|--help]

    Options:
        -h, --help                              show this
        -c, --config <FILE>                     config file with ID,BED-file
        -n, --n_threshold <INT>                 Threshold of N's to split scaffolds into contigs [default: 10]
        -o, --outprefix <STRING>                Output prefix
        -f, --format <STR>                      Plot format [default: png]

"""
from __future__ import division
import re
import sys
from docopt import docopt
import os
import itertools
import matplotlib as mat
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
mat.use("agg")
import numpy as np
import matplotlib.pyplot as plt
#mat.rcParams['font.size'] = 14

#plt.style.use('ggplot')

def read_file(infile):
    if not infile or not os.path.exists(infile):
        sys.exit("[ERROR] - File '%s' does not exist." % (infile))
    print "[+] Parsing %s ..." % (infile)
    with open(infile) as fh:
        for line in fh:
            line = line.replace(r'\r', '\n')
            if not line.startswith("#"):
                yield line.rstrip("\n")


def write_file(out_f, outprefix, header, strings):
    if outprefix:
        if outprefix.endswith("/"):
            if not os.path.exists(outprefix):
                os.mkdir(outprefix)
            out_f = "%s" % os.path.join(outprefix, out_f)
        else:
            out_f = "%s.%s" % (outprefix, out_f)
    print "[+] \t Writing file %s ..." % (out_f)
    with open(out_f, 'w') as out_fh:
        out_fh.write("%s\n" % (header))
        out_fh.write("%s\n" % "\n".join(strings))

class PlotObj():
    def __init__(self, label):
        self.label = label
        self.x = []
        self.y = []

class GenomeObj():
    def __init__(self, genome_id):
        self.genome_id = genome_id
        self.scaffold_lengths = []
        self.contig_lengths = []
        self.N_lengths = []

    def add_contig_length(self, length):
        if length:
            self.contig_lengths.append(length)

    def add_scaffold_length(self, length):
        if length:
            self.scaffold_lengths.append(length)

    def add_N_length(self, length):
        if length:
            self.N_lengths.append(length)


class MainObj():
    def __init__(self, args):
        self.config_f = args['--config']
        self.n_threshold = int(args['--n_threshold'])
        self.outprefix = args['--outprefix']
        self.format = args['--format']
        self.bed_f_by_genome_id = {}
        self.order = []
        self.genomeObj_by_genome_id = {}
        self.parse_config_f()
        self.parse_bed_fs()
        self.plot_cummulative_length('scaffolds')
        self.plot_cummulative_length('contigs')
        self.plot_cummulative_length('Ns')
        self.plot_n_bar()

    def plot_n_bar(self):
        plotObjs = []
        out_f, xlab, ylab = None, None, None
        xlab = "Length of runs of N's (in b)"
        ylab = "Count"
        out_f = 'n_bins'
        for genome_id in self.order:
            plotObj = PlotObj(genome_id)
            genomeObj = self.genomeObj_by_genome_id[genome_id]
            plotObj.x = [x for x in genomeObj.N_lengths]
            #for k, g in itertools.groupby(sorted(genomeObj.N_lengths), lambda x: x // 50 * 1):
            #    plotObj.x.append(k)
            #    plotObj.y.append(len(list(g)))
            plotObjs.append(plotObj)
        self.plot(plotObjs, 'n_bar', out_f, xlab, ylab)

    def plot(self, plotObjs, plot_type, out_f, xlab, ylab):
        out_f = "%s.%s" % (out_f, self.format)
        print "[+] Plotting \"%s\" ..." % (out_f)
        if plot_type == 'cum_length':
            f, ax = plt.subplots(figsize=(5,5), dpi=500)
            ax.set_facecolor('white')
            for plotObj in plotObjs:
                ax.plot(plotObj.x, plotObj.y, '--', alpha=0.8, ms=20, label=plotObj.label)
            ax.legend(loc=4, frameon=False)
            ax.set_xlabel(xlab)
            ax.set_ylabel(ylab)
        elif plot_type == 'n_bar':
            bins = 25
            f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True, dpi=500, figsize=(6, 6))
            ax1.hist(plotObjs[0].x, bins, alpha=0.8, color='tab:blue', label=plotObjs[0].label)
            ax1.set_yscale('log')
            ax1.legend(loc='best', frameon=False)
            ax2.hist(plotObjs[1].x, bins, alpha=0.8, color='tab:orange', label=plotObjs[1].label)
            ax2.set_yscale('log')
            ax2.legend(loc='best', frameon=False)
            ax3.hist(plotObjs[2].x, bins, alpha=0.8, color='tab:green', label=plotObjs[2].label)
            ax3.set_yscale('log')
            ax3.legend(loc='best', frameon=False)
            ax2.set_ylabel(ylab)
            ax3.set_xlabel(xlab)
        else:
            sys.exit("[X] Wrong plot type %s" % plot_type)
        #ax.grid(True, linewidth=0.5, which="minor", color="lightgrey")
        #ax.grid(True, linewidth=1, which="major", color="lightgrey")
        f.tight_layout()
        f.savefig(out_f, format=self.format)
        plt.close()

    def plot_cummulative_length(self, length_type):
        plotObjs = []
        out_f, xlab, ylab = None, None, None
        if length_type == 'scaffolds':
            xlab = "Scaffold index"
            ylab = "Cummulative length of scaffolds (in Mb)"
            out_f = 'scaffold_cummulative_length'
        elif length_type == 'contigs':
            xlab = "Contig index"
            ylab = "Cummulative length of contigs (in Mb)"
            out_f = 'contig_cummulative_length'
        elif length_type == 'Ns':
            xlab = "N-block index"
            ylab = "Cummulative length of N's (in Mb)"
            out_f = 'n_cummulative_length'
        else:
            sys.exit("[X] Wrong length type %s" % length_type)
        for genome_id in self.order:
            plotObj = PlotObj(genome_id)
            genomeObj = self.genomeObj_by_genome_id[genome_id]
            cum_length = 0
            if length_type == 'scaffolds':
                for idx, length in enumerate(sorted(genomeObj.scaffold_lengths, reverse=True)):
                    plotObj.x.append(idx)
                    cum_length += length
                    plotObj.y.append(cum_length/1e6)
                print "Scaffold length %s: %s Mb (%s b)" % (genome_id, cum_length/1e6, cum_length)
            elif length_type == 'contigs':
                for idx, length in enumerate(sorted(genomeObj.contig_lengths, reverse=True)):
                    plotObj.x.append(idx)
                    cum_length += length
                    plotObj.y.append(cum_length/1e6)
                print "Contig length %s: %s Mb (%s b)" % (genome_id, cum_length/1e6, cum_length)
            elif length_type == 'Ns':
                for idx, length in enumerate(sorted(genomeObj.N_lengths, reverse=True)):
                    plotObj.x.append(idx)
                    cum_length += length
                    plotObj.y.append(cum_length/1e6)
                print "N-length %s: %s Mb (%s b)" % (genome_id, cum_length/1e6, cum_length)
            else:
                sys.exit("[X] Wrong length type %s" % length_type)
            plotObjs.append(plotObj)
        self.plot(plotObjs, 'cum_length', out_f, xlab, ylab)

    def parse_config_f(self):
        for line in read_file(self.config_f):
            col = line.split(",")
            genome_id = col[0]
            self.order.append(genome_id)
            bed_f = col[1]
            self.bed_f_by_genome_id[genome_id] = bed_f

    def parse_bed_fs(self):
        for genome_id, bed_f in self.bed_f_by_genome_id.items():
            genomeObj = GenomeObj(genome_id)
            last_scaffold_id = None
            contig_length = 0
            scaffold_length = 0
            for line in read_file(bed_f):
                col = line.split()
                scaffold_id = col[0]
                nuc_type = col[3]
                nuc_length = int(col[4])
                # All nuc_types
                if not last_scaffold_id or last_scaffold_id == scaffold_id:  # same scaffold
                    scaffold_length += nuc_length
                else:
                    genomeObj.add_scaffold_length(scaffold_length)
                    scaffold_length = nuc_length
                # Only nuc_types == N
                if nuc_type == 'N':
                    genomeObj.add_N_length(nuc_length)
                    if nuc_length >= self.n_threshold:
                        genomeObj.add_contig_length(contig_length)
                        contig_length = 0
                else:
                    contig_length += nuc_length
                last_scaffold_id = scaffold_id
            if contig_length:
                genomeObj.add_contig_length(contig_length)
            if scaffold_length:
                genomeObj.add_scaffold_length(scaffold_length)
            self.genomeObj_by_genome_id[genome_id] = genomeObj
            # print "Ns: ", len(genomeObj.N_lengths)
            # print "Scaffolds: ", len(genomeObj.scaffold_lengths)
            # print "Contigs: ", len(genomeObj.contig_lengths)



if __name__ == "__main__":
    __version__ = 0.2
    args = docopt(__doc__)
    print "[+] Start ..."
    mainObj = MainObj(args)
