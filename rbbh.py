#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Usage:
        rbbh_proteins.py                         -1 <DIR> -2 <DIR> [-o <STRING>]

Options:
        -h --help                                show this
        -1, --input_dir1 <DIR>                   Input directory of proteome 1
        -2, --input_dir2 <DIR>                   Input directory of proteome 2
        -o, --outprefix <STRING>                 Output prefix

"""

from __future__ import division
from docopt import docopt
from collections import Counter
import sys
import os
import matplotlib as mat
mat.use("agg")
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')

"""[summary]
Analysis of RBBHs of proteins

[description]
- Reads FASTAs (*.faa) of proteins to determine counts, lengths
- Read BLAST hits (*.out)
- Calculates RBBH
- Write RBBH
- Output metrics:
    - for how many sequences in each proteome RBBHs could be found
"""


def read_file(infile):
    if not infile or not os.path.exists(infile):
        sys.exit("[X] - File '%s' does not exist." % (infile))
    print "[+] Parsing file %s ..." % (infile)
    with open(infile) as fh:
        for line in fh:
            yield line


def write_file(out_f, outprefix, header, string):
    if outprefix:
        if outprefix.endswith("/"):
            if not os.path.exists(outprefix):
                os.mkdir(outprefix)
            out_f = "%s" % os.path.join(outprefix, out_f)
        else:
            out_f = "%s.%s" % (outprefix, out_f)
    print "[+] \t Writing file %s ..." % (out_f)
    with open(out_f, 'w') as out_fh:
        if header:
            out_fh.write("%s\n" % (header))
        out_fh.write("%s\n" % "\n".join(string))

class ProteinObj():
    def __init__(self, protein_id, protein_seq):
        self.protein_id = protein_id
        self.seq = protein_seq.rstrip("\n")
        self.length = len(protein_seq)
        self.rbbh_protein_id = None

class ProteomeObj():
    def __init__(self, proteome_id, protein_f, protein_blast_f):
        self.proteome_id = proteome_id
        self.protein_f = protein_f
        self.blast_f = protein_blast_f
        self.protein_length = 0
        self.protein_count = 0
        self.protein_ids = []
        self.proteinObjs_by_protein_id = {}
        self.hitObjs = []
        self.hit_count = 0
        self.parse_fasta_f()
        self.parse_blast_f()

    def parse_fasta_f(self):
        header, seqs = '', []
        for line in read_file(self.protein_f):
            if line[0] == '>':
                if header:
                    proteinObj = ProteinObj(header, ''.join(seqs))
                    self.add_proteinObj(proteinObj)
                header, seqs = line[1:-1].split()[0], []  # Header is split at first whitespace
            else:
                seqs.append(line)
        proteinObj = ProteinObj(header, ''.join(seqs))
        self.add_proteinObj(proteinObj)
        print "[+] \t %s proteins parsed ..." % (self.protein_count)

    def add_proteinObj(self, proteinObj):
        self.protein_count += 1
        self.protein_length += proteinObj.length
        self.protein_ids.append(proteinObj.protein_id)
        self.proteinObjs_by_protein_id[proteinObj.protein_id] = proteinObj

    def parse_blast_f(self):
        last_hitObj_pair = None
        for line in read_file(self.blast_f):
            col = [x.strip() for x in line.split("\t")]  # cleaning leading and trailing whitespaces, and turning it into a list ...
            if not col[0] in self.proteinObjs_by_protein_id:
                sys.exit("[X] qseqid '%s' in BLAST file %s is not part of protein FASTA file %s." % (col[0], self.blast_f, self.protein_f))
            try:
                hitObj = HitObj(col)
            except TypeError:
                sys.exit("[X] BLAST outfmt should bed '6 std qlen slen qcovs qcovhsp'")
            if not hitObj.pair == last_hitObj_pair:
                self.hitObjs.append(hitObj)
                self.hit_count += 1
                last_hitObj_pair = hitObj.pair
        print "[+] \t %s hits parsed ..." % (self.hit_count)


class HitObj():
    def __init__(self, col):
        self.qseqid = col[0]
        self.sseqid = col[1]
        self.pair = frozenset([col[0], col[1]])
        self.pident = float(col[2])
        self.evalue = float(col[10])
        self.bitscore = float(col[11])
        self.qlen = int(col[12])
        self.slen = int(col[13])
        self.qcov = int(col[14])

class FastaObj():
    def __init__(self, proteinObj1, proteinObj2):
        self.header1 = proteinObj1.protein_id
        self.seq1 = proteinObj1.seq
        self.header2 = proteinObj2.protein_id
        self.seq2 = proteinObj2.seq

    def get_fasta_list(self):
        return [">%s" % self.header1, self.seq1, ">%s" % self.header2, self.seq2]

class AnalysisCollection():
    def __init__(self, proteomeObjs):
        self.proteomeObjs_by_proteome_id = {proteomeObj.proteome_id: proteomeObj for proteomeObj in proteomeObjs}
        self.proteome_ids = sorted([proteomeObj.proteome_id for proteomeObj in proteomeObjs])
        self.analyse_rbbhs()

    def analyse_rbbhs(self):
        print "[+] Calculate RBBHs (reciprocal best BLAST hits) ..."
        hitObjs = []
        for proteome_id in self.proteome_ids:
            for hitObj in self.proteomeObjs_by_proteome_id[proteome_id].hitObjs:
                hitObjs.append(hitObj)
        seen_pairs = set()
        seen_proteins = set()
        rbbh_protein_id_pairs = []
        length_difference_by_pair = {}
        bitscores_by_pair = {}
        for hitObj in sorted(hitObjs, key=lambda x: x.bitscore, reverse=True):
            if hitObj.evalue <= EVALUE and hitObj.qcov >= QCOV:
                if hitObj.pair in seen_pairs:
                    rbbh_protein_id_pairs.append(hitObj.pair)
                    bitscores_by_pair[hitObj.pair].append(hitObj.bitscore)
                    length_difference_by_pair[hitObj.pair] = abs(hitObj.qlen - hitObj.slen)
                else:
                    if hitObj.pair.intersection(seen_proteins):
                        seen_proteins.add(hitObj.qseqid)
                        seen_proteins.add(hitObj.sseqid)
                    else:
                        bitscores_by_pair[hitObj.pair] = [hitObj.bitscore]
                        seen_pairs.add(hitObj.pair)
        header = "#%s\tlength_diff\tmean_bitscore" % "\t".join(self.proteome_ids)
        body = []
        rbbh_fastaObjs = []
        for rbbh_protein_id_pair in rbbh_protein_id_pairs:
            protein_a, protein_b = rbbh_protein_id_pair
            if protein_a in self.proteomeObjs_by_proteome_id[self.proteome_ids[0]].proteinObjs_by_protein_id:
                proteinObj1 = self.proteomeObjs_by_proteome_id[self.proteome_ids[0]].proteinObjs_by_protein_id[protein_a]
                proteinObj2 = self.proteomeObjs_by_proteome_id[self.proteome_ids[1]].proteinObjs_by_protein_id[protein_b]
                rbbh_fastaObjs.append(FastaObj(proteinObj1, proteinObj2))
                body.append("%s\t%s\t%s\t%s" % (protein_a, protein_b, length_difference_by_pair[rbbh_protein_id_pair], sum(bitscores_by_pair[rbbh_protein_id_pair]) / 2))
            else:
                proteinObj2 = self.proteomeObjs_by_proteome_id[self.proteome_ids[1]].proteinObjs_by_protein_id[protein_a]
                proteinObj1 = self.proteomeObjs_by_proteome_id[self.proteome_ids[0]].proteinObjs_by_protein_id[protein_b]
                rbbh_fastaObjs.append(FastaObj(proteinObj1, proteinObj2))
                body.append("%s\t%s\t%s\t%s" % (protein_b, protein_a, length_difference_by_pair[rbbh_protein_id_pair], sum(bitscores_by_pair[rbbh_protein_id_pair]) / 2))
        rbbh_out_f = "%s.%s.%s.%s.txt" % (".vs.".join(self.proteome_ids), EVALUE, QCOV, 'rbbh')
        write_file(rbbh_out_f, outprefix, header, body)
        print "[+] %s RBBHs found ..." % len(rbbh_protein_id_pairs)
        for proteome_id in self.proteome_ids:
            proteomeObj = self.proteomeObjs_by_proteome_id[proteome_id]
            percentage_with_rbbh = 100 * (len(rbbh_protein_id_pairs) / proteomeObj.protein_count)
            print "[+]\t%s%% of proteins in %s have RBBHs ..." % (percentage_with_rbbh, proteomeObj.protein_f)
        for idx, rbbh_fastaObj in enumerate(rbbh_fastaObjs):
            out_f = "%s.%s.rbbh_%s.%s.%s.faa" % (self.proteome_ids[0], self.proteome_ids[1], idx, EVALUE, QCOV)
            write_file(out_f, outprefix, None, rbbh_fastaObj.get_fasta_list())

def generate_proteomeObjs(input_dirs):
    proteomeObjs = []
    for input_dir in input_dirs:
        proteome_id = input_dir
        paths = {'fasta_f': '', 'blast_f': ''}
        for root, dirs, files in os.walk(input_dir):
            for infile in files:
                if infile.endswith(".faa"):
                    paths['fasta_f'] = os.path.join(root, infile)
                elif infile.endswith(".out"):
                    paths['blast_f'] = os.path.join(root, infile)
                else:
                    pass
        proteomeObj = ProteomeObj(proteome_id, paths['fasta_f'], paths['blast_f'])
        proteomeObjs.append(proteomeObj)
    return proteomeObjs

if __name__ == "__main__":
    __version__ = 0.2

    NUCCOUNT = 2
    EVALUE = 1e-5
    QCOV = 25
    args = docopt(__doc__)
    input_dir_1 = args['--input_dir1']
    input_dir_2 = args['--input_dir2']
    input_dirs = [input_dir_1, input_dir_2]
    outprefix = args['--outprefix']
    proteomeObjs = generate_proteomeObjs(input_dirs)
    analysisCollection = AnalysisCollection(proteomeObjs)


