#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Usage:
        parse_rbbh_tables.py                     -f <FILE> --rbbh_gpal_lit <FILE> --rbbh_gros_lit <FILE> --rbbh_gpal_gros <FILE>

Options:
        -h --help                                show this
        -f, --function <FILE>                    list of effectors,function,species
        --rbbh_gpal_lit <FILE>                   GPALL.vs.literature_GPALL.1e-05.25.rbbh.txt
        --rbbh_gros_lit <FILE>                   GROST.vs.literature_GROST.1e-05.25.rbbh.txt
        --rbbh_gpal_gros <FILE>                  GPALL.vs.GROST.1e-05.25.rbbh.txt
"""

from __future__ import division
from docopt import docopt
from collections import Counter
import sys
import os

"""[summary]
Analysis of RBBHs results

[description]
- reads files
- for each PCN
    - write protein_id, orthologue in literature, function, known before, orthologue in other PCN
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
    def __init__(self, protein_id, species):
        self.protein_id = protein_id
        self.effector = False
        self.part_of_lit = False
        self.rbbh_lit_id = None
        self.rbbh_lit_species = None
        self.rbbh_pcn_id = None
        self.rbbh_ids = set()
        self.rbbh_types = set()
        self.species = species
        self.functions = set()

    def __str__(self):
        return "# ID: %s\n%s" % (self.protein_id, '\n'.join("%s: %s" % item for item in vars(self).items()))

    def add_rbbh_lit(self, litObj):
        self.rbbh_lit_id = litObj.protein_id
        self.rbbh_ids.add(litObj.protein_id)
        self.rbbh_lit_species = litObj.species
        self.functions.add(litObj.function)
        self.rbbh_types.add('lit')
        self.effector = True

    def add_rbbh_pcn(self, proteinObj):
        self.rbbh_pcn_id = proteinObj.protein_id
        self.rbbh_types.add('pcn')
        self.rbbh_ids.add(proteinObj.protein_id)
        if proteinObj.effector:
            self.effector = True
        if proteinObj.functions:
            for function in proteinObj.functions:
                self.functions.add(function)



class LitObj():
    def __init__(self, protein_id, species, function):
        self.protein_id = protein_id
        self.species = species
        self.function = function

    def __str__(self):
        return "ID: %s\n\tspecies: %s\n\tfunction: %s" % (self.protein_id, self.species, self.function)

class MainObj():
    def __init__(self, args):
        self.function_f = args['--function']  # convert IDs if gpal/gros
        self.proteinObjs_by_protein_id = {}
        self.litObjs_by_lit_id = {}
        self.protein_ids = []
        self.parse_function_f()
        self.rbbh_gpal_lit_f = args['--rbbh_gpal_lit']
        self.rbbh_gros_lit_f = args['--rbbh_gros_lit']
        self.parse_rbbh_pcn_lit()
        self.rbbh_gpal_gros_f = args['--rbbh_gpal_gros']
        self.parse_rbbh_pcn()
        self.evaluate_rbbhs()

    def evaluate_rbbhs(self):
        effector_output = []
        effector_header = "\t".join(['proteinID', 'species', 'effector', 'known', 'rbbh_literature', 'rbbh_literature_species', 'rbbh_pcn', 'class'])
        for protein_id in self.protein_ids:
            proteinObj = self.proteinObjs_by_protein_id[protein_id]
            if proteinObj.effector == True:
                effector_line = []
                effector_line.append(proteinObj.protein_id)
                effector_line.append(proteinObj.species)
                effector_line.append(proteinObj.effector)
                effector_line.append(proteinObj.part_of_lit)
                effector_line.append(proteinObj.rbbh_lit_id)
                effector_line.append(proteinObj.rbbh_lit_species)
                effector_line.append(proteinObj.rbbh_pcn_id)
                effector_line.append(";".join(proteinObj.functions))
                effector_output.append("\t".join([str(x) for x in effector_line]))
        effector_out_f = "rbbh_effector_table.txt"
        write_file(effector_out_f, None, effector_header, effector_output)

    def add_litObj(self, litObj):
        if not litObj.protein_id in self.litObjs_by_lit_id:
            self.litObjs_by_lit_id[litObj.protein_id] = litObj
        else:
            sys.exit("[X] - ID '%s' already exists." % (litObj.protein_id))

    def parse_rbbh_pcn(self):
        for line in read_file(self.rbbh_gpal_gros_f):
            if not line.startswith("#"):
                col = [x.strip() for x in line.split("\t")]
                gpal_id = col[0]
                proteinObj_1 = ProteinObj(gpal_id, 'Globodera pallida')
                gros_id = ".".join(col[1].split(".")[0:-1])
                proteinObj_2 = ProteinObj(gros_id, 'Globodera rostochiensis')
                self.add_proteinObjs(proteinObj_1, proteinObj_2)

    def parse_rbbh_pcn_lit(self):
        for line in read_file(self.rbbh_gpal_lit_f):
            if not line.startswith("#"):
                col = [x.strip() for x in line.split("\t")]
                protein_id = col[0]
                rbbh_lit_id = col[1]
                species = "Globodera pallida"
                proteinObj = ProteinObj(protein_id, species)
                litObj = self.litObjs_by_lit_id[rbbh_lit_id]
                self.add_proteinObj(proteinObj, litObj)
        for line in read_file(self.rbbh_gros_lit_f):
            if not line.startswith("#"):
                col = [x.strip() for x in line.split("\t")]
                protein_id = ".".join(col[0].split(".")[0:-1])
                rbbh_lit_id = col[1]
                species = "Globodera rostochiensis"
                proteinObj = ProteinObj(protein_id, species)
                litObj = self.litObjs_by_lit_id[rbbh_lit_id]
                self.add_proteinObj(proteinObj, litObj)

    def add_proteinObj(self, proteinObj, litObj):
        # operates only on lit_rbbhs
        if proteinObj.protein_id in self.litObjs_by_lit_id:
            proteinObj.part_of_lit = True
        proteinObj.add_rbbh_lit(litObj)
        self.protein_ids.append(proteinObj.protein_id)
        self.proteinObjs_by_protein_id[proteinObj.protein_id] = proteinObj

    def add_proteinObjs(self, proteinObj_1, proteinObj_2):
        if proteinObj_1.protein_id in self.litObjs_by_lit_id:
            proteinObj_1.part_of_lit = True
            proteinObj_1.effector = True
            proteinObj_1.functions.add(self.litObjs_by_lit_id[proteinObj_1.protein_id].function)
        if proteinObj_2.protein_id in self.litObjs_by_lit_id:
            proteinObj_2.part_of_lit = True
            proteinObj_2.effector = True
            proteinObj_2.functions.add(self.litObjs_by_lit_id[proteinObj_2.protein_id].function)
        if not proteinObj_1.protein_id in self.proteinObjs_by_protein_id:
            self.protein_ids.append(proteinObj_1.protein_id)
            self.proteinObjs_by_protein_id[proteinObj_1.protein_id] = proteinObj_1
        if not proteinObj_2.protein_id in self.proteinObjs_by_protein_id:
            self.protein_ids.append(proteinObj_2.protein_id)
            self.proteinObjs_by_protein_id[proteinObj_2.protein_id] = proteinObj_2
        self.proteinObjs_by_protein_id[proteinObj_1.protein_id].add_rbbh_pcn(proteinObj_2)
        self.proteinObjs_by_protein_id[proteinObj_2.protein_id].add_rbbh_pcn(proteinObj_1)

    def parse_function_f(self):
        count = 0
        for line in read_file(self.function_f):
            col = [x.strip() for x in line.split(",")]  # cleaning leading and trailing whitespaces, and turning it into a list ...
            protein_id = col[0]
            if protein_id.startswith("GPLIN"):
                protein_id = "GPALL.%s" % (protein_id)
            elif protein_id.startswith("GROS"):
                protein_id = "GROST.%s" % (protein_id)
            else:
                pass
            function = col[1]
            species = col[2]
            litObj = LitObj(protein_id, species, function)
            self.add_litObj(litObj)
            count += 1
        print "[+] \t %s sequences parsed ..." % (count)

if __name__ == "__main__":
    __version__ = 0.2

    args = docopt(__doc__)
    MainObj(args)


