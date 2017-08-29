#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: generate_popgenome_calls.py           -l LIST [-h|--help]

    Options:
        -h --help                       show this
        -l, --list FILE                 List file
"""

from __future__ import division
import sys
from docopt import docopt
import os

def read_file(infile):
    if not infile or not os.path.exists(infile):
        sys.exit("[ERROR] - File '%s' does not exist." % (infile))
    print "[+] Parsing %s ..." % (infile)
    with open(infile) as fh:
        for line in fh:
            line = line.replace(r'\r', '\n')
            if not line.startswith("#"):
                yield line.rstrip("\n")

class MainObj():
    def __init__(self, args):
        self.list = args['--list']
        self.length_by_contigs = {}
        self.parse_list()
        self.generate_calls()

    def generate_calls(self):
        for contig, length in self.length_by_contigs.items():
            self.write(contig, length)

    def parse_list(self):
        for line in read_file(self.list):
            col = line.split()
            contig = col[0]
            length = int(col[1])
            self.length_by_contigs[contig] = length

    def write(self, contig, length):
        outfile = "%s.MK.fisher.tsv" % contig
        script_file = "%s.MK.fisher.R" % contig
        vcf = "%s.vcf.gz" % contig
        gff3 = "%s.gff3" % contig
        fasta = "%s.fas" % contig
        output = []
        output.append('require(PopGenome)')
        output.append('snp <- readVCF("popgenome-vcf/%s", numcols=10000, tid="%s", from=1, to=%s, gffpath="popgenome-gff/%s")' % (vcf, contig, length, gff3))
        output.append('snp <- set.synnonsyn(snp, ref.chr="popgenome-fasta/%s")' % fasta)
        output.append('snp <- set.populations(snp, list(c("ERR123957.P4A", "ERR123955.Pa1", "ERR123954.Newton", "Gp4-8.WGS.Nextera_XT", "Gp24.WGS.Nextera_XT", "Gp4-8.WGA.Nextera", "Gp24.WGA.Nextera", "Gp4-8.WGA.Nextera_XT", "Gp12-17.WGA.Nextera", "Gp12-17.WGS.Nextera_XT", "ERR114517.Lindley", "Gp24.WGA.Nextera_XT", "Gp12-17.WGA.Nextera_XT", "ERR123953.Luffness", "ERR123952.Bedale"), c("ERR123956.P5A")), diploid=T)')
        output.append('snp <- set.outgroup(snp, c("ERR123956.P5A"), diploid=TRUE)')
        output.append('coding <- splitting.data(snp, subsites="coding", whole.data=T)')
        output.append('coding <- MKT(coding, do.fisher.test=T)')
        output.append('from.pos <- sapply(coding@region.names,function(x){return(as.numeric(strsplit(x," ")[[1]][1]))})')
        output.append('to.pos <- sapply(coding@region.names,function(x){return(as.numeric(strsplit(x," ")[[1]][3]))})')
        output.append('DATA <- cbind(from.pos, to.pos, coding@MKT)')
        output.append('DATA <- as.data.frame(DATA)')
        output.append('DATA$chrom <- "%s"' % contig)
        output.append('write.table(DATA, file="%s", quote=FALSE, sep="\\t", col.names = NA)' % outfile)
        with open(script_file, 'w') as fh:
            fh.write("\n".join(output) + "\n")

if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    print "[+] Start ..."
    main = MainObj(args)
