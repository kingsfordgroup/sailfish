"""AggregateToGeneLevel

Usage: AggregateToGeneLevel.py --gtf=<gtffile> --texp=<texpfile> --out=<ofile> --method=<method>

Options:
  -h                 Display help
  --gtf=<gtffile>    The gene feature file
  --texp=<texpfile>  The transcript level expression file
  --method=<method>  The method by which multiple transcript expression should be aggreaged
  --out=<out>        File where gene-level expressions should be written

"""
from docopt import docopt

import os
import cPickle
import ExpressionTools
from collections import defaultdict
from HTSeq import GFF_Reader


def main(args):


    texps = ExpressionTools.parseExpressionFile(args['--texp'])
    print(texps)
    texps.normalize(args['--method'])

    ### build the transcript => gene and gene => transcript dictionaries
    gtffile = args['--gtf']
    
    pkfile = gtffile+'.pk'
    if os.path.exists(pkfile):
        with open(pkfile, 'rb') as ifile:
            tran2gene = cPickle.load(ifile)
            gene2tran = cPickle.load(ifile)
    else:
        gfeats = GFF_Reader(args['--gtf'])
        tran2gene = {}
        gene2tran = defaultdict(list)
        for f in gfeats:
            tname = ExpressionTools.accession(f.attr['transcript_id'])
            gname = f.attr['gene_id']
            tran2gene[tname] = gname
            gene2tran[gname].append(tname)

        with open(pkfile, 'wb') as dfile:
            cPickle.dump(tran2gene, dfile)
            cPickle.dump(gene2tran, dfile)

    gexps = texps.aggregateWith(tran2gene)
    #gexps.normalize(args['--method'])

    with open(args['--out'],'wb') as ofile:
        for d in gexps.exps_:
            if d.length is None:
                ofile.write('{0}\t{1}\n'.format(d.name, d.expression))
            else:
                ofile.write('{0}\t{1}\t{2}\n'.format(d.name, d.length, d.expression))

if __name__ == "__main__":
    args = docopt(__doc__, version="AggregateToGeneLevel 1.0")
    main(args)