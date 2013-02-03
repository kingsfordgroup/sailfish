"""PlotCorrelation

Usage: PlotCorrelation.py EXP1 --name=<name> [--norm=<norm>] EXP2 --name=<name> [--norm=<norm>]

Arguments:
  EXP1     The first expression file
  EXP2     The second expression file

Options:
  -h                 Display help
  --name=<name>      The name of the expression estimation method
  --norm=<norm>      The normalization to perform

"""
from docopt import docopt

from matplotlib import pyplot as plt
import numpy as np
import scipy as sp

import ExpressionTools

def main(args):
    print(args)
    exp1 = ExpressionTools.parseExpressionFile(args['EXP1'])
    name1 = args['--name'][0]
    norm1 = args['--norm'][0]
    exp1.normalize(norm1)

    exp2 = ExpressionTools.parseExpressionFile(args['EXP2'])
    name2 = args['--name'][1]
    norm2 = args['--norm'][1]
    exp2.normalize(norm2)

    matches = exp1.zipWithMatching(exp2)
    x,y = zip(*[ (e[0].expression, e[1].expression) for e in matches ])
    print(min(x),max(x))
    print(min(y),max(y))
    
    import scipy.stats
    pr = sp.stats.pearsonr(x,y)[0]
    sr = sp.stats.spearmanr(x,y)[0]
    print("There were {0} datapoints".format(len(matches)))
    print("Pearson r = {0}".format(pr))
    print("Spearman r = {0}".format(sr))

    font = {'family' : 'normal',
    'weight' : 'bold',
    'size'   : 18}

    plt.rc('font', **font)

    nstr = {"id" : "", "log" : "($\\log_2$)", "rpkm" : "(RPKM)", "lrpkm" : "($\\log_2$ RPKM)"}

    plt.title("Corr. Pearson = {0:.3}, Spearman = {1:.3}".format(pr,sr))  
    plt.xlabel("{0}{1}".format(name1, nstr[norm1]))
    plt.ylabel("{0}{1}".format(name2, nstr[norm2]))
    plt.hexbin(x,y)
    plt.show()


if __name__ == "__main__":
    args = docopt(__doc__, version="PlotCorrelation 1.0")
    main(args)