from collections import namedtuple
from collections import defaultdict
import numpy as np

ExpressionDatum = namedtuple('ExpressionDatum', 'name length expression')


class ExpresionDataset(object):

  def __init__(self):
    self.exps_ = []


  def addElement(self, datum):
    self.exps_.append(datum)

  def aggregateWith(self, edict):
    aexp = ExpresionDataset()
    adict = defaultdict(list)

    haveNone = self.exps_[0].length is None
    for d in self.exps_:
        k = edict[d.name]
        adict[k].append(d)

    if haveNone:
        aggregatedExps = [ExpressionDatum(k, None, sum([v.expression for v in e])) for k, e in adict.iteritems()]
    else:
        aggregatedExps = [ExpressionDatum(k, sum([v.length for v in e]), sum([v.expression for v in e])) for k, e in adict.iteritems()]

    aexp.exps_ = aggregatedExps
    return aexp

  def normalize(self, method="id"):
    def rpkm(e, tot):
        mult = np.power(10, 9)
        if e.expression < 1.0:
            return 0.0
        return mult * (e.expression / (tot*e.length)) if e.length > 0.0 else 0.0

    if method == "id":
        pass
    elif method == "log":
        slen = len(self.exps_)
        self.exps_ = filter(lambda x: x.expression > 0.0, self.exps_)
        print("Removed {0} of {1} original expression values for log transform".format(slen-len(self.exps_), slen))
        self.exps_ = [ExpressionDatum(e.name, e.length, np.log2(e.expression)) for e in self.exps_]
        print("Minimum expression value is {0}".format(min([e.expression for e in self.exps_])))
    elif method == "rpkm":
        tot = np.sum(d.expression for d in self.exps_)
        self.exps_ = [ExpressionDatum(e.name, e.length, rpkm(e, tot)) for e in self.exps_]
    elif method == "lrpkm":
        self.normalize("rpkm")
        self.normalize("log")

  def zipWithMatching(self, other):
    self.exps_ = sorted(self.exps_, key=lambda x: x.name)
    other.exps_ = sorted(other.exps_, key=lambda x: x.name)
    sind, oind = 0, 0
    slen, olen = len(self.exps_), len(other.exps_)
    matches = []

    while sind < slen and oind < olen:
        sdat = self.exps_[sind]
        odat = other.exps_[oind]
        if sdat.name == odat.name:
            if np.isfinite(sdat.expression) and np.isfinite(odat.expression):
                matches.append((sdat, odat))
            oind += 1
            sind += 1
        elif sdat.name < odat.name:
            sind += 1
        else:
            oind += 1

    return matches

def accession(tname):
    return tname.split('.')[0]

def version(tname):
    return int(tname.split('.')[-1])

def parseExpressionFile(fname, epsilon=0.0):

    dset = ExpresionDataset()

    with open(fname, 'rb') as ifile:
        for l in ifile:
            toks = l.rstrip().split()
            ntoks = len(toks)
            hasLength = ntoks == 3

            dname = accession(toks[0])
            dlength = float(toks[1]) if hasLength else None
            dexp = float(toks[-1])

            dset.addElement(ExpressionDatum(dname, dlength, dexp + epsilon))

    return dset
