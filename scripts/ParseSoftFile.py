"""ParseSoftFile.py

Usage:
  ParseSoftFile.py (--soft=<sfile>) (--out=<ofile>) (--gfile=<gfile>) [--samps=<samps>]

Options:
  -h --help        Help
  --soft=<sfile>   SOFT format file
  --samps=<samps>  comma separated list of samples to consider
  --gfile=<gfile>  file containing the transcript to gene mapping
  --out=<ofile>    output file
"""
from docopt import docopt
from collections import defaultdict
from collections import namedtuple
import itertools 

import HTSeq

class SoftFile(object):
    
    def __init__(self, consideredSamples, transcriptMap, ifile):
        self._transcriptMap = transcriptMap
        self._geneMap = transcriptMap
        self._transcriptsForGene = defaultdict(set)
        self._platformTable = defaultdict(list)
        self._reversePlatformTable = {}
        self._sampleTables = []
        with open(ifile,'rb') as inputFile:
            self.__readFromFile(consideredSamples, inputFile)

    def __parsePlatformTable(self, ifile):
        l = ifile.readline()
        while True:
            l = ifile.readline()
            if l.startswith('!platform_table_end'): return
            ind, orf, assayID, posDescrMfr, gbList = l.rstrip().split()
            ids = set(gbList.split(','))
            self._transcriptsForGene[ orf ] |= ids
            for i in ids:
                self._platformTable[ orf ].append(ind)
            self._reversePlatformTable[ind] = ids

    def __parseSampleTable(self, ifile):
        SampleEntry = namedtuple('SampleEntry', ['val', 'expressed'], verbose=False)
        sampleTable = {}
        readData = False
        while True:
            l = ifile.readline()
            if l.startswith('!sample_table_begin'):
                l = ifile.readline() # eat one extra line of header info
                l = ifile.readline()
                readData = True

            elif l.startswith('!sample_table_end'): 
                self._sampleTables.append(sampleTable)
                return

            if readData:
                ind, val, flag = l.rstrip().split()
                sampleTable[ind] = SampleEntry( val = float(val), expressed = True if flag == 'P' else False )

    def __readFromFile(self, consideredSamples, ifile):
        while True:
            l = ifile.readline()
            if l == '': return

            if l.startswith('!platform_table_begin'):
                self.__parsePlatformTable(ifile)
                for k,v in self._platformTable.iteritems():
                    if k in self._geneMap:
                        GEOSet = self._transcriptsForGene[k]
                        RefSeqSet = self._geneMap[k]
                        print(GEOSet)
                        print(RefSeqSet)
                    else:
                        print("{0} is not in RefSeq".format(k))


            elif l.startswith('^SAMPLE'):
                toks = l.rstrip().split(' = ')
                if len(consideredSamples) == 0 or (toks[-1] in consideredSamples):
                    self.__parseSampleTable(ifile)

    def averageExpressions(self):
        geneExpressionDict = {}
        transcripts = {}
        for name,inds in self._platformTable.iteritems():
            numTabs = 0
            tot = 0.0
            numExpressed = 0.0
            numOcc = 0.0
            for t in self._sampleTables:
                for v in inds:
                    if v in t:
                        numOcc += 1
                        if t[v].expressed:
                            tot += t[v].val
                            numTabs += 1
                            numExpressed += 1
            precentExpressed = (numExpressed / numOcc) if numOcc > 0.0 else 0.0
            if precentExpressed >= 0.75 and name in self._geneMap: 
                res = set(itertools.chain(*[ self._reversePlatformTable[i] for i in inds ]))
                print(res)
                if res.issuperset(self._geneMap[name]):#len(res) >= len(self._geneMap[name]):
                    geneExpressionDict[name] = tot / float(numTabs)
                    transcripts[name] = set(itertools.chain(*[ self._reversePlatformTable[i] for i in inds ]))
                else:
                    print("annotation set from RefSeq was {0}".format(','.join(self._geneMap[name])))
                    print("annotation set from GEO was {0}\n\n".format(','.join(res)))

        return (geneExpressionDict, transcripts)

###
#  Functions for extracing information from transcript
#  names (Note: I hate all these naming schemes).
##
def getTranscriptAccession(transcriptName):
    return transcriptName.split('.')[0]

def getTranscriptVersion(transcriptName):
    return transcriptNames.split('.')[1]

def main(args):
    sfile = args['--soft']
    ofile = args['--out']
    gfile = args['--gfile']
    sline = args['--samps']

    t2g = {}
    gtfFile = HTSeq.GFF_Reader( gfile )
    for record in gtfFile:
        # Right now we're only interested in the transcript accession,
        # so strip off version information
        taccession = getTranscriptAccession(record.attr['transcript_id'])
        t2g[ taccession ] = record.attr['gene_id']

    g2t = defaultdict(set)
    for transcript, gene in t2g.iteritems():
        g2t[ gene ].add(transcript)

    consideredSamples = set(sline.split(','))
    sdata = SoftFile(consideredSamples, g2t, sfile)

    avg, tr = sdata.averageExpressions()
    with open(ofile,'wb') as output:
        for k,v in avg.iteritems():#sdata.averageExpressions().iteritems():
            output.write('{0}\t{1}\n'.format(k,v))
    
    # Print out the transcripts mapped to each gene
    # with open('test','wb') as output:
    #     for k,v in tr.iteritems():
    #         output.write('{0}\t{1}\n'.format(k, ','.join(v)))

if __name__ == "__main__":
  arguments = docopt(__doc__, version="ParseSoftFile.py v1.0")
  main(arguments)
