#
# loads scored lists with actives and
# decoys, and calculates different
# validation methods
#
# INPUT
# required:
# -m [] : file containing the methods 
#         implemented methods are: AUC, BEDROC ([alpha] optional),
#         RIE ([alpha] optional), EF ([percentage] optional)
# optional:
# -i [] : absolute input path (default: pwd/../scoring)
# -o [] : absolute output path (default: pwd)
# -r [] : file containing fingerprints to leave out
# --help : prints usage
#
# OUTPUT: for each target in each dataset
#         a dictionary with methods and for each method
#         a dictionary with fingerprints and for each fp
#         a list with the [score, rank]
#
#  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met: 
#
#     * Redistributions of source code must retain the above copyright 
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following 
#       disclaimer in the documentation and/or other materials provided 
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
#       nor the names of its contributors may be used to endorse or promote 
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

import gzip, cPickle, math, sys, os, os.path
from collections import defaultdict
from optparse import OptionParser
from rdkit.ML.Scoring import Scoring

# import configuration file with global variables
sys.path.insert(0, os.getcwd()+'/../')
import configuration_file as conf

# paths
cwd = os.getcwd()
parentpath = cwd+'/../'
path = cwd+'/'

##################### HELPER FUNCTIONS #########################

def checkPaths(filepaths):
    for f in filepaths:
        if not os.path.exists(f):
            raise IOError('path does not exist:', f)

# helper method for readMethods()
def _readMethods(line):
    if line: # if params are provided
        params = []
        for i in line: params.append(float(i))
    else:
        raise ValueError("Method requires parameters.")
    return params

# dictionary for readMethods()
read_dict = {}
read_dict['AUC'] = lambda l: EvalMethod(l[0])
read_dict['EF'] = lambda l: EFMethod(l[0], _readMethods(l[1:]), 100)
read_dict['BEDROC'] = lambda l: BEDROCMethod(l[0], _readMethods(l[1:]), 1)
read_dict['RIE'] = lambda l: RIEMethod(l[0], _readMethods(l[1:]), 1)

def readMethods(filepath):
    try:
        myfile = open(filepath, 'r')
    except:
        raise IOError('file does not exist:', filepath)
    else:
        method_dict = {}
        for line in myfile:
            if line[0] != "#": # ignore comments
                line = line.rstrip().split()
                method_dict[line[0]] = read_dict[line[0]](line)
        return method_dict

def readFPs(filepath):
    try:
        myfile = open(filepath, 'r')
    except:
        raise IOError('file does not exist:', filepath)
    else:
        fps = []
        for line in myfile:
            if line[0] != "#": # ignore comments
                line = line.rstrip().split()
                fps.append(line[0])
        return fps

def printInputParam(method_dict, inpath):
    print "-------------------------------"
    print "PARAMETERS USED"
    print "Validation methods: "
    for m in method_dict.keys():
        if isinstance(method_dict[m], ParamEvalMethod):
            print m, "- parameters:", method_dict[m].params
        else:
            print m
    print ""
    print "Input paths:"
    for inp in inpath:
        print inp
    print "-------------------------------"

def printFPs(fps):
    print "-------------------------------"
    print "FINGERPRINTS CONSIDERED"
    for fp in fps:
        print fp,
    print ""
    print "-------------------------------"

# class for handling of evaluation methods
class EvalMethod:
    def __init__(self, name):
        self.method_name = name
        self.names = name
    def addNames(self, results):
        results[self.method_name] = defaultdict(list)
    def calculate(self, score, index):
        return Scoring.CalcAUC(score,index)
    def runMethod(self, results, scores, query, index):
        tmp_list = []
        for k in scores.keys(): # fingerprints
            tmp = self.calculate(scores[k][query], index)
            tmp_list.append([tmp, k])
        # sort list according to the descending score
        tmp_list.sort(reverse=True)
        # store [score, rank]
        for i,l in enumerate(tmp_list):
            # l[1] = fp, l[0] = score, i+1 = rank
            results[self.method_name][l[1]].append([l[0], i+1])

class ParamEvalMethod(EvalMethod):
    def __init__(self, name, params, factor):
        EvalMethod.__init__(self, name)
        self.params = params
        self.names = []
        for p in self.params:
            self.names.append(name + str(int(factor*p)))
    def addNames(self, results):
        for n in self.names: 
            results[n] = defaultdict(list)
    def runMethod(self, results, scores, query, index):
        tmp_list = [[] for i in range(len(self.names))]
        # loop over fingerprints
        for k in scores.keys(): 
            tmp = self.calculate(scores[k][query], index)
            # loop over parameters
            for i in range(len(self.names)):
                tmp_list[i].append([tmp[i], k])
        # loop over parameters
        for i,n in enumerate(self.names):
            # sort list according to the descending score
            tmp_list[i].sort(reverse=True)
            # store [score, rank]
            for j,l in enumerate(tmp_list[i]):
                # l[1] = fp, l[0] = score, j+1 = rank
                results[n][l[1]].append([l[0], j+1])

class EFMethod(ParamEvalMethod):
    def calculate(self, score, index):
        return Scoring.CalcEnrichment(score,index,self.params)

class BEDROCMethod(ParamEvalMethod):
    def calculate(self, score, index):
        tmp = []
        for p in self.params:
            tmp.append(Scoring.CalcBEDROC(score,index,p))
        return tmp

class RIEMethod(ParamEvalMethod):
    def calculate(self, score, index):
        tmp = []
        for p in self.params:
            tmp.append(Scoring.CalcRIE(score,index,p))
        return tmp

# prepare command-line option parser
usage = "usage: %prog [options] arg"
parser = OptionParser(usage)
parser.add_option("-m", "--methods", dest="filename", metavar="FILE", help="FILE containing the methods")
parser.add_option("-i", "--inpath", action="append", dest="inpath", metavar="PATH", help="relative input PATH (default: pwd/../scoring)")
parser.add_option("-o", "--outpath", dest="outpath", metavar="PATH", help="relative output PATH (default: pwd)")
parser.add_option("-r", "--remove", dest="rm_file", metavar="FILE", help="FILE containing the fingerprints to be left out (default: all fingerprints are read)")

############## END OF HELPER FUNCTIONS #########################

######################## MAIN PART ###########################
if __name__=='__main__':

    # read in command line options
    (options, args) = parser.parse_args()
    # required arguments
    if options.filename:
        methods_file = path+options.filename
    else:
        raise RuntimeError('one or more of the required options was not given!')
    # read method file
    method_dict = readMethods(methods_file)
    if not method_dict: raise ValueError('No methods given in', methods_file)

    # optional arguments
    inpath = parentpath+'scoring/'
    if options.inpath: 
        inpath = [path+i for i in options.inpath]
        checkPaths(inpath)
    outpath = path
    if options.outpath: 
        outpath = path+options.outpath
        checkPaths([outpath])
    remove_fps = []
    if options.rm_file:
        remove_fps = readFPs(path+options.rm_file)

    # print input parameters
    printInputParam(method_dict, inpath)

    # print the fingerprint names
    printfp = True

    # loop over data-set sources
    for dataset in conf.set_data.keys():
        print dataset
        # output directory
        outdir = outpath+'/'+dataset
        if not os.path.exists(outdir): os.makedirs(outdir)

        # loop over targets
        for target in conf.set_data[dataset]['ids']:
            print target

            # load scored lists
            scores = {}
            for inp in inpath: # loop over input paths
                myfile = gzip.open(inp+'/list_'+dataset+'_'+str(target)+'_.pkl.gz', 'r')
                while 1:
                    try:
                        tmp = cPickle.load(myfile)
                    except (EOFError):
                        break
                    else:
                        # check that fp is not in remove_fps list
                        if tmp[0] not in remove_fps:
                            # check if fp already exists. if yes, add a number
                            if tmp[0] in scores:
                                suffix = 2
                                tmp_name = tmp[0]+'_'+str(suffix)
                                while tmp_name in scores:
                                    suffix += 1
                                    tmp_name = tmp[0]+'_'+str(suffix)
                                tmp[0] = tmp_name
                            # input line: [fp_name, list of scored lists]
                            scores[tmp[0]] = tmp[1]
            print "scored lists read in"
            if printfp:
                printFPs(scores.keys())
                printfp = False

            # prepare to store results
            results = {}
            for m in method_dict.keys():
                method_dict[m].addNames(results)

            # loop of repetitions
            for q in range(conf.num_reps):
                # loop over evaluation methods
                for m in method_dict.keys(): 
                    method_dict[m].runMethod(results, scores, q, -1)

            print "validation methods calculated"

            # write results
            outf = gzip.open(outdir+'/validation_'+str(target)+'.pkl.gz', 'wb+')
            cPickle.dump(results, outf, 2)
            outf.close()

            print "results written out"
