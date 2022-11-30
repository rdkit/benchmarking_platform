#
# $Id$
#
# file containing the functions for the validation step
#
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

import os
from collections import defaultdict
from rdkit.ML.Scoring import Scoring
from sklearn.metrics import roc_auc_score, average_precision_score


def checkPaths(filepaths):
    '''Checks if the given paths exist'''
    for f in filepaths:
        if not os.path.exists(f):
            raise IOError('path does not exist:', f)


def _readMethods(line):
    '''Helper function for readMethods()'''
    if line:  # if params are provided
        params = []
        for i in line:
            params.append(float(i))
    else:
        raise ValueError("Method requires parameters.")
    return params


# dictionary for readMethods()
read_dict = {}
read_dict['AUC'] = lambda l: EvalMethod(l[0])
read_dict['EF'] = lambda l: EFMethod(l[0], _readMethods(l[1:]), 100)
read_dict['BEDROC'] = lambda l: BEDROCMethod(l[0], _readMethods(l[1:]), 1)
read_dict['RIE'] = lambda l: RIEMethod(l[0], _readMethods(l[1:]), 1)
read_dict['AUROC'] = lambda l: AUROCEvalMethod(l[0])
read_dict['AUPRC'] = lambda l: AUPRCEvalMethod(l[0])


def readMethods(filepath):
    '''Reads the methods names and parameters from a file'''
    try:
        myfile = open(filepath, 'r')
    except:
        raise IOError('file does not exist:', filepath)
    else:
        method_dict = {}
        for line in myfile:
            if line[0] != "#":  # ignore comments
                line = line.rstrip().split()
                method_dict[line[0]] = read_dict[line[0]](line)
        return method_dict


def readFPs(filepath):
    '''Reads a list of fingerprints from a file'''
    try:
        myfile = open(filepath, 'r')
    except:
        raise IOError('file does not exist:', filepath)
    else:
        fps = []
        for line in myfile:
            if line[0] != "#":  # ignore comments
                line = line.rstrip().split()
                fps.append(line[0])
        return fps


def printInputParam(method_dict, inpath):
    '''Prints the input parameters'''
    print("-------------------------------")
    print("PARAMETERS USED")
    print("Validation methods: ")
    for m in method_dict.keys():
        if isinstance(method_dict[m], ParamEvalMethod):
            print(m, "- parameters:", method_dict[m].params)
        else:
            print(m)
    print("")
    print("Input paths:")
    for inp in inpath:
        print(inp)
    print("-------------------------------")


def printFPs(fps):
    '''Prints a list of fingerprints'''
    print("-------------------------------")
    print("FINGERPRINTS CONSIDERED")
    for fp in fps:
        print("   ", fp)
    print("")
    print("-------------------------------")


def getName(fp, fp_names):
    '''Determines the new name of a fingerprint in case
    multiple fingerprints with the same name'''
    # check if fp already exists. if yes, add a number
    if fp in fp_names:
        suffix = 2
        tmp_name = fp + '_' + str(suffix)
        while tmp_name in fp_names:
            suffix += 1
            tmp_name = fp + '_' + str(suffix)
        return tmp_name
    else:
        return fp


# class for handling of evaluation methods
class EvalMethod:
    def __init__(self, name):
        self.method_name = name
        self.names = name

    def addNames(self, results):
        results[self.method_name] = defaultdict(list)

    def calculate(self, score, index):
        return Scoring.CalcAUC(score, index)

    def runMethod(self, results, scores, query, index):
        tmp_list = []
        for k in scores.keys():  # fingerprints
            tmp = self.calculate(scores[k][query], index)
            tmp_list.append([tmp, k])
        # sort list according to the descending score
        tmp_list.sort(reverse=True)
        # store [score, rank]
        for i, l in enumerate(tmp_list):
            # l[1] = fp, l[0] = score, i+1 = rank
            results[self.method_name][l[1]].append([l[0], i + 1])


class AUROCEvalMethod(EvalMethod):
    def __init__(self, name):
        self.method_name = name
        self.names = name

    def calculate(self, score, index):
        scores = [x[0] for x in score]
        acts = [x[index] for x in score]
        return roc_auc_score(acts, scores)


class AUPRCEvalMethod(EvalMethod):
    def __init__(self, name):
        self.method_name = name
        self.names = name

    def calculate(self, score, index):
        scores = [x[0] for x in score]
        acts = [x[index] for x in score]
        return average_precision_score(acts, scores)


class ParamEvalMethod(EvalMethod):
    def __init__(self, name, params, factor):
        EvalMethod.__init__(self, name)
        self.params = params
        self.names = []
        for p in self.params:
            self.names.append(name + str(int(factor * p)))

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
        for i, n in enumerate(self.names):
            # sort list according to the descending score
            tmp_list[i].sort(reverse=True)
            # store [score, rank]
            for j, l in enumerate(tmp_list[i]):
                # l[1] = fp, l[0] = score, j+1 = rank
                results[n][l[1]].append([l[0], j + 1])


class EFMethod(ParamEvalMethod):
    def calculate(self, score, index):
        return Scoring.CalcEnrichment(score, index, self.params)


class BEDROCMethod(ParamEvalMethod):
    def calculate(self, score, index):
        tmp = []
        for p in self.params:
            tmp.append(Scoring.CalcBEDROC(score, index, p))
        return tmp


class RIEMethod(ParamEvalMethod):
    def calculate(self, score, index):
        tmp = []
        for p in self.params:
            tmp.append(Scoring.CalcRIE(score, index, p))
        return tmp
