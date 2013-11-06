#
# $Id$
#
# file containing the functions for the scoring step
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

import os, sys, operator
from rdkit import DataStructs

# import the fingerprint library
import fingerprint_lib

def checkFPFile(filepath):
    '''Checks if file containing fingerprint names exists
    and reads the fingerprints'''
    try:
        myfile = open(filepath, 'r')
    except:
        raise IOError('file does not exist:', filepath)
    else:
        fp_names = []
        for line in myfile:
            line = line.rstrip().split()
            fp_names.append(line[0])
        return fp_names

def checkPath(filepath, name):
    '''Checks if a filepath exists'''
    if not os.path.exists(filepath):
        raise IOError(name, 'path does not exist:', filepath)

def checkSimil(simil):
    '''Checks if the chosen similarity metric is supported'''
    simil_list = ['Dice', 'Tanimoto', 'Cosine', 'Russel', 'Kulczynski', 'McConnaughey', 'Manhattan', 'RogotGoldberg']
    if simil not in simil_list:
        raise ValueError('provided similarity metric not supported:', simil)

def checkQueryMols(num, list_num_query_mols):
    '''Checks if the chosen number of query molecules is supported'''
    if num not in list_num_query_mols:
        raise ValueError('provided number of query molecules not supported:', num)

def getFPDict(fp_names, smiles):
    '''Gets the fingerprints from the fingerprint library
    and stores them in a dictioanry'''
    fp_dict = {}
    for fp in fp_names:
        fp_dict[fp] = fingerprint_lib.CalculateFP(fp, smiles)
    return fp_dict

def getFP(fp_name, smiles):
    '''Gets fingerprint from fingerprint library'''
    return fingerprint_lib.CalculateFP(fp_name, smiles)

# dictionary for similarity measures
simil_dict = {}
simil_dict['Dice'] = lambda x,y: sorted(DataStructs.BulkDiceSimilarity(x,y), reverse=True)
simil_dict['Tanimoto'] = lambda x,y: sorted(DataStructs.BulkTanimotoSimilarity(x, y), reverse=True)
simil_dict['Cosine'] = lambda x,y: sorted(DataStructs.BulkCosineSimilarity(x,y), reverse=True)
simil_dict['Russel'] = lambda x,y: sorted(DataStructs.BulkRusselSimilarity(x,y), reverse=True)
simil_dict['Kulczynski'] = lambda x,y: sorted(DataStructs.BulkKulczynskiSimilarity(x,y), reverse=True)
simil_dict['McConnaughey'] = lambda x,y: sorted(DataStructs.BulkMcConnaugheySimilarity(x,y), reverse=True)
simil_dict['Manhattan'] = lambda x,y: sorted(DataStructs.BulkAllBitSimilarity(x,y), reverse=True)
simil_dict['RogotGoldberg'] = lambda x,y: sorted(DataStructs.BulkRogotGoldbergSimilarity(x,y), reverse=True)

def getBulkSimilarity(fp, fp_list, simil):
    '''Calculate the bulk similarity for a given list of fingerprints'''
    return simil_dict[simil](fp,fp_list)

# helper functions for the fusion
def printFPs(fps, fpname):
    '''Prints a list of fingerprints'''
    print "-------------------------------"
    print "FUSION DONE FOR:"
    for fp in fps:
        print fp,
    print ""
    print "Name of fusion:", fpname
    print "-------------------------------"

def getName(fp, fp_names):
    '''Determines the new name of a fingerprint in case
    multiple fingerprints with the same name'''
    # check if fp already exists. if yes, add a number
    if fp in fp_names:
        suffix = 2
        tmp_name = fp+'_'+str(suffix)
        while tmp_name in fp_names:
            suffix += 1
            tmp_name = fp+'_'+str(suffix)
        return tmp_name
    else:
        return fp

def readFPs(filepath):
    '''Reads a list of fingerprints from a file'''
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

def getRanks(probas):
    '''Add the ranks for a ranked list'''
    num_mol = len(probas)
    ranks = [[num_mol-i] + j for i,j in enumerate(probas)]
    # sort based on internal ID
    ranks.sort(key=operator.itemgetter(-2))
    return ranks
