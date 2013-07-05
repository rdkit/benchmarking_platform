#
# calculates fingerprints and scores lists
# based on similarity
#
# INPUT
# required:
# -n [] : number of query mols
# -f [] : file containing fingerprint names
# optional:
# -o [] : relative output path (default: pwd)
# -a : append to the output file (default: overwrite)
# -s [] : similarity metric (default: Dice, 
#         other options: Tanimoto, Cosine, Russel, Kulczynski, 
#         McConnaughey, Manhattan, RogotGoldberg)
# --help : prints usage
#
# OUTPUT: for each target in each data set
#         a file with a list of fingerprints
#         per fingerprint: [name, list of 50 scored lists]
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

from rdkit import Chem, DataStructs
import cPickle, gzip, sys, os, os.path
from collections import defaultdict
from optparse import OptionParser 

# import fingerprint library
import fingerprint_lib

# import configuration file with global variables
sys.path.insert(0, os.getcwd()+'/../')
import configuration_file as conf

# paths
cwd = os.getcwd()
parentpath = cwd+'/../'
inpath1 = parentpath+'compounds/'
inpath2 = parentpath+'query_lists/'
path = cwd+'/'

# flag to read in ChEMBL decoys only once
firstchembl = True

############### HELPER FUNCTIONS #################
def checkFPFile(filepath):
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

def checkOutputPath(filepath):
    if not os.path.exists(filepath):
        raise IOError('output path does not exist:', filepath)

def checkSimil(simil):
    simil_list = ['Dice', 'Tanimoto', 'Cosine', 'Russel', 'Kulczynski', 'McConnaughey', 'Manhattan', 'RogotGoldberg']
    if simil not in simil_list:
        raise ValueError('provided similarity metric not supported:', simil)

def checkQueryMols(num):
    if num not in conf.list_num_query_mols:
        raise ValueError('provided number of query molecules not supported:', num)

def getFPDict(fp_names, smiles):
    fp_dict = {}
    for fp in fp_names:
        fp_dict[fp] = fingerprint_lib.CalculateFP(fp, smiles)
    return fp_dict

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
    return simil_dict[simil](fp,fp_list)

# prepare command-line option parser
usage = "usage: %prog [options] arg"
parser = OptionParser(usage)
parser.add_option("-n", "--num", dest="num", type="int", metavar="INT", help="number of query mols")
parser.add_option("-f", "--fingerprints", dest="fp_file", metavar="FILE", help="FILE containing fingerprint names")
parser.add_option("-o", "--outpath", dest="outpath", metavar="PATH", help="relative output PATH (default: pwd)")
parser.add_option("-s", "--similarity", dest="simil", type="string", metavar="NAME", help="NAME of similarity metric to use (default: Dice, other options are: Tanimoto, Cosine, Russel, Kulczynski, McConnaughey, Manhattan, RogotGoldberg")
parser.add_option("-a", "--append", dest="do_append", action="store_true", help="append to the output file (default: False)")

############## END OF HELPER FUNCTIONS #################

############# MAIN PART ########################
if __name__=='__main__':

    # read in command line options
    (options, args) = parser.parse_args()
    # required arguments
    if options.num and options.fp_file: 
        num_query_mols = options.num
        fp_file = path+options.fp_file
    else:
        raise RuntimeError('one or more of the required options was not given!')

    # optional arguments
    do_append = False
    if options.do_append: do_append = options.do_append
    simil_metric = 'Dice'
    if options.simil: simil_metric = options.simil
    outpath = path
    outpath_set = False
    if options.outpath:
        outpath_set = True
        outpath = path+options.outpath

    # check for sensible input
    fp_names = checkFPFile(fp_file)
    if not fp_names: raise ValueError('No fingerprints given in', fp_file)
    if outpath_set: checkOutputPath(outpath)
    checkSimil(simil_metric)
    checkQueryMols(num_query_mols)

    # loop over data-set sources
    for dataset in conf.set_data.keys():
        print dataset
        # loop over targets
        for target in conf.set_data[dataset]['ids']:
            print target

            # read in actives and calculate fps
            actives = []
            for line in gzip.open(inpath1+dataset+'/cmp_list_'+dataset+'_'+str(target)+'_actives.dat.gz', 'r'):
                if line[0] != '#': 
                    # structure of line: [external ID, internal ID, SMILES]]
                    line = line.rstrip().split()
                    fp_dict = getFPDict(fp_names, line[2])
                    # store: [internal ID, dict with fps]
                    actives.append([line[1], fp_dict])
            num_actives = len(actives)
            num_test_actives = num_actives - num_query_mols

            # read in decoys and calculate fps
            if dataset == 'ChEMBL':
                if firstchembl:
                    decoys = []
                    for line in gzip.open(inpath1+dataset+'/cmp_list_'+dataset+'_zinc_decoys.dat.gz', 'r'):
                        if line[0] != '#': 
                            # structure of line: [external ID, internal ID, SMILES]]
                            line = line.rstrip().split()
                            fp_dict = getFPDict(fp_names, line[2])
                            # store: [internal ID, dict with fps]
                            decoys.append([line[1], fp_dict])
                    firstchembl = False
            else:
                decoys = []
                for line in gzip.open(inpath1+dataset+'/cmp_list_'+dataset+'_'+str(target)+'_decoys.dat.gz', 'r'):
                    if line[0] != '#': 
                        # structure of line: [external ID, internal ID, SMILES]]
                        line = line.rstrip().split()
                        fp_dict = getFPDict(fp_names, line[2])
                        # store: [internal ID, dict with fps]
                        decoys.append([line[1], fp_dict])
            num_decoys = len(decoys)
            print "molecules read in and fingerprints calculated"

            # open training lists
            training_input = open(inpath2+dataset+'/training_'+dataset+'_'+str(target)+'_'+str(num_query_mols)+'.pkl', 'r')
            # to store the scored lists
            scores = defaultdict(list)

            # loop over repetitions
            for q in range(conf.num_reps):
                training_list = cPickle.load(training_input)
                test_list = [i for i in range(num_actives) if i not in training_list[:num_query_mols]]
                test_list += [i for i in range(num_decoys) if i not in training_list[num_query_mols:]]
                # loop over fps
                single_score = defaultdict(list)
                for fp in fp_names:
                    query_fps = [actives[i][1][fp] for i in training_list[:num_query_mols]]
                    # test_list: first actives then decoys
                    test_fps = [[actives[i][0], actives[i][1][fp], 1] for i in test_list[:num_test_actives]]
                    test_fps += [[decoys[i][0], decoys[i][1][fp], 0] for i in test_list[num_test_actives:]]
                    for tmp_mol in test_fps:
                        tmp_score = getBulkSimilarity(tmp_mol[1], query_fps, simil_metric)
                        # use max fusion
                        # store : [similarity, internal ID, active/inactive]
                        single_score[fp].append([tmp_score[0], tmp_mol[0], tmp_mol[2]]) 
                    # rank list according to similarity
                    scores[fp].append(sorted(single_score[fp], reverse=True))

            # write scores to file
            if do_append:
                outfile = gzip.open(outpath+'/list_'+dataset+'_'+str(target)+'_'+'.pkl.gz', 'ab+') # binary format
            else:
                outfile = gzip.open(outpath+'/list_'+dataset+'_'+str(target)+'_'+'.pkl.gz', 'wb+') # binary format
            for fp in fp_names:
                cPickle.dump([fp, scores[fp]], outfile, 2)
            outfile.close()
            print "scoring done and scored lists written"
