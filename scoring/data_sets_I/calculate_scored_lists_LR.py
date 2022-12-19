#
# calculates fingerprints and scores lists
# based on the predicted probability
#
# INPUT
# required:
# -n [] : number of query mols
# -f [] : fingerprint to build the logistic regression with
# optional:
# -o [] : relative output path (default: pwd)
# -a : append to the output file (default: overwrite)
# -s [] : similarity metric (default: Dice,
#         other options: Tanimoto, Cosine, Russel, Kulczynski,
#         McConnaughey, Manhattan, RogotGoldberg)
# -r [] : file containing the logistic regression info
#          default parameters: penalty='l2', dual=0 (false), C=1.0,
#          fit_intercept=1 (true), intercept_scaling=1.0,
#          class_weight=None, tol=0.0001
# --help : prints usage
#
# OUTPUT: for each target in each data set
#         a file with a list (1 element) of LR prediction
#         per LR prediction: [name, list of 50 scored lists]
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
import pickle, gzip, sys, os, os.path, numpy
from collections import defaultdict
from optparse import OptionParser
from sklearn.linear_model import LogisticRegression

# import configuration file with global variables
sys.path.insert(0, os.getcwd()+'/../../')
import configuration_file_I as conf

# import functions for scoring step
sys.path.insert(0, os.getcwd()+'/../')
import scoring_functions as scor

# import ML functions
import ml_functions_13 as ml_func

# paths
cwd = os.getcwd()
parentpath = cwd+'/../../'
inpath_cmp = parentpath+'compounds/'
inpath_list = parentpath+'query_lists/data_sets_I/'
path = cwd+'/'

# flag to read in ChEMBL decoys only once
firstchembl = True

# dictionary for readMLFile()
read_dict = {}
read_dict['penalty'] = lambda x: x
read_dict['dual'] = lambda x: bool(x)
read_dict['C'] = lambda x: float(x)
read_dict['fit_intercept'] = lambda x: bool(x)
read_dict['intercept_scaling'] = lambda x: float(x)
read_dict['class_weight'] = lambda x: x
read_dict['tol'] = lambda x: float(x)

# prepare command-line option parser
usage = "usage: %prog [options] arg"
parser = OptionParser(usage)
parser.add_option("-n", "--num", dest="num", type="int", metavar="INT", help="number of query mols")
parser.add_option("-f", "--fingerprint", dest="fp", help="fingerprint to train logistic regression with")
parser.add_option("-o", "--outpath", dest="outpath", metavar="PATH", help="relative output PATH (default: pwd)")
parser.add_option("-s", "--similarity", dest="simil", type="string", metavar="NAME", help="NAME of similarity metric to use (default: Dice, other options are: Tanimoto, Cosine, Russel, Kulczynski, McConnaughey, Manhattan, RogotGoldberg")
parser.add_option("-m", "--ml", dest="ml", metavar="FILE", help="file containing the logistic regression info (default parameters: penalty=l2, dual=0 (false), C=1.0, fit_intercept=1 (true), intercept_scaling=1.0, class_weight=None, tol=0.0001)")
parser.add_option("-a", "--append", dest="do_append", action="store_true", help="append to the output file (default: False)")

############# MAIN PART ########################
if __name__=='__main__':

    # read in command line options
    (options, args) = parser.parse_args()
    # required arguments
    if options.num and options.fp:
        num_query_mols = options.num
        fp_build = options.fp
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
    if outpath_set: scor.checkPath(outpath, 'output')
    scor.checkSimil(simil_metric)
    scor.checkQueryMols(num_query_mols, conf.list_num_query_mols)

    # default machine-learning method variables
    ml_dict = dict(penalty='l2', dual=False, C=1.0, fit_intercept=True, intercept_scaling=1.0, class_weight=None, tol=0.0001)
    if options.ml:
        ml_dict = ml_func.readMLFile(ml_dict, read_dict, path+options.ml)

    # initialize machine-learning method
    ml = LogisticRegression(penalty=ml_dict['penalty'], dual=ml_dict['dual'], C=ml_dict['C'], fit_intercept=ml_dict['fit_intercept'], intercept_scaling=ml_dict['intercept_scaling'], class_weight=ml_dict['class_weight'], tol=ml_dict['tol'])

    # loop over data-set sources
    for dataset in conf.set_data.keys():
        print( dataset)
        # loop over targets
        for target in conf.set_data[dataset]['ids']:
            print( target)

            # read in actives and calculate fps
            actives = []
            for line in gzip.open(inpath_cmp+dataset+'/cmp_list_'+dataset+'_'+str(target)+'_actives.dat.gz', 'r'):
                line=line.decode('UTF-8')
                if line[0] != '#':
                    # structure of line: [external ID, internal ID, SMILES]]
                    line = line.rstrip().split()
                    fp_dict = scor.getFP(fp_build, line[2])
                    # store: [internal ID, dict with fps]
                    actives.append([line[1], fp_dict])
            num_actives = len(actives)
            num_test_actives = num_actives - num_query_mols
            # convert fps to numpy arrays
            np_fps_act = ml_func.getNumpy(actives)

            # read in decoys and calculate fps
            if dataset == 'ChEMBL':
                if firstchembl:
                    decoys = []
                    for line in gzip.open(inpath_cmp+dataset+'/cmp_list_'+dataset+'_zinc_decoys.dat.gz', 'r'):
                        line=line.decode('UTF-8')
                        if line[0] != '#':
                            # structure of line: [external ID, internal ID, SMILES]]
                            line = line.rstrip().split()
                            fp_dict = scor.getFP(fp_build, line[2])
                            # store: [internal ID, dict with fps]
                            decoys.append([line[1], fp_dict])
                    # convert fps to numpy arrays
                    np_fps_dcy = ml_func.getNumpy(decoys)
                    firstchembl = False
            else:
                decoys = []
                for line in gzip.open(inpath_cmp+dataset+'/cmp_list_'+dataset+'_'+str(target)+'_decoys.dat.gz', 'r'):
                    line=line.decode('UTF-8')
                    if line[0] != '#':
                        # structure of line: [external ID, internal ID, SMILES]]
                        line = line.rstrip().split()
                        fp_dict = scor.getFP(fp_build, line[2])
                        # store: [internal ID, dict with fps]
                        decoys.append([line[1], fp_dict])
                # convert fps to numpy arrays
                np_fps_dcy = ml_func.getNumpy(decoys)
            num_decoys = len(decoys)
            print( "molecules read in and fingerprints calculated")

            # open training lists
            training_input = open(inpath_list+dataset+'/training_'+dataset+'_'+str(target)+'_'+str(num_query_mols)+'.pkl', 'rb')
            # to store the scored lists
            scores = defaultdict(list)

            # loop over repetitions
            for q in range(conf.num_reps):
                print( q)
                training_list = pickle.load(training_input)
                test_list = [i for i in range(num_actives) if i not in training_list[:num_query_mols]]
                test_list += [i for i in range(num_decoys) if i not in training_list[num_query_mols:]]

                # list with active/inactive info
                ys_fit = [1]*num_query_mols + [0]*(len(training_list)-num_query_mols)
                # training fps
                train_fps = [np_fps_act[i] for i in training_list[:num_query_mols]]
                train_fps += [np_fps_dcy[i] for i in training_list[num_query_mols:]]
                # fit logistic regression
                ml.fit(train_fps, ys_fit)

                # test fps and molecule info
                test_fps = [np_fps_act[i] for i in test_list[:num_test_actives]]
                test_fps += [np_fps_dcy[i] for i in test_list[num_test_actives:]]
                test_mols = [[actives[i][0], 1] for i in test_list[:num_test_actives]]
                test_mols += [[decoys[i][0], 0] for i in test_list[num_test_actives:]]

                # rank based on probability
                single_score = ml.predict_proba(test_fps)
                # store: [probability, internal ID, active/inactive]
                single_score = [[s[1], m[0], m[1]] for s,m in zip(single_score, test_mols)]
                single_score.sort(reverse=True)
                scores['lr_'+fp_build].append(single_score)

            # write scores to file
            if do_append:
                outfile = gzip.open(outpath+'/list_'+dataset+'_'+str(target)+'.pkl.gz', 'ab+') # binary format
            else:
                outfile = gzip.open(outpath+'/list_'+dataset+'_'+str(target)+'.pkl.gz', 'wb+') # binary format
            for fp in ['lr_'+fp_build]:
                pickle.dump([fp, scores[fp]], outfile, 2)
            outfile.close()
            print( "scoring done and scored lists written")
