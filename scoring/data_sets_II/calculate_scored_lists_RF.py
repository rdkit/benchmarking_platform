#
# calculates fingerprints and scores lists
# based on the predicted probability
#
# INPUT
# required:
# -f [] : fingerprint to build the random forest with
# optional:
# -o [] : relative output path (default: pwd)
# -a : append to the output file (default: overwrite)
# -s [] : similarity metric (default: Dice,
#         other options: Tanimoto, Cosine, Russel, Kulczynski,
#         McConnaughey, Manhattan, RogotGoldberg)
# -r [] : file containing the random forest info
#          default parameters: criterion=gini, max_depth=10,
#          max_features=auto (=sqrt), num_estimators=100,
#          min_samples_split=2, min_samples_leaf=1, n_jobs=1
# --help : prints usage
#
# OUTPUT: for each target in each data set
#         a file with a list (1 element) of RF prediction
#         per RF prediction: [name, list of 50 scored lists]
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
import cPickle, gzip, sys, os, os.path, numpy
from collections import defaultdict
from optparse import OptionParser
from sklearn.ensemble import RandomForestClassifier, forest
from sklearn.tree import tree
from rdkit.ML.Data import DataUtils
from multiprocessing import Pool

# import configuration file with global variables
sys.path.insert(0, os.getcwd()+'/../../')
import configuration_file_II as conf

# import functions for scoring step
sys.path.insert(0, os.getcwd()+'/../')
import scoring_functions as scor

# import ML functions
import ml_functions_13 as ml_func

# paths
cwd = os.getcwd()
parentpath = cwd+'/../../'
inpath_cmp = parentpath+'compounds/'
inpath_list = parentpath+'query_lists/data_sets_II/ChEMBL/'
path = cwd+'/'

# flag to read in ChEMBL decoys only once
firstchembl = True

# dictionary for readMLFile()
read_dict = {}
read_dict['criterion'] = lambda x: x
read_dict['max_depth'] = lambda x: int(x)
read_dict['max_features'] = lambda x: x
read_dict['num_estimators'] = lambda x: int(x)
read_dict['min_samples_split'] = lambda x: int(x)
read_dict['min_samples_leaf'] = lambda x: int(x)
read_dict['n_jobs'] = lambda x: int(x)

# forest._parallel_build_trees = ml_func._balanced_parallel_build_trees

# prepare command-line option parser
usage = "usage: %prog [options] arg"
parser = OptionParser(usage)
parser.add_option("-f", "--fingerprint", dest="fp", help="fingerprint to train random forest with")
parser.add_option("-o", "--outpath", dest="outpath", metavar="PATH", help="relative output PATH (default: pwd)")
parser.add_option("-s", "--similarity", dest="simil", type="string", metavar="NAME", help="NAME of similarity metric to use (default: Dice, other options are: Tanimoto, Cosine, Russel, Kulczynski, McConnaughey, Manhattan, RogotGoldberg")
parser.add_option("-m", "--ml", dest="ml", metavar="FILE", help="file containing the random forest info (default parameters: criterion=gini, max_depth=10, max_features=auto (=sqrt), num_estimators=100, min_samples_split=2, min_samples_leaf=1, n_jobs=1)")
parser.add_option("-a", "--append", dest="do_append", action="store_true", help="append to the output file (default: False)")

############# MAIN PART ########################
if __name__=='__main__':

    # read in command line options
    (options, args) = parser.parse_args()
    # required arguments
    if options.fp:
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

    # default machine-learning method variables
    ml_dict = dict(criterion='gini', max_features='auto', n_jobs=1, max_depth=10, min_samples_split=2, min_samples_leaf=1, num_estimators=100)
    if options.ml:
        ml_dict = ml_func.readMLFile(ml_dict, read_dict, path+options.ml)

    # initialize machine-learning method
    ml = RandomForestClassifier(criterion=ml_dict['criterion'], max_features=ml_dict['max_features'], min_samples_split=ml_dict['min_samples_split'], max_depth=ml_dict['max_depth'], min_samples_leaf=ml_dict['min_samples_leaf'], n_estimators=ml_dict['num_estimators'], n_jobs=ml_dict['n_jobs'])

    # loop over targets
    for target in conf.set_data:
        print target

        # read in training actives and calculate fps
        actives = cPickle.load(open(inpath_cmp+'ChEMBL_II/Target_no_'+str(target)+'.pkl', 'r'))
        for k in actives.keys():
            for i,m in enumerate(actives[k]):
                fp_dict = scor.getFP(fp_build, m[1])
                actives[k][i] = [str(target)+'_'+str(k)+'_A_'+str(i+1), fp_dict]

        # read in test actives and calculate fps
        div_actives = []
        for line in gzip.open(inpath_cmp+'ChEMBL/cmp_list_ChEMBL_'+str(target)+'_actives.dat.gz', 'r'):
            if line[0] != '#':
                # structure of line: [external ID, internal ID, SMILES]]
                line = line.rstrip().split()
                fp_dict = scor.getFP(fp_build, line[2])
                # store: [internal ID, dict with fps]
                div_actives.append([line[1], fp_dict])
        num_test_actives = conf.num_div_act - 1
        # convert fps to numpy arrays
        np_fps_div_act = ml_func.getNumpy(div_actives)

        # read in decoys and calculate fps
        if firstchembl:
            decoys = []
            for line in gzip.open(inpath_cmp+'ChEMBL/cmp_list_ChEMBL_zinc_decoys.dat.gz', 'r'):
                if line[0] != '#':
                    # structure of line: [external ID, internal ID, SMILES]]
                    line = line.rstrip().split()
                    fp_dict = scor.getFP(fp_build, line[2])
                    # store: [internal ID, dict with fps]
                    decoys.append([line[1], fp_dict])
            # convert fps to numpy arrays
            np_fps_dcy = ml_func.getNumpy(decoys)
            firstchembl = False
            num_decoys = len(decoys)
        print "molecules read in and fingerprints calculated"

        # open training and test lists
        training_input = open(inpath_list+'/training_'+str(target)+'.pkl', 'r')
        test_input = open(inpath_list+'/test_'+str(target)+'.pkl', 'r')
        # to store the scored lists
        scores = defaultdict(list)
        # loop over repetitions
        for q in actives.keys():
            print q
            num_actives = len(actives[q])
            np_fps_act = ml_func.getNumpy(actives[q])
            training_list = cPickle.load(training_input)
            test_list = cPickle.load(test_input)
            test_list += [i for i in range(num_decoys) if i not in training_list[num_actives:]]

            # list with active/inactive info
            ys_fit = [1]*num_actives + [0]*(len(training_list)-num_actives)
            # training fps
            train_fps = [actives[q][i][1] for i in range(num_actives)]
            np_train_fps = np_fps_act + [np_fps_dcy[i] for i in training_list[num_actives:]]
            # fit random forest
            ml.fit(np_train_fps, ys_fit)

            # test fps and molecule info
            test_fps = [div_actives[i][1] for i in test_list[:num_test_actives]]
            test_fps += [decoys[i][1] for i in test_list[num_test_actives:]]
            np_test_fps = [np_fps_div_act[i] for i in test_list[:num_test_actives]]
            np_test_fps += [np_fps_dcy[i] for i in test_list[num_test_actives:]]
            test_mols = [[div_actives[i][0], 1] for i in test_list[:num_test_actives]]
            test_mols += [[decoys[i][0], 0] for i in test_list[num_test_actives:]]

            # calculate similarity with standard fp
            std_simil = []
            for fp in test_fps:
                tmp_simil = scor.getBulkSimilarity(fp, train_fps, simil_metric)
                tmp_simil.sort(reverse=True)
                std_simil.append(tmp_simil[0])

            # rank based on probability (and second based on similarity)
            single_score = ml.predict_proba(np_test_fps)
            # store: [probability, similarity, internal ID, active/inactive]
            single_score = [[m[1], s, t[0], t[1]] for m,s,t in zip(single_score,std_simil,test_mols)]
            single_score.sort(reverse=True)
            scores['rf_'+fp_build].append(single_score)

        # write scores to file
        if do_append:
            outfile = gzip.open(outpath+'/list_'+str(target)+'.pkl.gz', 'ab+') # binary format
        else:
            outfile = gzip.open(outpath+'/list_'+str(target)+'.pkl.gz', 'wb+') # binary format
        for fp in ['rf_'+fp_build]:
            cPickle.dump([fp, scores[fp]], outfile, 2)
        outfile.close()
        print "scoring done and scored lists written"
