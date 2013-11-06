#
# $Id:
#
# file containing functions for the scoring step with machine-learning methods
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
from rdkit.ML.Data import DataUtils
import numpy
from multiprocessing import Pool
from sklearn.ensemble import RandomForestClassifier, forest
from sklearn.tree import tree
from sklearn.naive_bayes import BernoulliNB

### FOR SKLEARN VERSION 0.13 ###

# HELPER FUNCTIONS FOR RANDOM FOREST
def _balanced_parallel_build_trees(n_trees, forest, X, y, sample_weight, sample_mask, X_argsorted, seed, verbose):
    """Private function used to build a batch of trees within a job"""
    from sklearn.utils import check_random_state
    from sklearn.utils.fixes import bincount
    import random
    MAX_INT = numpy.iinfo(numpy.int32).max
    random_state = check_random_state(seed)

    trees = []
    for i in xrange(n_trees):
        if verbose > 1:
            print("building tree %d of %d" % (i+1, n_trees))
        seed = random_state.randint(MAX_INT)

        tree = forest._make_estimator(append = False)
        tree.set_params(compute_importances=forest.compute_importances)
        tree.set_params(random_state = check_random_state(seed))

        if forest.bootstrap:
            n_samples = X.shape[0]
            if sample_weight is None:
                curr_sample_weight = numpy.ones((n_samples,), dtype=numpy.float64)
            else:
                curr_sample_weight = sample_weight.copy()

            ty = list(enumerate(y))
            indices = DataUtils.FilterData(ty, val=1, frac=0.5, col=1, indicesToUse=0, indicesOnly=1)[0]
            indices2 = random_state.randint(0, len(indices), len(indices))
            indices = [indices[j] for j in indices2]
            sample_counts = bincount(indices, minlength=n_samples)

            curr_sample_weight *= sample_counts
            curr_sample_mask = sample_mask.copy()
            curr_sample_mask[sample_counts==0] = False

            tree.fit(X, y, sample_weight=curr_sample_weight, sample_mask=curr_sample_mask, X_argsorted=X_argsorted, check_input=False)
            tree.indices = curr_sample_mask
        else:
            tree.fit(X, y, sample_weight=sample_weight, sample_mask=sample_mask, X_argsorted=X_argsorted, check_input=False)
        trees.append(tree)
    return trees

def getNumpy(inlist):
    outlist = []
    for i in inlist:
        arr = numpy.zeros((3,), tree.DTYPE)
        DataStructs.ConvertToNumpyArray(i[1], arr)
        outlist.append(arr)
    return outlist

def readMLFile(ml_dict, read_dict, filepath):
    '''Reads file with the parameters of the machine-learning method
    and stores it in a dictionary'''
    try:
        myfile = open(filepath, 'r')
    except:
        raise IOError('file does not exist:', filepath)
    else:
        for line in myfile:
            l = line.rstrip().split()
        if len(l) != 2:
            raise ValueError('Wrong number of arguments in ML file:', line)
        if l[0] in read_dict:
            ml_dict[l[0]] = read_dict[l[0]](l[1])
        else:
            raise KeyError('Wrong parameter in naive Bayes file:', line)
    return ml_dict

