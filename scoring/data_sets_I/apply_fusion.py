#
# $Id$
#
# loads ranked lists from different
# models and/or fingerprints and
# apply rank-based fusion
#
# INPUT
# required:
# -i [] : relative input path(s)
# optional:
# -m [] : fusion method (max or ave, default: max)
# -r [] : file containing fingerprints to leave out
# -o [] : relative output path (default: pwd)
# -a : append to the output file (default: overwrite)
# --help : prints usage
#
# OUTPUT: for each target in each data set
#         a file with a list of fusion predictions
#         per fusion prediction: [name, list of 50 scored lists]
#         names: fusion_[fp name] (ml model)
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

import gzip, pickle, math, sys, os, os.path
from collections import defaultdict
from optparse import OptionParser

# import configuration file with global variables
sys.path.insert(0, os.getcwd()+'/../../')
import configuration_file_I as conf

# import functions for scoring step
sys.path.insert(0, os.getcwd()+'/../')
import scoring_functions as scor

# paths
cwd = os.getcwd()
parentpath = cwd+'/../../'
path = cwd+'/'

# prepare command-line option parser
usage = "usage: %prog [options] arg"
parser = OptionParser(usage)
parser.add_option("-i", "--inpath", action="append", dest="inpath", metavar="PATH", help="relative input PATHs")
parser.add_option("-m", "--method", dest="method", help="method for data fusion (max or ave, default: max)")
parser.add_option("-r", "--remove", dest="rm_file", metavar="FILE", help="FILE containing the fingerprints to be left out (default: all fingerprints are read)")
parser.add_option("-o", "--outpath", dest="outpath", metavar="PATH", help="relative output PATH (default: pwd)")
parser.add_option("-a", "--append", dest="do_append", action="store_true", help="append to the output file (default: False)")


######################## MAIN PART ###########################
if __name__=='__main__':

    # read in command line options
    (options, args) = parser.parse_args()
    # required arguments
    if options.inpath:
        inpath = [path+i for i in options.inpath]
        for inp in inpath:
            scor.checkPath(inp, 'input')
    else:
        raise RuntimeError('one or more of the required options was not given!')

    # optional arguments
    method = 'max'
    if options.method:
        if options.method not in ['max', 'ave']:
            raise ValueError('method is unkown. supported methods are: max and ave')
        else:
            method = options.method
    remove_fps = []
    if options.rm_file:
        remove_fps = scor.readFPs(path+options.rm_file)
    outpath = path
    if options.outpath:
        outpath = path+options.outpath
        scor.checkPath(outpath, 'output')
    do_append = False
    if options.do_append: do_append = options.do_append

    # print the fingerprint names
    printfp = True

    # loop over data-set sources
    for dataset in conf.set_data.keys():
        print( dataset)

        # loop over targets
        for target in conf.set_data[dataset]['ids']:
            print( target)

            # load scored lists
            scores = {}
            for inp in inpath: # loop over input paths
                myfile = gzip.open(inp+'/list_'+dataset+'_'+str(target)+'_.pkl.gz', 'rb')
                while 1:
                    try:
                        tmp = pickle.load(myfile)
                    except (EOFError):
                        break
                    else:
                        # check that fp is not in remove_fps list
                        if tmp[0] not in remove_fps:
                            tmp[0] = scor.getName(tmp[0], scores.keys())
                            # input line: [fp_name, list of scored lists]
                            scores[tmp[0]] = tmp[1]
            print( "scored lists read in")
            if len(scores.keys()) < 2:
                print( "number of fingerprints/models < 2, nothing to be done")
                break
            if printfp:
                # determine the name of the fusion
                fpname = 'fusion'
                for k in scores.keys():
                    fpname += '_'+k
                scor.printFPs(scores.keys(), fpname)
                printfp = False

            # FUSION
            # loop over repetitions
            new_scores = []
            for q in range(conf.num_reps):
                # get ranks
                ranks = []
                for k in scores.keys():
                    ranks.append(scor.getRanks(scores[k][q]))

                # do fusion
                new_list = []
                for i in range(len(ranks[0])):
                    current_ranks = [r[i][0] for r in ranks]
                    current_scores = [r[i][1] for r in ranks]
                    if method == 'max': # MAX fusion
                        fused_rank = max(current_ranks)
                        fused_score = max(current_scores)
                    else: # AVE fusion
                        fused_rank = numpy.average(current_ranks)
                        fused_score = numpy.average(current_scores)
                    # store: [fused rank, fused score, internal ID, active/inactive]
                    new_list.append([fused_rank, fused_score, ranks[0][i][-2], ranks[0][i][-1]])
                # sort the new list based on fused ranks
                new_list.sort(reverse=True)
                new_scores.append(new_list)

            # write out the new scores
            if do_append:
                outfile = gzip.open(outpath+'/list_'+dataset+'_'+str(target)+'_'+'.pkl.gz', 'ab+') # binary format
            else:
                outfile = gzip.open(outpath+'/list_'+dataset+'_'+str(target)+'_'+'.pkl.gz', 'wb+') # binary format
            pickle.dump([fpname, new_scores], outfile, 2)
            outfile.close()
        print( "fusion ranking done and ranked list written")
