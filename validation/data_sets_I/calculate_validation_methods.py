#
# $Id$
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
# -i [] : relative input path (default: pwd/../scoring)
# -o [] : relative output path (default: pwd)
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
sys.path.insert(0, os.getcwd()+'/../../')
import configuration_file_I as conf

# import validation functions
sys.path.insert(0, os.getcwd()+'/../')
import validation_functions as vfunc

# paths
cwd = os.getcwd()
parentpath = cwd+'/../'
path = cwd+'/'

# prepare command-line option parser
usage = "usage: %prog [options] arg"
parser = OptionParser(usage)
parser.add_option("-m", "--methods", dest="filename", metavar="FILE", help="FILE containing the methods")
parser.add_option("-i", "--inpath", action="append", dest="inpath", metavar="PATH", help="relative input PATH (default: pwd/../scoring)")
parser.add_option("-o", "--outpath", dest="outpath", metavar="PATH", help="relative output PATH (default: pwd)")
parser.add_option("-r", "--remove", dest="rm_file", metavar="FILE", help="FILE containing the fingerprints to be left out (default: all fingerprints are read)")


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
    method_dict = vfunc.readMethods(methods_file)
    if not method_dict: raise ValueError('No methods given in', methods_file)

    # optional arguments
    inpath = parentpath+'scoring/'
    if options.inpath:
        inpath = [path+i for i in options.inpath]
        vfunc.checkPaths(inpath)
    outpath = path
    if options.outpath:
        outpath = path+options.outpath
        vfunc.checkPaths([outpath])
    remove_fps = []
    if options.rm_file:
        remove_fps = vfunc.readFPs(path+options.rm_file)

    # print input parameters
    vfunc.printInputParam(method_dict, inpath)

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
                myfile = gzip.open(inp+'/list_'+dataset+'_'+str(target)+'.pkl.gz', 'r')
                while 1:
                    try:
                        tmp = cPickle.load(myfile)
                    except (EOFError):
                        break
                    else:
                        # check that fp is not in remove_fps list
                        if tmp[0] not in remove_fps:
                            tmp[0] = vfunc.getName(tmp[0], scores.keys())
                            # input line: [fp_name, list of scored lists]
                            scores[tmp[0]] = tmp[1]
            print "scored lists read in"
            if printfp:
                vfunc.printFPs(scores.keys())
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
