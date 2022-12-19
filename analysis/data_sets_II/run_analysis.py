#
# $Id$
#
# loads scored results from different
# evaluation methods and averages them
#
# INPUT
# optional:
# -i [] : absolute input path (default: pwd/../validation)
# -o [] : absolute output path (default: pwd)
# --help : prints usage
#
# OUTPUT: for each target in each dataset
#         a summary for each fp and method;
#         for each method a list with ranks
#         for each fingerprint and repetition
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
import numpy as np
from scipy import special, stats
from collections import defaultdict
from optparse import OptionParser

# import configuration file with global variables
sys.path.insert(0, os.getcwd()+'/../../')
import configuration_file_II as conf

# import analysis function
sys.path.insert(0, os.getcwd()+'/../')
import analysis_functions as ana_func

#paths
cwd = os.getcwd()
path = cwd+'/'


# prepare command-line option parser
usage = "usage: %prog [options] arg"
parser = OptionParser(usage)
parser.add_option("-i", "--inpath", dest="inpath", metavar="PATH", help="relative input PATH (default: pwd/../validation)")
parser.add_option("-o", "--outpath", dest="outpath", metavar="PATH", help="relative output PATH (default: pwd)")

######################## MAIN PART ###########################
if __name__=='__main__':

    # read in command line options
    (options, args) = parser.parse_args()

    # optional arguments
    inpath = path+'../validation/'
    if options.inpath:
        inpath = path+options.inpath
        ana_func.checkPath(inpath)
    outpath = path
    if options.outpath:
        outpath = path+options.outpath
        ana_func.checkPath(outpath)
    tmppath = inpath

    ranks = {}
    header = False # header written to csv file

    # output directories and input directory
    outdir = outpath+'/ChEMBL'
    if not os.path.exists(outdir): os.makedirs(outdir)
    outdir_csv = outpath+'/stat_analysis'
    if not os.path.exists(outdir_csv): os.makedirs(outdir_csv)
    inpath = tmppath+'/ChEMBL'

    # loop over targets
    for target in conf.set_data:
        print( target)

        # load results
        validation = pickle.load(gzip.open(inpath+'/validation_'+str(target)+'.pkl.gz', 'rb'))
        methodkeys = list(validation.keys())
        fpkeys = validation[methodkeys[0]].keys()

        # if ranks is not yet set: prepare it
        if len(ranks) == 0:
            for m in methodkeys:
                ranks[m] = defaultdict(list)

        # prepare writing out
        # if header is not yet written to csv file
        if not header:
            for m in methodkeys:
                outfile_csv = open(outdir_csv+'/fp_ranking_'+m+'.csv', 'w')
                outfile_csv.write("\"fp\",\"target\",\"rank\",\"repeat\"\n")
                outfile_csv.close()
            header = True

        # header of average score per target
        outfile = open(outdir+'/target_'+str(target)+'.txt', 'w')
        ana_func.writeHeader(outfile, fpkeys)

        # calculate average score
        # loop over methods
        for m in methodkeys:
            outfile.write("%s\t" % m) # average score file
            outfile_csv = open(outdir_csv+'/fp_ranking_'+m+'.csv', 'a') # csv file: append!
            # loop over fingerprints
            for k in fpkeys:
                # csv rank file
                for i,v in enumerate(validation[m][k]):
                    outfile_csv.write("\"%s\",\"%s\",\"%i\",\"%i\"\n" % (k, str(target), v[1], i+1))
                # average score file
                tmp = []
                for i in validation[m][k]:
                    ranks[m][k].append(i[1]) # store ranks
                    tmp.append(i[0])
                tmp = np.array(tmp)
                ave = np.average(tmp)
                std = np.std(tmp)
                outfile.write("%.3f %.3f " % (ave, std))
            outfile.write("\r\n") # for Windows

            outfile_csv.close()
        outfile.close()

    # calculate average ranks and write out
    outfile = open(outpath+'/average_rank_fps.txt', 'w')
    ana_func.writeHeader(outfile, fpkeys)
    averages = defaultdict(list)
    for m in methodkeys:
        outfile.write("%s " % m)
        for k in fpkeys:
            tmp = np.array(ranks[m][k])
            ave = np.average(tmp)
            std = np.std(tmp)
            outfile.write("%.2f %.2f " % (ave, std))
            averages[m].append([ave, k])
        outfile.write("\r\n") # for Windows
    outfile.close()
