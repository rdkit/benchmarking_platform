#
# $Id$
#
# loads average performance per target
# and summarizes them in a single file
# for each fingerprint
#
# INPUT
# optional:
# -i [] : relative input path (default: pwd)
# --help : prints usage
#
# OUTPUT: for each fingerprint
#         a summary for each method and target
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

import gzip, cPickle, sys, os, os.path
from collections import defaultdict
from optparse import OptionParser

# import configuration file with global variables
sys.path.insert(0, os.getcwd()+'/../')
import configuration_file as conf

# paths
cwd = os.getcwd()
path = cwd+'/'

##################### HELPER FUNCTIONS #########################

def checkPath(filepath):
    if not os.path.exists(filepath):
        raise IOError('path does not exist:', filepath)

def readFile(myfile):
    input_data = {}
    fpkeys = []
    for line in myfile:
        if line.startswith("#"): # read fp names
            line = line.rstrip().split()
            for i in range(1,len(line),2): # the first element is #
                fpkeys.append(line[i]) 
        else:
            line = line.rstrip().split()
            num = len(line)
            if num > 0: # not an empty line
                input_data[line[0]] = [float(line[i]) for i in range(1,num)]
    return input_data, fpkeys

def writeHeader(outfile, mk):
    outfile.write("# ")
    outfile.writelines("%s std " % k for k in mk)
    outfile.write("\r\n") # for Windows

# prepare command-line option parser
usage = "usage: %prog [options] arg"
parser = OptionParser(usage)
parser.add_option("-i", "--inpath", dest="inpath", metavar="PATH", help="relative input PATH (default: pwd)")

############## END OF HELPER FUNCTIONS #########################

######################## MAIN PART ###########################
if __name__=='__main__':

    # read in command line options
    (options, args) = parser.parse_args()

    # optional arguments
    inpath = path
    if options.inpath:
        inpath = path+options.inpath
        checkPath(inpath)
    outpath = inpath

    summary = {}

    # loop over dataset sources
    for dataset in conf.set_data.keys():
        print dataset
        # input path
        inpath = outpath+'/'+dataset

        # loop over targets
        for target in conf.set_data[dataset]['ids']:
            print target

            # load results
            results, fpkeys = readFile(open(inpath+'/target_'+str(target)+'.txt', 'r'))
            methodkeys = results.keys()

            # if summary is not yet set: prepare it
            if len(summary) == 0:
                for m in methodkeys:
                    summary[m] = defaultdict(list)

            # fill results into summary dictionary
            for m in methodkeys:
                for i,k in enumerate(fpkeys):
                    j = i*2
                    summary[m][k].append([target, results[m][j], results[m][j+1]])

    # write out
    outdir = outpath+'/fp_summary/'
    if not os.path.exists(outdir): os.makedirs(outdir)
    # loop over fingerprints
    for k in fpkeys:
        outfile = open(outdir+'summary_'+k+'.txt', 'w')
        writeHeader(outfile, methodkeys)
        # loop over targets
        for t in range(len(summary[methodkeys[0]][k])):
            curr_id = summary[methodkeys[0]][k][t][0]
            outfile.write("%s " % curr_id)
            # loop over methods
            for m in methodkeys:
                element = summary[m][k][t]
                outfile.write("%.3f %.3f " % (element[1], element[2]))
            outfile.write("\r\n") # for Windows
        outfile.close()
