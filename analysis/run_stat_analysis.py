# 
# $Id$
#
# executes an R script for each method;
# reads the generated csv files with 
# fp pairs, mean rank, and 100 resampled p-values
#
# INPUT
# optional:
# -i [] : relative input path (default: pwd/../validation/[data set name])
# --help : prints usage
#
# OUTPUT: correlation tables for each method
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

import csv, os, os.path, sys
import operator, subprocess
from optparse import OptionParser

# import configuration file with global variables
sys.path.insert(0, os.getcwd()+'/../')
import configuration_file as conf

# paths
cwd = os.getcwd()
path = cwd + '/'
Rpath = '/usr/prog/R/2.15.2_gcc/bin/Rscript'

# list of evaluation methods
methods = ['AUC', 'EF1', 'EF5', 'BEDROC20', 'BEDROC100', 'RIE20', 'RIE100']

##################### HELPER FUNCTIONS #########################

def checkPath(filepath):
    if not os.path.exists(filepath):
        raise IOError('path does not exist:', filepath)

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
    inpath = path+'stat_analysis/'
    if options.inpath:
        inpath = path+options.inpath+'/stat_analysis/'
        checkPath(inpath)
    outpath = inpath

    # output file(s)
    outfile = open(outpath+'correlation_tables.dat', 'w')

    # loop over methods
    for m in methods:
        print m

        # run R script
        R_input = inpath+'fp_ranking_'+m+'.csv'
        R_output1 = outpath+'friedman_fp_ranking_'+m+'.csv'
        R_output2 = outpath+'friedman_resampled_p_values_'+m+'.csv'
        proc = subprocess.call([Rpath, path+'Friedman_Test.Target_x_FP_Data.R', R_input, R_output1, R_output2])

        if not proc: # 0 if successfull

            outfile.write("%s\n" % m)

            # read in csv file
            friedman = []
            with open(R_output2, 'r') as csvfile:
                myreader = csv.reader(csvfile, delimiter=',')
                for row in myreader: 
                    if row[0] != "": # header leave out
                        # row has the format:
                        # [index, fp1, fp2, mean rank 1, mean rank 2, p-value1, ..., p-value100]

                        # check if mean rank 1 < mean rank 2, else exchange
                        if float(row[3]) > float(row[4]):
                            # exchange fingerprint names
                            tmp = row[1]
                            row[1] = row[2]
                            row[2] = tmp
                            # exchange mean ranks
                            tmp = row[3]
                            row[3] = row[4]
                            row[4] = tmp

                        # change format of the row to:
                        # [mean rank 1, mean rank 2, fp1, fp2, [p-values 1..100]]
                        tmp = [float(row[3]), float(row[4]), row[1], row[2]]
                        tmp.append([float(x) for x in row[5:]])
                        friedman.append(tmp)

            # sort the table (first column, than second column)
            friedman = sorted(friedman)
            # order list of fps with respect to mean rank
            fps = []
            for f in friedman:
                if f[2] not in fps: fps.append(f[2])
            fps.append(friedman[-1][3])

            # column description
            outfile.write("\t")
            outfile.writelines("%s " % x for x in fps[1:])
            outfile.write("\n")

            # loop over table
            nextfp = 1
            outfile.write("%s\t" % fps[0])
            for f in friedman:
                if f[2] == fps[nextfp]: # new fingerprint
                    outfile.write("\n") # first end previous line
                    outfile.write("%s\t" % f[2])
                    for i in range(nextfp):
                        outfile.write("  ")
                    nextfp += 1

                # check if all p-values are below confidence level
                all = True
                none = True
                for i in f[4]:
                    if none and i < conf.p_value: 
                        none = False
                    if all and i > conf.p_value:
                        all = False
                if all and not none: # all p-values smaller than 0.05
                    outfile.write("- ") # significant difference
                elif not all and not none:
                    outfile.write("o ") # p-values below and above 0.05
                elif none and not all:
                    outfile.write("X ") # no significant difference
            outfile.write("\n")
            outfile.write("%s\n\n" % fps[-1])

    outfile.close()
