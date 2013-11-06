#
# $Id$
#
# file containing functions for analysis step
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

import os

##################### HELPER FUNCTIONS #########################

def checkPath(filepath):
    '''Checks if a path exists'''
    if not os.path.exists(filepath):
        raise IOError('path does not exist:', filepath)

def writeHeader(outfile, mk):
    '''Writes a header to a file'''
    outfile.write("# ")
    outfile.writelines("%s std " % k for k in mk)
    outfile.write("\r\n") # for Windows

def readFile(myfile):
    '''Reads a given file and stores the data in a dictionary'''
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

def readMethods(filepath):
    '''Reads method names from a file'''
    try:
        myfile = open(filepath, 'r')
    except:
        raise IOError('file does not exist:', filepath)
    else:
        methods = []
        for line in myfile:
            if line[0] != "#": # ignore comments
                line = line.rstrip().split()
                methods.append(line[0])
        return methods
