#
# $Id$
#
# configuration file for benchmarking platform
#
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

# hardcoded global variables
list_num_query_mols = [5, 10, 20]
num_reps = 50 # number of repetitions
percent_dcy = 0.2 # percentage of decoys used for training
p_value = 0.05 # confidence level for statistical analysis

# collection of data sets
muv_ids = [466, 548, 600, 644, 652, 689, 692, 712, 713, 733, 737, 810, 832, 846, 852, 858, 859]
dud_ids = ['ace', 'ache', 'ar', 'cdk2', 'cox2', 'dhfr', 'egfr', 'er_agonist', 'fgfr1', 'fxa', 'gpb', 'gr', 'hivrt', 'inha', 'na', 'p38', 'parp', 'pdgfrb', 'sahh', 'src', 'vegfr2']
chembl_ids = [11359, 28, 11536, 8, 10434, 12670, 20014, 234, 12261, 12209, 25, 36, 43, 219, 130, 105, 11336, 20174, 126, 11225, 12252, 11682, 134, 116, 11265, 10475, 12679, 10579, 11575, 18061, 237, 276, 11534, 10198, 10498, 12911, 12968, 100579, 100126, 10378, 10417, 10752, 10773, 11631, 10927, 11085, 11442, 11279, 11488, 12840]

set_data = {}
set_data['MUV'] = dict(fullname='MUV', ids=muv_ids, prefix='aid', suffix='_actives.sdf', dcy_prefix='aid', dcy_suffix='_decoys.sdf', propName='PUBCHEM_COMPOUND_CID', dcy_propName='PUBCHEM_COMPOUND_CID')
set_data['DUD'] = dict(fullname='DUD', ids=dud_ids, prefix='', suffix='_clustered_3D_MM.sdf', dcy_prefix='DUD_', dcy_suffix='_decoys_ID_pass_MWPass_I_MM.sdf', propName='id', dcy_propName='Mol_Title')
set_data['ChEMBL'] = dict(fullname='ChEMBL-sereina/diverse_100', ids=chembl_ids, prefix='Target_no_', suffix='.sdf', dcy_name='decoys_10000_zinc.sdf', propName='Name', dcy_propName='_Name')
