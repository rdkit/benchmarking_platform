Benchmarking Platform
=====================
presented in S. Riniker, G. Landrum, J. Cheminf., 5, 26 (2013).

GENERAL USAGE NOTES
-------------------
The virtual-screening process implemented by the benchmarking
platform is divided into three steps:
1) Scoring
2) Validation
3) Analysis.
The three steps are run separately and read in the output of the
previous step. In the scoring step, the data from the directories
compounds and query_lists is read in.

The directory compounds contains lists of compounds for 88 targets
from three public data sources: MUV, DUD and ChEMBL. The compound
lists contain the external ID, the internal ID and the SMILES of
each compound.

The directory query_lists contains training lists for each target
with the indices of randomly selected active and inactive molecules.
Training lists with 5, 10 or 20 active molecules are available.
The number of training decoys is 20 % of the decoys in all cases.

The scripts are written in Python and use the open-source
cheminformatics library RDKit (www.rdkit.org).

Running a script with the option [--help] gives a description of the 
required and optional input parameters of the script.
