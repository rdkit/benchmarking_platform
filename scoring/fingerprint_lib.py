#
# $Id$
#
# module to calculate a fingerprint from SMILES

from rdkit import Chem
from rdkit.Chem import MACCSkeys, AllChem
from rdkit.Avalon import pyAvalonTools as fpAvalon
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.ChemicalFeatures import BuildFeatureFactory
from rdkit.Chem.Pharm2DFP import Pharm2DFPCalculator
from rdkit.Chem import rdMolDescriptors

# implemented fingerprints:
# ECFC0 (ecfc0), ECFP0 (ecfp0), MACCS (maccs), 
# atom pairs (ap), atom pairs bit vector (apbv), topological torsions (tt)
# hashed atom pairs (hashap), hashed topological torsions (hashtt) --> with 1024 bits
# ECFP4 (ecfp4), ECFP6 (ecfp6), ECFC4 (ecfc4), ECFC6 (ecfc6) --> with 1024 bits
# FCFP4 (fcfp4), FCFP6 (fcfp6), FCFC4 (fcfc4), FCFC6 (fcfc6) --> with 1024 bits
# Avalon (avalon) --> with 1024 bits
# long Avalon (laval) --> with 16384 bits
# long ECFP4 (lecfp4), long ECFP6 (lecfp6), long FCFP4 (lfcfp4), long FCFP6 (lfcfp6) --> with 16384 bits
# RDKit with path length = 5 (rdk5), with path length = 6 (rdk6), with path length = 7 (rdk7)
# 2D pharmacophore (pharm) ?????????????

nbits = 1024
longbits = 16384

# dictionary
fpdict = {}
fpdict['ecfp0'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 0, nBits=nbits)
fpdict['ecfp2'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 1, nBits=nbits)
fpdict['ecfp4'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=nbits)
fpdict['ecfp6'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 3, nBits=nbits)
fpdict['ecfc0'] = lambda m: AllChem.GetMorganFingerprint(m, 0)
fpdict['ecfc2'] = lambda m: AllChem.GetMorganFingerprint(m, 1)
fpdict['ecfc4'] = lambda m: AllChem.GetMorganFingerprint(m, 2)
fpdict['ecfc6'] = lambda m: AllChem.GetMorganFingerprint(m, 3)
fpdict['fcfp2'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 1, useFeatures=True, nBits=nbits)
fpdict['fcfp4'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, useFeatures=True, nBits=nbits)
fpdict['fcfp6'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 3, useFeatures=True, nBits=nbits)
fpdict['fcfc2'] = lambda m: AllChem.GetMorganFingerprint(m, 1, useFeatures=True)
fpdict['fcfc4'] = lambda m: AllChem.GetMorganFingerprint(m, 2, useFeatures=True)
fpdict['fcfc6'] = lambda m: AllChem.GetMorganFingerprint(m, 3, useFeatures=True)
fpdict['lecfp4'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=longbits)
fpdict['lecfp6'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 3, nBits=longbits)
fpdict['lfcfp4'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, useFeatures=True, nBits=longbits)
fpdict['lfcfp6'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 3, useFeatures=True, nBits=longbits)
fpdict['maccs'] = lambda m: MACCSkeys.GenMACCSKeys(m)
fpdict['ap'] = lambda m: Pairs.GetAtomPairFingerprint(m)
fpdict['tt'] = lambda m: Torsions.GetTopologicalTorsionFingerprintAsIntVect(m)
fpdict['hashap'] = lambda m: rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(m, nBits=nbits)
fpdict['hashtt'] = lambda m: rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(m, nBits=nbits)
fpdict['avalon'] = lambda m: fpAvalon.GetAvalonFP(m, nbits)
fpdict['laval'] = lambda m: fpAvalon.GetAvalonFP(m, longbits)
fpdict['rdk5'] = lambda m: Chem.RDKFingerprint(m, maxPath=5, fpSize=nbits, nBitsPerHash=2)
fpdict['rdk6'] = lambda m: Chem.RDKFingerprint(m, maxPath=6, fpSize=nbits, nBitsPerHash=2)
fpdict['rdk7'] = lambda m: Chem.RDKFingerprint(m, maxPath=7, fpSize=nbits, nBitsPerHash=2)


def CalculateFP(fp_name, smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError('SMILES cannot be converted to a RDKit molecules:', smiles)

    return fpdict[fp_name](m)
