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
from rdkit.Chem import rdMolDescriptors


def FoldedRDKFingerprintCountBased(mol, fpSize=1024, **kwargs):
    bitInfo = {}
    unfolded = Chem.UnfoldedRDKFingerprintCountBased(mol,
                                                     branchedPaths=False,
                                                     minPath=3,
                                                     maxPath=3,
                                                     bitInfo=bitInfo,
                                                     **kwargs)
    res = {}
    newBitInfo = defaultdict(list)
    for k, b in unfolded.GetNonzeroElements().items():
        res[k % fpSize] = b
        newBitInfo[k % fpSize].extend(bitInfo[k])
    return res, newBitInfo


# implemented fingerprints:
# mfc0 (mfc0), mfp0 (mfp0), MACCS (maccs),
# atom pairs (ap), atom pairs bit vector (apbv), topological torsions (tt)
# hashed atom pairs (hashap), hashed topological torsions (hashtt) --> with 1024 bits
# mfp2 (mfp2), mfp3 (mfp3), mfc2 (mfc2), mfc3 (mfc3) --> with 1024 bits
# fmfp2 (fmfp2), fmfp3 (fmfp3), fmfc2 (fmfc2), fmfc3 (fmfc3) --> with 1024 bits
# Avalon (avalon) --> with 1024 bits
# long Avalon (laval) --> with 16384 bits
# long mfp2 (lmfp2), long mfp3 (lmfp3), long fmfp2 (lfmfp2), long fmfp3 (lfmfp3) --> with 16384 bits
# RDKit with path length = 5 (rdk5), with path length = 6 (rdk6), with path length = 7 (rdk7)
# 2D pharmacophore (pharm) ?????????????

nbits = 1024
longbits = 16384

# dictionary
fpdict = {}
fpdict['mfp0'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 0, nBits=nbits)
fpdict['mfp1'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 1, nBits=nbits)
fpdict['mfp2'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 2, nBits=nbits)
fpdict['mfp3'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 3, nBits=nbits)
fpdict['mfc0'] = lambda m: AllChem.GetMorganFingerprint(m, 0)
fpdict['mfc1'] = lambda m: AllChem.GetMorganFingerprint(m, 1)
fpdict['mfc2'] = lambda m: AllChem.GetMorganFingerprint(m, 2)
fpdict['mfc3'] = lambda m: AllChem.GetMorganFingerprint(m, 3)
fpdict['fmfp1'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 1, useFeatures=True, nBits=nbits)
fpdict['fmfp2'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 2, useFeatures=True, nBits=nbits)
fpdict['fmfp3'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 3, useFeatures=True, nBits=nbits)
fpdict['fmfc1'] = lambda m: AllChem.GetMorganFingerprint(
    m, 1, useFeatures=True)
fpdict['fmfc2'] = lambda m: AllChem.GetMorganFingerprint(
    m, 2, useFeatures=True)
fpdict['fmfc3'] = lambda m: AllChem.GetMorganFingerprint(
    m, 3, useFeatures=True)
fpdict['lmfp2'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 2, nBits=longbits)
fpdict['lmfp3'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 3, nBits=longbits)
fpdict['lfmfp2'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 2, useFeatures=True, nBits=longbits)
fpdict['lfmfp3'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(
    m, 3, useFeatures=True, nBits=longbits)
fpdict['maccs'] = lambda m: MACCSkeys.GenMACCSKeys(m)
fpdict['ap'] = lambda m: rdMolDescriptors.GetAtomPairFingerprint(m)
fpdict['tt'] = lambda m: rdMolDescriptors.GetTopologicalTorsionFingerprint(m)
fpdict[
    'hashap'] = lambda m: rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
        m, nBits=nbits)
fpdict[
    'hashtt'] = lambda m: rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(
        m, nBits=nbits)
fpdict[
    'lhashap'] = lambda m: rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
        m, nBits=longbits)
fpdict[
    'lhashtt'] = lambda m: rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(
        m, nBits=longbits)
fpdict['avalon'] = lambda m: fpAvalon.GetAvalonFP(m, nbits)
fpdict['laval'] = lambda m: fpAvalon.GetAvalonFP(m, longbits)
fpdict['rdk5'] = lambda m: Chem.RDKFingerprint(
    m, maxPath=5, fpSize=nbits, nBitsPerHash=2)
fpdict['rdk6'] = lambda m: Chem.RDKFingerprint(
    m, maxPath=6, fpSize=nbits, nBitsPerHash=2)
fpdict['rdk7'] = lambda m: Chem.RDKFingerprint(
    m, maxPath=7, fpSize=nbits, nBitsPerHash=2)
fpdict['lrdk5'] = lambda m: Chem.RDKFingerprint(
    m, maxPath=5, fpSize=longbits, nBitsPerHash=2)
fpdict['lrdk6'] = lambda m: Chem.RDKFingerprint(
    m, maxPath=6, fpSize=longbits, nBitsPerHash=2)
fpdict['lrdk7'] = lambda m: Chem.RDKFingerprint(
    m, maxPath=7, fpSize=longbits, nBitsPerHash=2)


def CalculateFP(fp_name, smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError('SMILES cannot be converted to a RDKit molecules:',
                         smiles)

    return fpdict[fp_name](m)
