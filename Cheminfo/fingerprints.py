#!/usr/bin/env python

from rdkit.Chem.Fingerprints import FingerprintMols #Linear (Daylight-like) fingerprint
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs #For atom-pairs fingerprint
from rdkit.Chem.AtomPairs import Torsions #For topological torsion fingerprint
from rdkit.Chem import AllChem #For Mrogan (also known as circular) fingerprint

def Fptlinear(rdkmol, minPath=1, maxPath=7, fpSize=2048, bitsPerHash=2, useHs=0, tgtDensity=0.3, minSize=64):
    '''When changing any one argument from default, we need to pass all 7 arguments to RDKit FingerprintMol'''
    if fpSize == 2048:
        return FingerprintMols.FingerprintMol(rdkmol)
    else:
        return FingerprintMols.FingerprintMol(rdkmol,minPath=minPath, maxPath=maxPath, fpSize=fpSize,
                                              bitsPerHash=bitsPerHash, useHs=useHs, tgtDensity=tgtDensity, minSize=minSize)

def FptMACCS(rdkmol):
    return MACCSkeys.GenMACCSKeys(rdkmol)
    

def FptAtomPairs(rdkmol, fptype='bit'):
    if fptype.lower() == 'bit':
        return Pairs.GetAtomPairFingerprintAsBitVect(rdkmol)
    else:
        return Pairs.GetAtomPairFingerprint(rdkmol)

def FptTorsion(rdkmol):
    return Torsions.GetTopologicalTorsionFingerprintAsIntVect(rdkmol)

def FptCircular(rdkmol, radius=2 , fptype='bit', feature=False, nBits=2048):
    if fptype.lower() == 'bit':
        return AllChem.GetMorganFingerprintAsBitVect(rdkmol,int(radius),nBits=nBits, useFeatures=feature)
    else:
        return AllChem.GetMorganFingerprint(rdkmol,int(radius), useFeatures=feature)
    

