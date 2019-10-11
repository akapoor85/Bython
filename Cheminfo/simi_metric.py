#!/usr/bin/env python

from rdkit import DataStructs
import sys

def GetSimilarity(rdkmol1, rdkmol2, metric='tanimoto'):
    '''
    mol1 and mol2 are RDKit fingerprint objects for molecule
    '''
    valid_metric = ('tanimoto', 'dice', 'cosine', 'sokal', 'russel', 'kulczynski', 'mcconnaughey')
    if metric.lower() == 'tanimoto':
        return DataStructs.TanimotoSimilarity(rdkmol1, rdkmol2)
    elif metric.lower() == 'dice':
        return DataStructs.DiceSimilarity(rdkmol1, rdkmol2)
    elif metric.lower() == 'cosine':
        return DataStructs.CosineSimilarity(rdkmol1, rdkmol2)
    elif metric.lower() == 'sokal':
        return DataStructs.SokalSimilarity(rdkmol1, rdkmol2)
    elif metric.lower() == 'russel':
        return DataStructs.RusselSimilarity(rdkmol1, rdkmol2)
    elif metric.lower() == 'kulczynski':
        return DataStructs.KulczynskiSimilarity(rdkmol1, rdkmol2)
    elif metric.lower() == 'mcconnaughey':
        return DataStructs.McConnaugheySimilarity(rdkmol1, rdkmol2)
    #elif metric.lower() == 'tversky': #Was returning error
    #    return DataStructs.TverskySimilarity(rdkmol1, rdkmol2)
    else:
        sys.exit('***ERROR: Unrecognized similarity metric: %s***.\n Use one of %s' % (metric,', '.join(valid_metric)))
    
    