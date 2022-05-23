"""
    generate a list of conformer structures for a given input structure
"""

from  rdkit.Chem import AllChem as Chem
from  rdkit.Chem import rdDistGeom

def  getListConformers( molecule, 
        numConformers, 
        maxAttempts,
        pruneRMSThresh
        ):
    moleculeH  =  Chem.AddHs( molecule )
    randomSeed  =  -1
    clearConfs  =  True
    useRandomCoords  =  False
    boxSizeMult  =  2.0
    randNegEig  =  True
    numZeroFail  =  1
    listIDs  =  rdDistGeom.EmbedMultipleConfs( moleculeH, numConformers, maxAttempts, randomSeed, clearConfs, useRandomCoords, boxSizeMult, randNegEig, numZeroFail, pruneRMSThresh )
    for cID in list( listIDs ):
        Chem.MMFFOptimizeMolecule( moleculeH, confId = cID )
    return  [ list( listIDs ), moleculeH ]
