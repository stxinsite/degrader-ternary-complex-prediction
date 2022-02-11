"""
    script to move a specified structure 
    (either POI with attached ligand or a ligase with attached ligand)
    in a way that both ligands are within close proximity to each other.
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms
#from rdkit.Chem import Draw
#from rdkit import rdBase

from rdkit.Chem import rdchem
from rdkit.Chem import rdmolops



def rot_ar_x( radi, transX ):
    return  np.array([[1, 0, 0, transX],
                      [0, np.cos(radi), -np.sin(radi), 0],
                      [0, np.sin(radi), np.cos(radi), 0],
                     [0, 0, 0, 1]], dtype=np.double)
 
def rot_ar_y( radi, transY ):
    return  np.array([[np.cos(radi), 0, np.sin(radi), 0],
                      [0, 1, 0, transY],
                      [-np.sin(radi), 0, np.cos(radi), 0],
                     [0, 0, 0, 1]], dtype=np.double)
 
def rot_ar_z( radi, transZ ):
    return  np.array([[np.cos(radi), -np.sin(radi), 0, 0],
                      [np.sin(radi), np.cos(radi), 0, 0],
                      [0, 0, 1, transZ],
                     [0, 0, 0, 1]], dtype=np.double)


def  getCenter( molecule ):
    center  =  [ 0, 0, 0 ]
    count  =  0
    for atoms in molecule.GetConformers():
        for i in range( 0, len( atoms.GetPositions()[ 0 ] ) ):
            center[ 0 ]  =  center[ 0 ] + atoms.GetPositions()[ i ][ 0 ]
            center[ 1 ]  =  center[ 1 ] + atoms.GetPositions()[ i ][ 1 ]
            center[ 2 ]  =  center[ 2 ] + atoms.GetPositions()[ i ][ 2 ]
            count  =  count + 1
    center[ 0 ]  =  center[ 0 ] / count
    center[ 1 ]  =  center[ 1 ] / count
    center[ 2 ]  =  center[ 2 ] / count
    return  center

def  getGyrationRadius( molecule ):
    center  =  getCenter( molecule )
    count  =  0
    result  =  0
    for atoms in molecule.GetConformers():
        for i in range( 0, len( atoms.GetPositions()[ 0 ] ) ):
            dx  =  center[ 0 ] - atoms.GetPositions()[ i ][ 0 ]
            dy  =  center[ 1 ] - atoms.GetPositions()[ i ][ 1 ]
            dz  =  center[ 2 ] - atoms.GetPositions()[ i ][ 2 ]
            result  =  result + ( dx * dx + dy * dy + dz * dz )
            count  =  count + 1
    result  =  result / count
    result  =  np.sqrt( result )
    return  result


tforms = {0: rot_ar_x, 1: rot_ar_y, 2: rot_ar_z}

def  moveToMakeLigandsMeet( partner1, partner2, outputFile ):
    p1  =  Chem.MolFromPDBFile( partner1 )
    p2  =  Chem.MolFromPDBFile( partner2 )
    #for atom in p1.GetAtoms():
    #    print( atom )#.GetName() )
    dictChainIDMol1  =  rdmolops.SplitMolByPDBChainId( p1 )
    dictChainIDMol2  =  rdmolops.SplitMolByPDBChainId( p2 )
    l1  =  dictChainIDMol1[ 'Y' ]
    l2  =  dictChainIDMol2[ 'X' ]
    center1  =  getCenter( l1 )
    center2  =  getCenter( l2 )
    rGyration1  =  getGyrationRadius( l1 )
    rGyration2  =  getGyrationRadius( l2 )
    margin  =  ( rGyration1 + rGyration2 ) #/ 3
    translation  =  [ center1[ 0 ] - center2[ 0 ] + margin, center1[ 1 ] - center2[ 1 ] + margin, center1[ 2 ] - center2[ 2 ] + margin ]
    print( translation )
    rdMolTransforms.TransformConformer( p2.GetConformer(0), tforms[0]( 0, translation[ 0 ] ) )
    rdMolTransforms.TransformConformer( p2.GetConformer(0), tforms[1]( 0, translation[ 1 ] ) )
    rdMolTransforms.TransformConformer( p2.GetConformer(0), tforms[2]( 0, translation[ 2 ] ) )
    # save PDB file:
    print( outputFile )
    success  =  Chem.rdmolfiles.MolToPDBFile( p2, outputFile )
