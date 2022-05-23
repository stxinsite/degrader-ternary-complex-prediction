"""
    employ RDKit to get RMSD value from two aligned structures.
    assumption: atom names and ordering in both PDB files must match.
"""

import  os
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from PROTAC_ternary import PDBContainer
#from pdb_manipulations.getChainID import getListLines


def  getAtomMap( pdb1 ):
    listLines  =  getListLines( pdb1 )
    atom_ids  =  {}
    listAtomIDs  =  []
    result  =  []
    i  =  0
    for line in listLines:
        if "HETATM" == line[ 0 : 6 ]:
            atomName  =  line[ 11 : 16 ].strip()
            #residueName  =  line[ 17 : 20 ].strip()
            atom_ids[ atomName ]  =  i
            i  =  i + 1
            listAtomIDs.append( atom_ids )
            result.append( ( atom_ids, atom_ids ) )
            atom_ids  =  {}
    #print( listAtomIDs )
    #result  =  list( ( listAtomIDs, listAtomIDs ) )
    return  result



def  getRMSD( pdb1, pdb2, decoy_atoms, linker_atoms, aln_iters ):
    if False == os.path.exists( pdb1 ):
        print( f'file {pdb1} does not exist' )
        exit()
    if False == os.path.exists( pdb2 ):
        print( f'file {pdb2} does not exist' )
        exit()
    atom_map  =  getAtomMap( pdb1 )
    struct1  =  Chem.rdmolfiles.MolFromPDBFile( pdb1,
                                                sanitize = True,
                                                removeHs = True,
                                                flavor = 0,
                                                proximityBonding = False
                                                )
    struct2  =  Chem.rdmolfiles.MolFromPDBFile( pdb2,
                                                sanitize = True,
                                                removeHs = True,
                                                flavor = 0,
                                                proximityBonding = False
                                                )
    #print( pdb1, pdb2 )
#    atom_map  =  getAtomMap( pdb1 )
    print( atom_map )
#    rmsd = Chem.rdMolAlign.AlignMol( struct1,
#                                    struct2,
#                                    atomMap = atom_map)
    decoy  =  PDBContainer( pdb1 )
    linker  =  PDBContainer( pdb2 )
    rmsd  =  linker.align( decoy, targ_atoms=linker_atoms, ref_atoms=decoy_atoms, aln_iters=aln_iters )
    return  rmsd
#                                        atomMap=atom_map,
#                                        maxIters=aln_iters)

