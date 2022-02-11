"""
    This reads all the *_lig.pdb files
    and tranfsorms them into the corresponding SDF files
    that will be further pipelined into the
    ROSETTA: molecule_to_params.py
"""

from rdkit import Chem
from rdkit.Chem import rdmolfiles

from pdb_manipulations.getChainID import getListLines
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms

import  os
import  sys
import  shutil


def  makeSDF( inputFile, outputFile ):
    listLines  =  getListLines( inputFile )
    with open( outputFile, "a" ) as out:
        count  =  0
        for line in listLines:
            if "ATOM  " == line[ 0 : 6 ]\
            or "HETATM" == line[ 0 : 6 ]:
                if "H" != line[ 13 ]\
                and "H" != line[ 77 ]\
                and 55 != count:
                    out.write( line + "\n" )
                    count  =  count + 1
    shutil.move( outputFile, inputFile )
    mol  =  Chem.rdmolfiles.MolFromPDBFile( inputFile)
    if None == mol:
        print( f'could not read file {inputFile}' )
    writer = Chem.rdmolfiles.SDWriter( outputFile )
    writer.write( mol )


if "__main__" == __name__:
    if 3 != len( sys.argv ):
        print( "please provide:\n(i). input PDB file;\n(ii). output SDF file;" )
        sys.exit()
    inputFile  =  sys.argv[ 1 ]
    outputFile  =  sys.argv[ 2 ]
    if False == os.path.exists( inputFile ):
        print( f'folder {folderPath} does not exist' )
        sys.exit()
    if True == os.path.exists( outputFile ):
        os.remove( outputFile )
    makeSDF( inputFile, outputFile )

