"""
    parses the output of RosettaDock: docking_protocol file score.sc
    ranks by the specified parameter and outputs a text file
    with filenames of resulting top-N structures
"""

import  os
import  sys
from collections import namedtuple
from  pdb_manipulations.getChainID import getListLines


ScoreEntry  =  namedtuple( 'ScoreEntry', 'totalScore_ rms_ CAPRIrank_ Fnat_ I_sc_ Irms_ filename_' )


def  main( inputFile, outputFile, outputPDBs, pathPDBFiles ):
    listLinesInput  =  getListLines( inputFile )
    listOutput  =  []
    for line in listLinesInput:
        if False == ( "SEQUENCE:" in line )\
        and False == ( "total_score" in line ):
            totalScore  =  line[ 8:20 ].strip()
            rms  =  line[ 23 :31 ].strip()
            CAPRI_rank  =  line[ 34:44 ].strip()
            Fnat  =  line[ 45:54 ].strip()
            I_sc  =  ( line[ 55:66 ].strip() )
            Irms  =  line[ 67:75 ].strip()
            filename  =  line.split( " " )[ len( line.split( " " ) ) - 1 ]
            entry  =  ScoreEntry( totalScore_ = totalScore, rms_ = rms, CAPRIrank_ = CAPRI_rank, Fnat_ = Fnat, I_sc_ = I_sc, Irms_ = Irms, filename_ = filename )
            listOutput.append( entry )
    listOutputSorted  =  sorted( listOutput, key = lambda x: x.I_sc_ )
    with open( outputFile, "a" ) as out:
        with open( outputPDBs, "a" ) as outPDB:
            for i in range( 0, len( listOutputSorted ) ):
                entry  =  listOutputSorted[ i ]
                if 0 == i:
                    lineToSave = "total_score,lrmsd,category,fnat,I_sc,irmsd,filename\n"
                    out.write(lineToSave)
                lineToSave  =  entry.totalScore_ + "," + entry.rms_ + "," + entry.CAPRIrank_ + "," + entry.Fnat_ + "," + str( entry.I_sc_ ) + "," + entry.Irms_ + "," + entry.filename_ + "\n"
                out.write( lineToSave )
                filePath  =  os.path.join( pathPDBFiles, entry.filename_ + ".pdb" )
                if False == os.path.exists( filePath ):
                    print( f'file {filePath} does not exist' )
                    continue#return
                lineToSave  =  filePath + "\n"
                outPDB.write( lineToSave )


if "__main__" == __name__:
    if 5 != len( sys.argv ):
        print( "please provide:\n(i). input score.sc file;\n(ii). output score.sc file to save results;\n(iii). output file for list of PDB files from Rosetta;\n(iv). path to folder where the PDB files are located;\n" )
        exit()
    inputFile  =  sys.argv[ 1 ]
    outputFile  =  sys.argv[ 2 ]
    outputPDBs  =  sys.argv[ 3 ]
    pathPDBFiles  =  sys.argv[ 4 ]
    if False == os.path.exists( inputFile ):
        print( f'file {inputFile} does not exist' )
        exit()
    if False == os.path.exists( pathPDBFiles )\
    or False == os.path.isdir( pathPDBFiles ):
        print( f'folder {pathPDBFiles} does not exist' )
        exit()
    if True == os.path.exists( outputFile ):
        os.remove( outputFile )
    if True == os.path.exists( outputPDBs ):
        os.remove( outputPDBs )
    main( inputFile, outputFile, outputPDBs, pathPDBFiles )
