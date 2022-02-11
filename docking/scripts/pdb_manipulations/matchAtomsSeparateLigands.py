"""
    This module is used to prepare PDB files to calculate CAPRI criteria:
    First, the atoms in the reference structure and the model structure are matched;
    Second, the ligands, warheads and linkers are assigned specific chain IDs;
"""

CONST_CHAIN_ID_PROTEIN_1__  =  "A"
CONST_CHAIN_ID_PROTEIN_2__  =  "C"
CONST_CHAIN_ID_LIGAND_1__   =  "X"
CONST_CHAIN_ID_LIGAND_2__   =  "Y"
CONST_CHAIN_ID_LINKER__     =  "D"


def  getListLines( file1 ):
    listLines  =  [ line.rstrip( "\t\n" ) for line in open( file1, "r" ) ]
    return  listLines


def  matchAtomsProteinsOnly( referenceFile, modelFile, referenceFileOut, modelFileOut ):
    listLinesReference  =  getListLines( referenceFile )
    listLinesModel      =  getListLines( modelFile )
    print( "matchAtomsProteinsOnly: modelFile = ", modelFile )
    print( "matchAtomsProteinsOnly: referenceFile = ", referenceFile )
    with open( referenceFileOut, "a" ) as outRef:
        with open( modelFileOut, "a" ) as outModel:
            j  =  0
            for i in range( 0, len( listLinesReference ) ):
                lineRef  =  listLinesReference[ i ]
                if "TER" == lineRef[ 0 : 3 ]\
                or "ATOM  " != lineRef[ 0 : 6 ]:
                    continue
                jTmp  =  j
                for k in range( 0, 10 ):
                    if len( listLinesModel ) > j + k:
                        lineModel  =  listLinesModel[ j + k ]
                        if 21 >= len( lineModel ):
                            j  =  j + 1
                            continue
                        atomNameRef  =  lineRef[ 12:17 ].strip()
                        resNameRef   =  lineRef[ 17:21 ].strip()
                        chainIDRef   =  lineRef[ 21 ]
                        atomNameMod  =  lineModel[ 12:17 ].strip()
                        resNameMod   =  lineModel[ 17:21 ].strip()
                        chainIDMod   =  lineModel[ 21 ]
                        rsnRef  =  lineRef[ 22 : 26 ].strip()
                        rsnModel  =  lineModel[ 22 : 26 ].strip()
                        if atomNameRef == atomNameMod\
                        and chainIDRef == chainIDMod:
                            if CONST_CHAIN_ID_PROTEIN_1__ == chainIDRef\
                            or CONST_CHAIN_ID_PROTEIN_2__ == chainIDRef:
                                if int( rsnRef ) != int( rsnModel ):
                                    outRef.write( lineRef[ 0:22 ] + rsnModel.rjust(4) + lineRef[ 26 : len( lineRef ) ] + "\n" )
                                else:
                                    outRef.write( lineRef + "\n" )
                                outModel.write( lineModel + "\n" )
                                j  =  j + k + 1
                                break
                            else:
                                j  =  j + 1
                                break



def  matchAtomsSeparateLigands( referenceFile, modelFile, referenceFileOut, modelFileOut ):
    listLinesReference  =  getListLines( referenceFile )
    listLinesModel      =  getListLines( modelFile )
    print( "mathAtomsSeparateLigands: modelFile = ", modelFile )
    print( "mathAtomsSeparateLigands: referenceFile = ", referenceFile )
    with open( referenceFileOut, "a" ) as outRef:
        with open( modelFileOut, "a" ) as outModel:
            j  =  0
            linker  =  False
            for i in range( 0, len( listLinesReference ) ):
                lineRef  =  listLinesReference[ i ]
                #print( "ref = ", lineRef[ 0 : 6 ] )
                if "TER" == lineRef[ 0 : 3 ]\
                or "ATOM  " != lineRef[ 0 : 6 ]:
                    continue
                #print( "ref", lineRef )
                jTmp  =  j
                for k in range( 0, 10 ):
                    if len( listLinesModel ) > j + k:
                        lineModel  =  listLinesModel[ j + k ]
                        #print( "model", lineModel )
                        #print( "model = ", lineModel[ 0 : 6 ] )
                        #if "TER" == lineModel[ 0 : 3 ]\
                        #or "ATOM  " != lineModel[ 0 : 6 ]:
                        #    continue
                        if 21 >= len( lineModel ):
                            j  =  j + 1
                            continue
                        atomNameRef  =  lineRef[ 12:17 ].strip()
                        resNameRef   =  lineRef[ 17:21 ].strip()
                        chainIDRef   =  lineRef[ 21 ]
                        atomNameMod  =  lineModel[ 12:17 ].strip()
                        resNameMod   =  lineModel[ 17:21 ].strip()
                        chainIDMod   =  lineModel[ 21 ]
                        #print( atomNameRef, atomNameMod )
                        #print( resNameRef, resNameMod )
                        #print( chainIDRef, chainIDMod )
                        rsnRef  =  lineRef[ 22 : 26 ].strip()
                        rsnModel  =  lineModel[ 22 : 26 ].strip()
                        if atomNameRef == atomNameMod\
                        and chainIDRef == chainIDMod:
                            if CONST_CHAIN_ID_LIGAND_1__ == chainIDRef\
                            and "LG1" == resNameMod:
                                if int( rsnRef ) != int( rsnModel ):
                                    outRef.write( lineRef[ 0:17 ] + "LG1" + lineRef[ 20 : 22 ] + rsnModel.rjust(4) + lineRef[ 26 : len( lineRef ) ] + "\n" )
                                else:
                                    outRef.write( lineRef[ 0:17 ] + "LG1" + lineRef[ 20 : len( lineRef ) ] + "\n" )
                            else:
                                if int( rsnRef ) != int( rsnModel ):
                                    outRef.write( lineRef[ 0:22 ] + rsnModel.rjust(4) + lineRef[ 26 : len( lineRef ) ] + "\n" )
                                else:
                                    outRef.write( lineRef + "\n" )
                            outModel.write( lineModel + "\n" )
                            j  =  j + k + 1
                            break
                        elif atomNameRef == atomNameMod\
                        and CONST_CHAIN_ID_LIGAND_2__ == chainIDRef\
                        and CONST_CHAIN_ID_LIGAND_1__ == chainIDMod:
                            if int( rsnRef ) != int( rsnModel ):
                                outRef.write( lineRef[ 0:17 ] + "LG2" + lineRef[ 20 : 22 ] + rsnModel.rjust(4) + lineRef[ 26 : len( lineRef ) ] + "\n" )
                            else:
                                outRef.write( lineRef[ 0:17 ] + "LG2" + lineRef[ 20 : len( lineRef ) ] + "\n" )
                            outModel.write( lineModel[ 0:17 ] + "LG2 Y" + lineModel[ 22 : len( lineModel ) ] + "\n" )
                            j  =  j + k + 1
                            if i + 1 < len( listLinesReference ):
                                if 22 <= len( listLinesReference[ i + 1 ] )\
                                and CONST_CHAIN_ID_LINKER__ == listLinesReference[ i + 1 ][ 21 ]:
                                    linker  =  True
                            break
                        elif True == linker\
                        and atomNameRef == atomNameMod\
                        and CONST_CHAIN_ID_LINKER__ == chainIDRef\
                        and CONST_CHAIN_ID_LIGAND_1__ == chainIDMod:
                            if int( rsnRef ) != int( rsnModel ):
                                outRef.write( lineRef[ 0:17 ] + "LIN Z" + rsnModel.rjust(4) + lineRef[ 26 : len( lineRef ) ] + "\n" )
                            else:
                                outRef.write( lineRef[ 0:17 ] + "LIN Z" + lineRef[ 22 : len( lineRef ) ] + "\n" )
                            outModel.write( lineModel[ 0:17 ] + "LIN Z" + lineModel[ 22 : len( lineModel ) ] + "\n" )
                            j  =  j + k + 1
                            break
                #if jTmp == j:
                #    j  =  j + 1
            #print( "j = ", j )


def  matchAtomsProteinLigand( referenceFile, modelFile, referenceFileOut, modelFileOut, chainIDProtein, chainIDLigand ):
    listLinesReference  =  getListLines( referenceFile )
    listLinesModel      =  getListLines( modelFile )
    print( "matchAtomsProteinLigand: modelFile = ", modelFile )
    print( "matchAtomsProteinLigand: referenceFile = ", referenceFile )
    with open( referenceFileOut, "a" ) as outRef:
        with open( modelFileOut, "a" ) as outModel:
            j  =  0
            for i in range( 0, len( listLinesReference ) ):
                lineRef  =  listLinesReference[ i ]
                if "TER" == lineRef[ 0 : 3 ]\
                or "ATOM  " != lineRef[ 0 : 6 ]:
                    continue
                for k in range( 0, 10 ):
                    if len( listLinesModel ) > j + k:
                        lineModel  =  listLinesModel[ j + k ]
                        if 21 >= len( lineModel ):
                            j  =  j + 1
                            continue
                        atomNameRef  =  lineRef[ 12:17 ].strip()
                        resNameRef   =  lineRef[ 17:21 ].strip()
                        chainIDRef   =  lineRef[ 21 ]
                        atomNameMod  =  lineModel[ 12:17 ].strip()
                        resNameMod   =  lineModel[ 17:21 ].strip()
                        chainIDMod   =  lineModel[ 21 ]
                        rsnRef  =  lineRef[ 22 : 26 ].strip()
                        rsnModel  =  lineModel[ 22 : 26 ].strip()
                        if atomNameRef == atomNameMod\
                        and chainIDRef == chainIDMod:
                            if chainIDProtein == chainIDRef\
                            or chainIDLigand == chainIDRef:
                                if int( rsnRef ) != int( rsnModel ):
                                    if chainIDLigand == chainIDRef:
                                        outRef.write( lineRef[ 0:17 ] + "LG1" + lineRef[ 20:22 ] + rsnModel.rjust(4) + lineRef[ 26 : len( lineRef ) ] + "\n" )
                                    else:
                                        outRef.write( lineRef[ 0:22 ] + rsnModel.rjust(4) + lineRef[ 26 : len( lineRef ) ] + "\n" )
                                else:
                                    if chainIDLigand == chainIDRef:
                                        outRef.write( lineRef[ 0:17 ] + "LG1" + lineRef[ 20:len( lineRef ) ] + "\n" )
                                    else:
                                        outRef.write( lineRef + "\n" )
                                if chainIDLigand == chainIDRef:
                                    outModel.write( lineModel[ 0:17 ] + "LG1" + lineModel[ 20:len( lineModel ) ] + "\n" )
                                else:
                                    outModel.write( lineModel + "\n" )
                                j  =  j + k + 1
                                break
                            else:
                                j  =  j + 1
                                break

