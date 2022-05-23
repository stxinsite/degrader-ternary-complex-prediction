from pdb2sql.StructureSimilarity import StructureSimilarity
from collections import namedtuple

CAPRI  =  namedtuple( 'CAPRI', 'fnat lrmsd irmsd category categoryNum dockQ' )



def  getCAPRIPairFiles( modelFile, referenceFile ):
    # create the class instance
    sim = StructureSimilarity( modelFile, referenceFile )

    # compute the irmsd with the two different methods
    irmsd_fast = sim.compute_irmsd_fast(method='svd' )#,izone='capri.izone')
    #irmsd = sim.compute_irmsd_pdb2sql(method='svd',izone='capri.izone')

    # compute the lrmsd with the two different methods
    lrmsd_fast = sim.compute_lrmsd_fast(method='svd' )#,lzone='capri.lzone',check=True)
    #lrmsd = sim.compute_lrmsd_pdb2sql(exportpath=None,method='svd')

    # compute the Fnat with the two different methods
    Fnat_fast = sim.compute_fnat_fast()
    #Fnat = sim.compute_fnat_pdb2sql()

    # compute the DOCKQ
    dockQ = sim.compute_DockQScore(Fnat_fast,lrmsd_fast,irmsd_fast)

    category = sim.compute_CapriClass( Fnat_fast, lrmsd_fast, irmsd_fast )

    capriVal  =  "0"
    if "acceptable" == category:
        capriVal  =  "1"
    elif "medium" == category:
        capriVal  =  "2"
    elif "high" == category:
        capriVal  =  "3"
    capri  =  CAPRI( fnat = Fnat_fast, lrmsd = lrmsd_fast, irmsd = irmsd_fast, category = category, categoryNum = capriVal, dockQ = dockQ )
    return  capri


def  getListLines( file1 ):
    listLines  =  [ line.rstrip( "\t\n" ) for line in open( file1, "r" ) ]
    return  listLines

def  getListAtomNamesPROTACPart( file1, chainIDPROTACPart  =  "X" ):
    listLines  =  getListLines( file1 )
    result  =  []
    for line in listLines:
        if chainIDPROTACPart == line[ 21 ]:
            atomName  =  line[ 12:16 ].strip()
            result.append( atomName )
    return  result


def  getCAPRIPairFilesPROTACTParts( modelFile, referenceFile, chainIDPROTACPart  =  "X" ):
    # create the class instance
    sim = StructureSimilarity( modelFile, referenceFile )

    # compute the irmsd with the two different methods
    irmsd_fast = sim.compute_irmsd_fast(method='svd' )#,izone='capri.izone')
    #irmsd = sim.compute_irmsd_pdb2sql(method='svd',izone='capri.izone')

    # compute the lrmsd with the two different methods
    lrmsd_fast = sim.compute_lrmsd_fast(method='svd' )#,lzone='capri.lzone',check=True)
    #lrmsd = sim.compute_lrmsd_pdb2sql(exportpath=None,method='svd')

    # compute the Fnat with the two different methods
    Fnat_fast = sim.compute_fnat_fast()
    #Fnat = sim.compute_fnat_pdb2sql()

    # compute the DOCKQ
    dockQ = sim.compute_DockQScore(Fnat_fast,lrmsd_fast,irmsd_fast)

    category = sim.compute_CapriClass( Fnat_fast, lrmsd_fast, irmsd_fast )

    capriVal  =  "0"
    if "acceptable" == category:
        capriVal  =  "1"
    elif "medium" == category:
        capriVal  =  "2"
    elif "high" == category:
        capriVal  =  "3"
    capri  =  CAPRI( fnat = Fnat_fast, lrmsd = lrmsd_fast, irmsd = irmsd_fast, category = category, categoryNum = capriVal, dockQ = dockQ )
    return  capri
