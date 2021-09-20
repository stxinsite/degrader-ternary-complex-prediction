#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT/init_files

COMMAND1="parm $WEST_SIM_ROOT/common_files/complex.parm7\n"
COMMAND1="${COMMAND1} trajin $WEST_STRUCT_DATA_REF \n"
COMMAND1="${COMMAND1} strip :WAT:K+:Cl- \n"
COMMAND1="${COMMAND1} trajout solutes.rst"
COMMAND1="${COMMAND1} go\n"
echo -e ${COMMAND1} | $CPPTRAJ

COMMAND2="parm $WEST_SIM_ROOT/common_files/reference.parm7\n"
COMMAND2="${COMMAND2} trajin solutes.rst\n"
COMMAND2="${COMMAND2} reference $WEST_SIM_ROOT/common_files/reference.rst7 \n"
COMMAND2="${COMMAND2} nativecontacts name SVI  :1-115&!@H* :116-264&!@H*  byresidue mindist distance 4.5  reference  out Smarca_VHL.dat \n"
COMMAND2="${COMMAND2} nativecontacts name SPI  :1-115&!@H* :FWZ&!@H*      byresidue mindist distance 4.5  reference  out Smarca_Protac.dat \n"
COMMAND2="${COMMAND2} rms iRMS_Smarca  '@4418,4419,4420,4421,4422,4423,4424,4425,4426,4427,4428,4429,4430,4431,4432,4433,4434,4435,4436,4437 | :33,34,36,37,38,43,46,54,55,81,84,85,88,89,95&!@H*'  reference  out iRMSD.dat \n"
COMMAND2="${COMMAND2} go\n"
echo -e ${COMMAND2} | $CPPTRAJ

python $WEST_SIM_ROOT/common_files/combine_cvs.py Smarca_VHL.dat  Smarca_Protac.dat  CV1.dat

sleep 1  # THIS IS CRUCIAL;  OTHERWISE YOU GET A FORMATTING ERROR (AND DEBUGGING IS A PAIN)

paste <(cat CV1.dat) <(cat iRMSD.dat | tail -n 1 | gawk '{print $2}')  >>$WEST_PCOORD_RETURN
paste <(cat CV1.dat) <(cat iRMSD.dat | tail -n 1 | gawk '{print $2}')  >>pcoord.dat

if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi
