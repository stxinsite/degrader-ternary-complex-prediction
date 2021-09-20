#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF

ln -sv $WEST_SIM_ROOT/common_files/{complex.parm7,rest.in} .

if [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_CONTINUES" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/md.in > md.in
  ln -sv $WEST_PARENT_DATA_REF/seg.rst ./parent.rst
elif [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_NEWTRAJ" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/md.in > md.in
  ln -sv $WEST_PARENT_DATA_REF ./parent.rst
fi


export CUDA_DEVICES=(`echo $CUDA_VISIBLE_DEVICES_ALLOCATED | tr , ' '`)
export CUDA_VISIBLE_DEVICES=${CUDA_DEVICES[$WM_PROCESS_INDEX]}

echo "RUNSEG.SH: CUDA_VISIBLE_DEVICES_ALLOCATED = " $CUDA_VISIBLE_DEVICES_ALLOCATED
echo "RUNSEG.SH: WM_PROCESS_INDEX = " $WM_PROCESS_INDEX
echo "RUNSEG.SH: CUDA_VISIBLE_DEVICES = " $CUDA_VISIBLE_DEVICES


$PMEMD -O -i md.in   -p complex.parm7  -c parent.rst -r seg.rst -x seg.nc  -o seg.log  -inf seg.nfo

COMMAND0="parm $WEST_SIM_ROOT/common_files/complex.parm7\n"
COMMAND0="${COMMAND0} trajin parent.rst \n"
COMMAND0="${COMMAND0} strip :WAT:K+:Cl- \n"
COMMAND0="${COMMAND0} trajout parent_solutes.rst"
COMMAND0="${COMMAND0} go\n"
sleep 1
echo -e ${COMMAND0} | $CPPTRAJ
sleep 1

COMMAND1="parm $WEST_SIM_ROOT/common_files/complex.parm7\n"
COMMAND1="${COMMAND1} trajin seg.nc \n"
COMMAND1="${COMMAND1} strip :WAT:K+:Cl- \n"
COMMAND1="${COMMAND1} trajout solutes.rst"
COMMAND1="${COMMAND1} go\n"
sleep 1
echo -e ${COMMAND1} | $CPPTRAJ
sleep 1

COMMAND2="parm $WEST_SIM_ROOT/common_files/reference.parm7\n"
COMMAND2="${COMMAND2} trajin parent_solutes.rst\n"
COMMAND2="${COMMAND2} trajin solutes.rst\n"
COMMAND2="${COMMAND2} autoimage \n"
COMMAND2="${COMMAND2} reference $WEST_SIM_ROOT/common_files/reference.rst7 \n"
COMMAND2="${COMMAND2} nativecontacts name SVI  :1-115&!@H* :116-264&!@H*  byresidue mindist distance 4.5  reference  out Smarca_VHL.dat \n"
COMMAND2="${COMMAND2} nativecontacts name SPI  :1-115&!@H* :FWZ&!@H*      byresidue mindist distance 4.5  reference  out Smarca_Protac.dat \n"
COMMAND2="${COMMAND2} rms iRMS_Smarca  '@4418,4419,4420,4421,4422,4423,4424,4425,4426,4427,4428,4429,4430,4431,4432,4433,4434,4435,4436,4437 | :33,34,36,37,38,43,46,54,55,81,84,85,88,89,95&!@H*'  reference  out iRMSD.dat \n"
COMMAND2="${COMMAND2} go\n"
#sleep 1
#echo -e $COMMAND2 | $CPPTRAJ
#sleep 1
# This wasn't working at one random iteration, thus we want to monitor

sleep 2
echo -e $COMMAND2 | $CPPTRAJ > cpptraj1.out
if [ ! -f "Smarca_VHL.dat" ] || [ ! -f "Smarca_Protac.dat" ] || [ ! -f "iRMSD.dat" ] ; then
sleep 2
echo -e $COMMAND2 | $CPPTRAJ > cpptraj2.out
fi

python $WEST_SIM_ROOT/common_files/combine_cvs.py Smarca_VHL.dat  Smarca_Protac.dat  CV1.dat

sleep 1

paste <(cat CV1.dat) <(cat iRMSD.dat | tail -n 2 | gawk '{print $2}')  >$WEST_PCOORD_RETURN
paste <(cat CV1.dat) <(cat iRMSD.dat | tail -n 2 | gawk '{print $2}')  >pcoord.dat

# Clean up
rm -f md.in rest.in parent.rst complex.parm7 seg.log seg.nfo  parent_solutes.rst mden 

