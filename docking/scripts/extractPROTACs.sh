#!/bin/bash

path=$1
pathOut=$2

for file in ${path}/*.pdb; do
  echo $file
  fileNoH=${pathOut}/$(basename "$file")".noH.pdb"
  python ${BASE_FOLDER}/scripts/removeH.py $file $fileNoH
  fileOut=$fileNoH"_protac.pdb"
  grep 'LG1' $fileNoH > $fileOut
  fileSDF=$fileOut".sdf"
  obabel -ipdb $fileOut -O $fileSDF
done
  
