#!/bin/bash

curDir=$(pwd)
if [[ -d "OUTPUT_INVERSION" ]]; then
id=$(ls output_*.log | grep -oe "[0-9]\+" | head -n 1)
echo "Creation of directory OUTPUT_INVERSION_"$id
mv output_*.log OUTPUT_INVERSION/
mv esl*.err OUTPUT_INVERSION/
mv OUTPUT_INVERSION OUTPUT_INVERSION_$id
#tar -czvf OUTPUT_INVERSION_$(id).tar.gz OUTPUT_INVERSION_$id
fi

