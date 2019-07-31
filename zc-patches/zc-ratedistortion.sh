#!/bin/bash

if [[ $# < 3 ]]; then
	echo "| Usage - option 1: $0 [data type (-f or -d)] [errBoundMode] [data directory] [extension] [dimension sizes....]"
	echo "|       - option 2: $0 [data type (-f or -d)] [errBoundMode] [varLisetFile]"
	echo "| Example: $0 -f ABS /home/shdi/CESM-testdata/1800x3600 3600 1800"
	exit
fi

datatype=$1
errBoundMode=$2
if [ -d $3 ]; then
	option=1
else
	option=2
fi

if [[ $option == 1 ]]; then
	dataDir=$3
	extension=$4
	dim1=$5
	dim2=$6
	dim3=$7
	dim4=$8
else
	varListFile=$3
fi

#Note: If you run this script by z-checker-installer, SZ_Err_Bounds will be overwritten by ../../errBounds.cfg as follows.
ZZ_Err_Bounds="1E-1 1E-2 1E-3 1E-4"

if [ -f ./errBounds.cfg ]; then
	z_err_env="`cat ./errBounds.cfg | grep -v "#" | grep ZZ_ERR_BOUNDS`"
	echo "export $z_err_env" > env.tmp
	source env.tmp
	rm env.tmp
	ZZ_Err_Bounds="`echo $ZZ_ERR_BOUNDS`"
fi


for errBound in $ZZ_Err_Bounds
do
	if [[ $option == 1 ]]; then
		./test_CompDecomp.sh $datatype $errBound "$dataDir" $extension $dim1 $dim2 $dim3 $dim4
	else
		./test_CompDecomp.sh $datatype $errBound "$varListFile"
	fi
done
