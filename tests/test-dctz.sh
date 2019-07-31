#!/bin/bash
echo "deleting *.z*, *.r, ..."
rm -f *.z* *.r 

for i in '1E-3' '1E-4' '1E-5'
do
    for f in 'sedov-chk139-dens.bin 31040' 'cellular-0249.bin 32768' 'rlds.bin 12960' 'mrsos.bin 12960' 'eddy-chk50-pres.bin 16384' 'vortex-chk50-pres.bin 37024'
    do
	X=`echo $f|cut -c1-4`
	echo running dctz-ec-test $f $i
	../dctz-ec-test -d $i var $f 2>&1 | tee dctz-ec-comp.$X.$i.log
	echo running dctz-qt-test $f $i
	../dctz-qt-test -d $i var $f 2>&1 | tee dctz-qt-comp.$X.$i.log
    done
done
