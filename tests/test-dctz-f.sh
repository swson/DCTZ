#!/bin/bash
echo "deleting *.z*, *.r, ..."
rm -f *.z* *.r 

#wget http://www.mcs.anl.gov/~shdi/download/CESM-ATM-tylor.tar.gz
#tar xvfz CESM-ATM-tylor.tar.gz

for i in '1E-3' '1E-4' '1E-5'
do
    for f in 'CESM-ATM-tylor/1800x3600/CLDHGH_1_1800_3600.dat 3600 1800' 'CESM-ATM-tylor/1800x3600/CLDLOW_1_1800_3600.dat 3600 1800' 'CESM-ATM-tylor/1800x3600/FLDSC_1_1800_3600.dat 3600 1800' 'CESM-ATM-tylor/1800x3600/FREQSH_1_1800_3600.dat 3600 1800' 'CESM-ATM-tylor/1800x3600/PHIS_1_1800_3600.dat 3600 1800'
    do
	X=`echo $f|cut -c1-4`
	echo running dctz-ec-test $f $i
	../dctz-ec-test -f $i var $f 2>&1 | tee dctz-ec-comp.$X.$i.log
	echo running dctz-qt-test $f $i
	../dctz-qt-test -f $i var $f 2>&1 | tee dctz-qt-comp.$X.$i.log
    done
done
