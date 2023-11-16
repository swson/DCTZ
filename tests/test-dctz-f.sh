#!/bin/bash
echo "deleting *.z*, *.r, ..."
rm -f *.z* *.r dct_result.bin

FILE=dct_result.bin
BIN_FILE=bin_index.bin

#wget http://www.mcs.anl.gov/~shdi/download/CESM-ATM-tylor.tar.gz
#tar xvfz CESM-ATM-tylor.tar.gz

while IFS="" read -r f || [ -n "$f" ]
do
    for i in '1E-3' '1E-4' '1E-5'
    do
	X=`echo $f|cut -c26-29`
	echo running dctz-ec-test $f $i
	../dctz-ec-test -f $i var $f 2>&1 | tee dctz-ec-comp.$X.$i.log
	if [ -f "$FILE" ]; then
	    echo "renaming dct_result.bin..."
	    mv dct_result.bin dct_result.f.ec.$X.$i.bin
	else
	    echo "$FILE does not exist"
	fi
	if [ -f "$BIN_FILE" ]; then
	    echo "renaming bin_file.bin..."
	    mv bin_index.bin bin_index.f.ec.$X.$i.bin
	else
	    echo "$FILE does not exist"
	fi

	echo ""

	echo running dctz-qt-test $f $i
	../dctz-qt-test -f $i var $f 2>&1 | tee dctz-qt-comp.$X.$i.log
	if [ -f "$FILE" ]; then
	    echo "renaming dct_result.bin..."
	    mv dct_result.bin dct_result.f.qt.$X.$i.bin
	else
	    echo "$FILE does not exist"
	fi
	if [ -f "$BIN_FILE" ]; then
	    echo "renaming bin_file.bin..."
	    mv bin_index.bin bin_index.f.qt.$X.$i.bin
	else
	    echo "$FILE does not exist"
	fi

	echo ""
    done
done < list-CESM-ATM-tylor.txt
