#!/bin/bash
echo "deleting *.z*, *.r, ..."
rm -f *.z* *.r dct_result.bin

DCT_FILE=dct_result.bin
BIN_FILE=bin_index.bin
QTABLE_FILE=qtable.bin
AC_EXACT_FILE=AC_exact.bin
DC_FILE=DC.bin

# wget https://sites.uml.edu/seungwoo-son/files/2019/07/dctz-test-data.zip

while IFS="" read -r f || [ -n "$f" ]
do
    for i in '1E-3' '1E-4' '1E-5'
    do
	X=`echo $f|cut -c1-4`
#	echo running dctz-ec-test $f $i
#	../dctz-ec-test -d $i var $f 2>&1 | tee dctz-ec-comp.$X.$i.log
#	if [ -f "$DCT_FILE" ]; then
#	    echo "renaming dct_result.bin..."
#	    mv dct_result.bin dct_result.ec.$X.$i.bin
#	fi
#	if [ -f "$BIN_FILE" ]; then
#	    echo "renaming bin_file.bin..."
#	    mv bin_index.bin bin_index.ec.$X.$i.bin
#	fi
	
#	echo ""
	
	echo running dctz-qt-test $f $i
	../dctz-qt-test -d $i var $f 2>&1 | tee dctz-qt-comp.$X.$i.log
	if [ -f "$DCT_FILE" ]; then
	    echo "renaming dct_result.bin..."
	    mv dct_result.bin dct_result.qt.$X.$i.bin
	fi
	if [ -f "$BIN_FILE" ]; then
	    echo "renaming bin_file.bin..."
	    mv bin_index.bin bin_index.qt.$X.$i.bin
	fi
	if [ -f "$QTABLE_FILE" ]; then
	    echo "renaming qtable.bin..."
	    mv qtable.bin qtable.qt.$X.$i.bin
	fi
	if [ -f "$AC_EXACT_FILE" ]; then
	    echo "renaming AC_exact.bin..."
	    mv AC_exact.bin AC_exact.qt.$X.$i.bin
	fi
	if [ -f "$DC_FILE" ]; then
	    echo "renaming DC.bin..."
	    mv DC.bin DC.qt.$X.$i.bin
	fi

	echo ""
    done
done < list-msst19.txt
