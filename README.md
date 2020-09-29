## DCTZ

DCTZ is a lossy compression software for numeric datasets in double (or single) precision. It is based on a block decomposition mechanism combined with well-known discrete cosine transformation (DCT) on floating-point datasets. It also uses an adaptive quantization with two specific task-oriented quantizers: guaranteed user-defined error bounds and higher compression ratios. DCTZ is written in C and is an implementation of the overall idea described in the following paper.

````
Jialing Zhang, Xiaoyan Zhuo, Aekyeung Moon, Hang Liu, Seung Woo Son
"Efficient Encoding and Reconstruction of HPC Datasets for Checkpoint/Restart"
International Conference on Massive Storage Systems and Technology (MSST), May 2019.
````

## Depencencies

DCTZ requires the following libraries to build:
- [FFTW](http://www.fftw.org/): DCT (including inverse DCT) routines in DCTZ are currently based on Fourier transform (DFT) routines in FFTW.
- [zlib](https://www.zlib.net/): compressing bin indices and AC coefficients that need to be saved as it is.
- [Z-checker-installer](https://github.com/CODARcode/z-checker-installer) (optional): when testing DCTZ within Z-Checker. Note that Z-checker-installer will install SZ and zfp automatically. 

## Building 
Retrieve the source code of DCTZ by either git cloning (shown below) or downloading.
````
$ git clone https://github.com/swson/DCTZ
````

Then, make:
```
$ cd DCTZ
$ make
```

The above command will generate `dctz-ec-test` and `dctz-qt-test`. DCTZ has an option to build test programs that can be used in conjunction with Z-checker, a library to characterize the data and check the compression results of lossy compressors. To build that, install Z-checker first and set `ZC_INSTALL_DIR` environment variable:
````
$ export ZC_INSTALL_DIR=path/to/z-checker/install/dir
````

Making after setting `ZC_INSTALL_DIR` will generate two additional test programs, `dctz-ec-zc-test` and `dctz-qt-zc-test`. 

## Testing

To test the compiled DCTZ is working properly, download the test datasets from [here](https://sites.uml.edu/seungwoo-son/files/2019/07/dctz-test-data.zip) and run the following command:
````
$ cd tests
$ wget https://sites.uml.edu/seungwoo-son/files/2019/07/dctz-test-data.zip
$ unzip dctz-test-data.zip
$ ./test-dctz.sh
````

## Testing with Z-checker
To test DCTZ with the state-of-the-art lossy compressors, SZ and zfp, use the Z-checker framework developed by Argonne National Laboratory.
1. Install Z-checker-installer, which can be cloned or downloaded from
[Z-checker-installer](http://github.com/CODARcode/z-checker-installer).

After installing z-checker-installer, add the installed zfp's library path to LD_LIBRARY_PATH.
```
$ export LD_LIBRARY_PATH=${ZC_INSTALL_HOME}/zfp/zfp-install/lib:$LD_LIBRARY_PATH
```

2. Copy the compiled DCTZ test binaries to `{ZC_INSTALL_HOME}/DCTZ`.
```
$ mkdir ${ZC_INSTALL_HOME}/DCTZ
$ cp dctz-ec-zc-test ${ZC_INSTALL_HOME}/DCTZ/.
$ cp dctz-qt-zc-test ${ZC_INSTALL_HOME}/DCTZ/.
```

Note that `$ZC_INSTALL_HOME` is not same as `$ZC_INSTALL_DIR`. Typically, `$ZC_INSTALL_DIR` is equal to `$ZC_INSTALL_HOME/Z-checker/zc-install`.

3. Add the DCTZ configuration files to Z-Checker.
```
$ cp zc-patches/manageCompressor-dctz-ec.cfg ${ZC_INSTALL_HOME}/.
$ cp zc-patches/manageCompressor-dctz-qt.cfg ${ZC_INSTALL_HOME}/.
```

Modify the `${ZC_INSTALL_HOME}` variable in `manageCompressor-dctz-ec.cfg` and `manageCompressor-dctz-qt.cfg`, and execute the following commands:
```
$ ./manageCompressor -a DCTZ -m ec -c manageCompressor-dctz-ec.cfg
$ ./manageCompressor -a DCTZ -m qt -c manageCompressor-dctz-qt.cfg
```

Then, open `errBounds.cfg` under `${ZC_INSTALL_HOME}` and modify the error bounds (e.g., 1E-3, 1E-4, etc.) and compression cases (`dctz_ec` and `dctz_qt`) for DCTZ. Refer to an example `erroBounds.cfg` under the `zc-patches` folder.


4. Create a new test-case, by executing `createZCCase.sh [test-case-name]`.


5. Modify listing files for all the test data.
```
$ cp zc-patches/varInfo.txt ${ZC_INSTALL_HOME}/.
```

Change file paths in `varInfo.txt` to use test data (`$ cd tests, $ unzip dctz-test-data.zip`).


6. Run Z-checker.
```
$ ./runZCCase.sh -d REL [test-case-name] varInfo.txt
```

