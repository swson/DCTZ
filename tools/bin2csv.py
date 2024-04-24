#!/usr/bin/env python3
import sys
import numpy as np

def bin2csv(filename, dtype):
    with open(filename, 'rb') as fid:
        if dtype == 'uint8':
            input_data = np.fromfile(fid, dtype=np.uint8)
        elif dtype == 'float':
            input_data = np.fromfile(fid, dtype=np.float32)
        elif dtype == 'double':
            input_data = np.fromfile(fid, dtype=np.float64)
    token = filename.split('.')[0]
    #print(token)
    if (dtype == 'uint8'):
        np.savetxt(f"{token}.csv", input_data, delimiter=",", fmt="%d")
    elif (dtype == 'float' or dtype == 'double'):
        np.savetxt(f"{token}.csv", input_data, delimiter=",", fmt="%e")

def main():
    if len(sys.argv) != 3 :
        #print("Usage: {sys.argv[0]} filename.bin dtype")
        sys.exit("Usage: %s filename.bin dtype={uint8,double}" % sys.argv[0])

    if sys.argv[2] not in ['uint8', 'float', 'double']:
        sys.exit("%s: unsupported data type" % sys.argv[2])
    bin2csv(sys.argv[1], sys.argv[2]);

if __name__ == "__main__":
    main()
