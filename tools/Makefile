CC=gcc
CFLAGS=-std=c99 -Wall -O3 -Wno-unused -g
LDFLAGS= # -lm -lfftw3 -lfftw3f -lz -lpthread
LDFLAGS_EXTRA = -lnetcdf

.SECONDEXPANSION:
%.o: $$(basename $$*).c
	$(CC) $(CFLAGS) $(INCFOLDER) $(DEFINES) -c $< -o $@

dctz-dump: dctz-dump.c
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -f *~ *.o dctz-dump
