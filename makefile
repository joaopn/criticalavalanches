CC=g++
CFLAGS = -std=c++11 -O3 -Wno-unused-result
# adjust path to find hdf5, for example:
IFLAGS = -L /usr/local/lib -L /usr/lib/x86_64-linux-gnu/hdf5/serial -I /usr/local/include -I /usr/include/hdf5/serial -I ./src/
LFLAGS = -Wl,-rpath,$(LIBRARY_PATH) -lstdc++ -lz -lhdf5 $(LDFLAGS)

all: cc cc_orlandi cc_random

more: nocc nocc_orlandi nocc_random

# default (coalesence compensated) local gauss branching network
cc: ./src/main.cpp
	@mkdir -p ./exe
	$(CC) $(CFLAGS) $(IFLAGS) -DTPGAUSS ./src/main.cpp -o ./exe/cc $(LFLAGS)

# disables coalesence compensation
nocc: ./src/main.cpp
	@mkdir -p ./exe
	$(CC) $(CFLAGS) $(IFLAGS) -DNOCC -DTPGAUSS  ./src/main.cpp -o ./exe/nocc $(LFLAGS)

# orlandi topology
cc_orlandi: ./src/main.cpp
	@mkdir -p ./exe
	$(CC) $(CFLAGS) $(IFLAGS) -DTPORLANDI ./src/main.cpp -o ./exe/cc_orlandi $(LFLAGS)

nocc_orlandi: ./src/main.cpp
	@mkdir -p ./exe
	$(CC) $(CFLAGS) $(IFLAGS) -DNOCC -DTPORLANDI ./src/main.cpp -o ./exe/nocc_orlandi $(LFLAGS)

# random graph topology.
cc_random: ./src/main.cpp
	@mkdir -p ./exe
	$(CC) $(CFLAGS) $(IFLAGS) -DTPRANDOM ./src/main.cpp -o ./exe/cc_random $(LFLAGS)

nocc_random: ./src/main.cpp
	@mkdir -p ./exe
	$(CC) $(CFLAGS) $(IFLAGS) -DNOCC -DTPRANDOM ./src/main.cpp -o ./exe/nocc_random $(LFLAGS)

clean:
	@rm -rf ./exe/*
