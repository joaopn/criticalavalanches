CC=g++
CFLAGS = -std=c++11 -O3 -Wno-unused-result -DNDEBUG
# adjust path to find hdf5, for example:
IFLAGS = -L /usr/local/lib -L /usr/lib/x86_64-linux-gnu/hdf5/serial -I /usr/local/include -I /usr/include/hdf5/serial -I ./src/
LFLAGS = -lstdc++ -lz -lhdf5

all: cc nocc cc_orlandi nocc_orlandi

# default (coalesence compensated) local gauss branching network
cc: ./src/main.cpp
	@mkdir -p ./exe
	$(CC) $(CFLAGS) $(IFLAGS) -DTPGAUSS ./src/main.cpp -o ./exe/cc $(LFLAGS)

# disables coalesence compensation
nocc: ./src/main.cpp
	@mkdir -p ./exe
	$(CC) $(CFLAGS) $(IFLAGS) -DNOCC -DTPGAUSS  ./src/main.cpp -o ./exe/nocc $(LFLAGS)

# default (coalesence compensated) local gauss branching network
cc_orlandi: ./src/main.cpp
	@mkdir -p ./exe
	$(CC) $(CFLAGS) $(IFLAGS) -DTPORLANDI ./src/main.cpp -o ./exe/cc_orlandi $(LFLAGS)

# disables coalesence compensation
nocc_orlandi: ./src/main.cpp
	@mkdir -p ./exe
	$(CC) $(CFLAGS) $(IFLAGS) -DNOCC -DTPORLANDI ./src/main.cpp -o ./exe/nocc_orlandi $(LFLAGS)

clean:
	@rm -rf ./exe/*
