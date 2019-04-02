CC=g++
CFLAGS = -std=c++11 -O3 -Wno-unused-result -DNDEBUG
# adjust to path where hdf5 is installed
IFLAGS = -L /usr/local/lib -I /usr/local/include
LFLAGS = -lstdc++ -lz -lhdf5
DEBUGFLAGS = -std=c++11 -g

all: findpar cc nocc

# runs a loop to find drive that matches A=1Hz see main.cpp for arguments
# WITHOUT coalescence compensation
findpar: ./src/main.cpp
	@mkdir -p ./exe
	$(CC) $(CFLAGS) -DFINDPAR -DNOCC $(IFLAGS) ./src/main.cpp -o ./exe/findpar $(LFLAGS)

# default (coalesence compensated) local branching network
cc: ./src/main.cpp
	@mkdir -p ./exe
	$(CC) $(CFLAGS) $(IFLAGS) ./src/main.cpp -o ./exe/cc $(LFLAGS)

# disables coalesence compensation
nocc: ./src/main.cpp
	@mkdir -p ./exe
	$(CC) $(CFLAGS) -DNOCC $(IFLAGS) ./src/main.cpp -o ./exe/nocc $(LFLAGS)

debug: ./src/main.cpp
	@mkdir -p ./exe
	$(CC) $(DEBUGFLAGS) $(IFLAGS) ./src/main.cpp -o ./exe/debug $(LFLAGS)

clean:
	@rm -rf ./exe/*
