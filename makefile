CC=g++
CFLAGS = -std=c++11 -O3 -Wno-unused-result -DNDEBUG -DFINDPAR
# adjust to path where hdf5 is installed
IFLAGS = -L /usr/local/lib -I /usr/local/include
LFLAGS = -lstdc++ -lz -lhdf5
DEBUGFLAGS = -std=c++11 -g -DFINDPAR

all: findpar


findpar: ./src/main.cpp
	@mkdir -p ./exe
	$(CC) $(CFLAGS) $(IFLAGS) ./src/main.cpp -o ./exe/findpar $(LFLAGS)

nocc: ./src/main.cpp
	@mkdir -p ./exe
	$(CC) $(CFLAGS) -DNOCC $(IFLAGS) ./src/main.cpp -o ./exe/nocc $(LFLAGS)

debug: ./src/main.cpp
	@mkdir -p ./exe
	$(CC) $(DEBUGFLAGS) $(IFLAGS) ./src/main.cpp -o ./exe/debug $(LFLAGS)

clean:
	@rm -rf ./exe/*
