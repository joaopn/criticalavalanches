CC=g++
CFLAGS = -std=c++11 -O3 -Wno-unused-result -DNDEBUG #-fdiagnostics-color=always
LFLAGS = -lstdc++ -lz
DEBUGFLAGS = -std=c++11 -g

all: ./src/gh_test_src.cpp
	@mkdir -p ./exe
	$(CC) $(CFLAGS) ./src/gh_test_src.cpp -o ./exe/gh_test_src $(LFLAGS)

debug: ./src/gh_test_src.cpp
	@mkdir -p ./exe
	$(CC) $(DEBUGFLAGS) ./src/gh_test_src.cpp -o ./exe/gh_test_src $(LFLAGS)

clean:
	rm -f ./exe/*
