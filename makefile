# CC=g++
# CFLAGS = -std=c++14 -O3 -Wno-unused-result -fdiagnostics-color=always
# LFLAGS = -lstdc++

# all: main

CC=g++
CFLAGS = -std=c++11 -O3 -Wno-unused-result -fdiagnostics-color=always
DEBUGFLAGS = -std=c++11 -g -Wno-unused-result -fdiagnostics-color=always

all: ./src/gh_test_src.cpp
	$(CC) $(CFLAGS) ./src/gh_test_src.cpp -o ./exe/gh_test_src

debug: ./src/gh_test_src.cpp
	$(CC) $(DEBUGFLAGS) ./src/gh_test_src.cpp -o ./exe/gh_test_src

clean:
	rm -f ./exe/*
