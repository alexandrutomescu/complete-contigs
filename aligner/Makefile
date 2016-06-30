CC = g++
FILES = main.cpp OptionParser.cpp
OUT_EXE = aligner

CPPFLAGS = -Wall -O3 -std=c++0x -DMASSIVE_DATA_RLCSA

build:
	$(CC) $(CPPFLAGS) -o $(OUT_EXE) $(FILES) rlcsa/rlcsa.a

clean:
	rm -f aligner