CC = g++
FILES = main.cpp rlcsa/rlcsa.a
OUT_EXE = maximality

CPPFLAGS = -O3 -fopenmp -std=c++0x -I ~/local/include/ -DMASSIVE_DATA_RLCSA -Wall

build: $(FILES)
	$(CC) $(CPPFLAGS) -o ../bin/$(OUT_EXE) $(FILES)

clean:
	rm -f $(OUT_EXE)

rebuild: clean build