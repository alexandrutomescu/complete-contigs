CC = g++
FILES = main.cpp
OUT_EXE = maximality

CPPFLAGS = -O3 -fopenmp -std=c++0x -DMASSIVE_DATA_RLCSA -Wall

COMPILE_RLCSA = rm -rf ../bin/rlcsa; mkdir -p ../bin/rlcsa; cp -r rlcsa/ ../bin/; cd ../bin/rlcsa; make build_rlcsa

build: $(FILES)
	mkdir -p ../../bin/; $(COMPILE_RLCSA); cd ../../src/; cd ../src/maximality/; $(CC) $(CPPFLAGS) -o ../bin/$(OUT_EXE) $(FILES) ../bin/rlcsa/rlcsa.a

quick: $(FILES)
	mkdir -p ../../bin/; cd ../src/maximality/; $(CC) $(CPPFLAGS) -o ../bin/$(OUT_EXE) $(FILES) ../bin/rlcsa/rlcsa.a

clean:
	rm -f $(OUT_EXE); rm -rf ../bin/rlcsa;

rebuild: clean build