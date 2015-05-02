CC = g++
FILES = main.cpp utils.cpp OptionParser.cpp 
OUT_EXE = complete-contigs

CPPFLAGS = -O3 -fopenmp -std=c++0x -I ~/local/include/
# -Wall

build: $(FILES)
	mkdir -p ../bin; cd ../src/; $(CC) $(CPPFLAGS) -o ../bin/$(OUT_EXE) $(FILES) -Wall

clean:
	rm -f ../bin/$(OUT_EXE)

rebuild: clean build