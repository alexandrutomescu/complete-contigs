CC = g++
FILES_ASSEMBLER = main.cpp utils.cpp OptionParser.cpp 
FILES_MAXIMALITY = maximality.cpp
OUT_EXE_ASSEMBLER = complete-contigs
OUT_EXE_MAXIMALITY = maximality
PATH_TO_LEMON = ./lemon_binaries_linux/include/

CPPFLAGS_ASSEMBLER = -O3 -fopenmp -std=c++0x -I $(PATH_TO_LEMON)

CPPFLAGS_MAXIMALITY = -O3 -fopenmp -std=c++0x -DMASSIVE_DATA_RLCSA

COMPILE_RLCSA = mkdir -p ../bin/rlcsa; cp -r rlcsa/ ../bin/; cd ../bin/rlcsa; make build_rlcsa; cd ../../src/

.PHONY: build maximality clean

build:
	mkdir -p ../bin; cd ../src/; $(CC) $(CPPFLAGS_ASSEMBLER) -o ../bin/$(OUT_EXE_ASSEMBLER) $(FILES_ASSEMBLER)

maximality:
	$(COMPILE_RLCSA); $(CC) $(CPPFLAGS_MAXIMALITY) -o ../bin/$(OUT_EXE_MAXIMALITY) $(FILES_MAXIMALITY) ../bin/rlcsa/rlcsa.a

clean:
	rm -f ../bin/$(OUT_EXE_ASSEMBLER); rm -f ../bin/$(OUT_EXE_MAXIMALITY); rm -fr ../bin/rlcsa