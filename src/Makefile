CC = g++
CPPFLAGS = -O3 -DNDEBUG -fopenmp -std=c++0x -I $(PATH_TO_LEMON) -DMASSIVE_DATA_RLCSA
LDFLAGS = -fopenmp

PATH_TO_LEMON = ./lemon_binaries_linux/include/

OBJECTS_ASSEMBLER = complete-contigs.o utils.o OptionParser.o 
OBJECTS_MAXIMALITY = maximality.o rlcsa/rlcsa.a

EXECUTABLES = complete-contigs maximality

# use default make rules
# see https://www.gnu.org/software/make/manual/html_node/Catalogue-of-Rules.html#Catalogue-of-Rules

build: $(EXECUTABLES)

complete-contigs: $(OBJECTS_ASSEMBLER)
maximality: $(OBJECTS_MAXIMALITY)

rlcsa/rlcsa.a:
	make -C rlcsa build_rlcsa
	
clean:
	rm -f $(OBJECTS_ASSEMBLER) $(OBJECTS_MAXIMALITY) $(EXECUTABLES)

.PHONY: build clean
