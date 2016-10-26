

# make the subdir first
CXX = g++
#FLAG = -O3 -static 
FLAG = -O3
SHELL :=/bin/bash

ComplexToolDir = ./Complex_Tool_src
ComplexToolEXE = /Complex_Tool

all: Complex_Tool PDB_To_SEQ_miss

Complex_Tool:
	cd $(ComplexToolDir)&&make
	mv ${ComplexToolDir}$(ComplexToolEXE) ../

PDB_To_SEQ_miss:
	${CXX} PDB_To_SEQ_miss.cpp ${FLAG} -o PDB_To_SEQ_miss
	mv PDB_To_SEQ_miss ../ 

# clean all the temp files
clean:
