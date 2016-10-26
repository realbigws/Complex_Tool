

# make the subdir first
CXX = g++
#FLAG = -O3 -static 
FLAG = -O3
SHELL :=/bin/bash

PDBToolDir = ./Complex_Tool_src
PDBToolEXE = /Complex_Tool

all: PDB_Tool

PDB_Tool:
	cd $(PDBToolDir)&&make
	mv ${PDBToolDir}$(PDBToolEXE) ../


# clean all the temp files
clean:
