all: samToBed viterbi hmmcnc

hmmcnc.o: hmmcnc.cpp
	g++ -c -g2 -I $(CONDA_PREFIX)/include $^ -o $@

hmcnc_defs.o: hmcnc_def.cpp
	g++ -c -g2 -I $(CONDA_PREFIX)/include $^ -o $@

hmcnc_io.o: hmcnc_io.cpp
	g++ -c -g2 -I $(CONDA_PREFIX)/include $^ -o $@

hmmcnc: main.cpp hmmcnc.o hmcnc_defs.o hmcnc_io.o
	g++ -g2 -I $(CONDA_PREFIX)/include $^ -o $@ -L $(CONDA_PREFIX)/lib -lhts -lpthread -ffast-math
#	g++ -fsanitize=address -g -I $(CONDA_PREFIX)/include $^ -o $@ -L $(CONDA_PREFIX)/lib -lhts -lpthread -ffast-math

viterbi: viterbi.cpp
	g++ -O2 -I $(CONDA_PREFIX)/include $^ -o $@

samToBed: SamToBed.cpp
	g++ $^ -o $@
