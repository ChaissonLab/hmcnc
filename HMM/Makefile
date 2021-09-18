all: samToBed viterbi hmmcnc

hmmcnc: hmmcnc.cpp
	g++ -std==c++14 -O2 -I $(CONDA_PREFIX)/include $^ -o $@ -L $(CONDA_PREFIX)/lib -lhts -lpthread -ffast-math

viterbi: viterbi.cpp
	g++ -O2 -I $(CONDA_PREFIX)/include $^ -o $@

samToBed: SamToBed.cpp
	g++ $^ -o $@
