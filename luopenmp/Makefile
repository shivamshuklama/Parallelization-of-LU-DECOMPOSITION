all: a b c d e f seq


run:
	g++ -fopenmp luopenmp.cpp
#	./a.out $A $B
	
a:
	g++ -fopenmp luopenmp.cpp
	./a.out 2000 2
b:
	g++ -fopenmp luopenmp.cpp
	./a.out 2000 4
	
c:
	g++ -fopenmp luopenmp.cpp
	./a.out 2000 6
d:
	g++ -fopenmp luopenmp.cpp
	./a.out 2000 8

e:
	g++ -fopenmp luopenmp.cpp
	./a.out 2000 16
	
f:
	g++ -fopenmp luopenmp.cpp
	./a.out 2000 32
	
seq:
	g++ -fopenmp seqopenmp.cpp
	./a.out 2000 1
