all: a b c d e f seq


run:
	g++ -pthread lupthread.cpp
	#./a.out 1000 2
	
a:
	g++ -pthread lupthread.cpp
	./a.out 2000 2
b:
	g++ -pthread lupthread.cpp
	./a.out 2000 4
	
c:
	g++ -pthread lupthread.cpp
	./a.out 2000 6
d:
	g++ -pthread lupthread.cpp
	./a.out 2000 8

e:
	g++ -pthread lupthread.cpp
	./a.out 2000 16

f:
	g++ -pthread lupthread.cpp
	./a.out 2000 32
	
seq:
	g++ -pthread seqpthread.cpp
	./a.out 2000

