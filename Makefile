all:
	g++ -g -O2 -fPIC -std=c++11 -c ./src/2drft.cpp -o ./obj/2drft.o
	cython -2 --cplus ./src/cpy_2drft.pyx --embed -o ./src/cpy_2drft.cpp
	g++ -g -O2 -fPIC -std=c++11 -c ./src/cpy_2drft.cpp -o ./obj/cpy_2drft.o `python-config --includes`
	g++ -g -O2 -std=c++11 -shared -o cpy_2drft.so ./obj/cpy_2drft.o ./obj/2drft.o `python-config --libs` -lfftw3

