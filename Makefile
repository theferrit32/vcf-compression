
default: main

all: main release

main: src/main.cpp src/utils.hpp
	g++ -Wall -std=c++11 -DDEBUG -o $@ $^

release: main_release

main_release: src/main.cpp src/utils.hpp
	g++ -Wall -std=c++11 -O3 -o $@ $^

huffman.so: huffman.cpp
	g++ -fPIC -shared -o huffman.so huffman.cpp

clean:
	rm -f main main_release
