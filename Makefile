default: main

main: main.cpp
	g++ -Wall -o main main.cpp

huffman.so: huffman.cpp
	g++ -fPIC -shared -o huffman.so huffman.cpp
