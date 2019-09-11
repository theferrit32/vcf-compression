default: main

main: src/main.cpp
	g++ -Wall -o main src/main.cpp

huffman.so: huffman.cpp
	g++ -fPIC -shared -o huffman.so huffman.cpp
