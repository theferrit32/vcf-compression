default: debug

debug: src/main.cpp
	g++ -Wall -o -DDEBUG main src/main.cpp

release: src/main.cpp
	g++ -Wall -O3 -o main src/main.cpp

huffman.so: huffman.cpp
	g++ -fPIC -shared -o huffman.so huffman.cpp
