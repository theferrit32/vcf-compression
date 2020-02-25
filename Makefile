
default: main

all: main release

main: src/main.cpp src/utils.hpp src/string_t.c src/string_t.h
	g++ -Wall -std=c++11 -D_GNU_SOURCE -DDEBUG -o $@ $^

release: main_release

main_release: src/main.cpp src/utils.hpp src/string_t.c src/string_t.h
	g++ -Wall -std=c++11 -D_GNU_SOURCE -O3 -o $@ $^

huffman.so: huffman.cpp
	g++ -fPIC -shared -o huffman.so huffman.cpp

clean:
	rm -f main main_release
