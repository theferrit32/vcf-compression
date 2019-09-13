default: main

main: src/main.cpp
	g++ -Wall -DDEBUG -o $@ $^

main_release: src/main.cpp
	g++ -Wall -O3 -o $@ $^

huffman.so: huffman.cpp
	g++ -fPIC -shared -o huffman.so huffman.cpp

clean:
	rm -f main main_release