SOURCE = src/main.cpp src/utils.cpp \
	src/compress.cpp src/sparse.cpp \
	src/string_t.c src/split_iterator.cpp

default: main

all: main release

main: $(SOURCE)
	g++ -Wall -std=c++11 -D_GNU_SOURCE -DDEBUG -o $@ $^

release: main_release

main_release: $(SOURCE)
	g++ -Wall -std=c++11 -D_GNU_SOURCE -O3 -o $@ $^

clean:
	rm -f main main_release
