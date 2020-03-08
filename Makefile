CPP_FLAGS = -Wall -std=c++11 -D_GNU_SOURCE

SOURCE = src/main.cpp src/utils.cpp \
	src/compress.cpp src/sparse.cpp \
	src/string_t.c src/split_iterator.cpp


default: main

all: main release timing

main: $(SOURCE)
	g++ $(CPP_FLAGS) -DDEBUG -o $@ $^

release: main_release

main_release: $(SOURCE)
	g++ $(CPP_FLAGS) -O3 -o $@ $^

timing: main_timing

main_timing: $(SOURCE)
	g++ $(CPP_FLAGS) -O3 -DTIMING -o $@ $^

clean:
	rm -f main main_release main_timing
