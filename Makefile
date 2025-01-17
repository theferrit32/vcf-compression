CPP_FLAGS = -Wall -std=c++11 -D_GNU_SOURCE
flag = -D_FILE_OFFSET_BITS=64

SOURCE = src/main.cpp src/utils.cpp \
	src/compress.cpp src/sparse.cpp \
	src/string_t.c src/split_iterator.cpp


default: debug

all: debug release timing

debug: main_debug

main_debug: $(SOURCE)
	g++ $(CPP_FLAGS) -DDEBUG -DTIMING -o $@ $^

release: main_release

main_release: $(SOURCE)
	g++ $(CPP_FLAGS) -O3 -o $@ $^

timing: main_timing

main_timing: $(SOURCE)
	g++ $(CPP_FLAGS) -O3 -DTIMING -o $@ $^

uniqc: src/uniqc.cpp
	g++ $(CPP_FLAGS) -O3 -o $@ $^

clean:
	rm -f main main_debug main_release main_timing

test: test.cpp
	g++ $(CPP_FLAGS) -DDEBUG -DTIMING -o test test.cpp