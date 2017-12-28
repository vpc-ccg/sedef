CXX=icpc
CFLAGS=-c -I fmt -I . -std=c++14 -I src -fopenmp
LDFLAGS=-lrt -lz -fopenmp

GIT_VERSION:=$(shell git describe --dirty --always --tags)
SOURCES:=$(wildcard src/*.cc) $(wildcard extern/*.cc) 
OBJECTS=$(SOURCES:.cc=.o) patterns.o
EXECUTABLE=sedef

all: CFLAGS+=-g -O2
all: $(SOURCES) $(EXECUTABLE)

sanitize: CXX=g++
sanitize: CFLAGS+= -g -O1 -fno-omit-frame-pointer -fsanitize=address
sanitize: LDFLAGS+=-fno-omit-frame-pointer -fsanitize=address
sanitize: $(SOURCES) $(EXECUTABLE)

release: CFLAGS+=-g -O3 -DNDEBUG 
release: $(SOURCES) $(EXECUTABLE)

debug: CFLAGS+=-g
debug: $(SOURCES) $(EXECUTABLE)

superdebug: CFLAGS+=-g -O0 -fno-inline
superdebug: $(SOURCES) $(EXECUTABLE)

profile: CFLAGS+=-g -pg -O2
profile: LDFLAGS+=-pg
profile: $(SOURCES) $(EXECUTABLE)

gprofile: LDFLAGS+=-ltcmalloc -lprofiler
gprofile: CFLAGS+=-g -O2
gprofile: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@

.cc.o:	
	$(CXX) $(CFLAGS) -DGITVER=\"$(GIT_VERSION)\" $< -o $@

clean:
	find . -name '*.o' -delete
	rm -rf $(EXECUTABLE) $(TESTEXE) gmon.out* 

patterns.o:
	ld -r -b binary -o patterns.o patterns.bin