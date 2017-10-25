CC=g++
CFLAGS=-c -std=c++14  -I edlib/edlib/include
LDFLAGS=

GIT_VERSION:=$(shell git describe --dirty --always --tags)
SOURCES:=$(wildcard src/*.cc) fmt/fmt/format.cc edlib/edlib/src/edlib.cc
OBJECTS=$(SOURCES:.cc=.o)
EXECUTABLE=sedef-jaccard

all: CFLAGS+=-g -O2
all: $(SOURCES) $(EXECUTABLE)

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
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cc.o:	
	$(CC) $(CFLAGS) -DGITVER=\"$(GIT_VERSION)\" $< -o $@

clean:
	find . -name '*.o' -delete
	rm -rf $(EXECUTABLE) $(TESTEXE) gmon.out* 
