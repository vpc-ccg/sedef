CXX=icpc
CPPFLAGS=-c -I fmt -I . -std=c++14 -I src -fopenmp
LDFLAGS=-lrt -lz -fopenmp

GIT_VERSION:=$(shell git describe --dirty --always --tags)
SOURCES:=$(wildcard src/*.cc) $(wildcard extern/*.cc) 
OBJECTS=$(SOURCES:.cc=.o) patterns.o
DEP := $(OBJECTS:.o=.d)

CPPFLAGS += -MMD -MP -I.

.PHONY: all clean

EXECUTABLE=sedef

all: CPPFLAGS+=-g -O2
all: $(SOURCES) $(EXECUTABLE)

sanitize: CXX=g++
sanitize: CPPFLAGS+= -g -O1 -fno-omit-frame-pointer -fsanitize=address
sanitize: LDFLAGS+=-fno-omit-frame-pointer -fsanitize=address
sanitize: $(SOURCES) $(EXECUTABLE)

release: CPPFLAGS+=-g -O3 -DNDEBUG 
release: $(SOURCES) $(EXECUTABLE)

debug: CPPFLAGS+=-g
debug: $(SOURCES) $(EXECUTABLE)

superdebug: CPPFLAGS+=-g -O0 -fno-inline
superdebug: $(SOURCES) $(EXECUTABLE)

profile: CPPFLAGS+=-g -pg -O2
profile: LDFLAGS+=-pg
profile: $(SOURCES) $(EXECUTABLE)

gprofile: CXX=g++
gprofile: LDFLAGS=-Wl,--no-as-needed,-lprofiler,--as-needed -ltcmalloc -lrt -lz -fopenmp
gprofile: CPPFLAGS+=-g -O1
gprofile: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@

.cc.o:	
	$(CXX) $(CPPFLAGS) -DGITVER=\"$(GIT_VERSION)\" $< -o $@

-include $(DEP)

clean:
	rm -rf $(EXECUTABLE) $(TESTEXE) $(OBJECTS) $(DEP)

patterns.o:
	ld -r -b binary -o patterns.o patterns.bin