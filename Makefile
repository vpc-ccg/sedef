CXX=icpc
CPPFLAGS=-c -I fmt -I . -std=c++14 -I src -fopenmp
LDFLAGS=-lrt -lz -fopenmp

GIT_VERSION:=$(shell git describe --dirty --always --tags)
SOURCES:=$(wildcard src/*.cc) $(wildcard extern/*.cc) 
OBJECTS:=$(SOURCES:.cc=.o) patterns.o
DEP := $(OBJECTS:.o=.d)

CPPFLAGS += -MMD -MP -I.

.PHONY: all clean

EXECUTABLE=sedef
EXE=sedef

all: CPPFLAGS+=-g -O2
all: $(SOURCES) $(EXECUTABLE)

sanitize: CXX=g++
sanitize: CPPFLAGS+= -g -O1 -fno-omit-frame-pointer -fsanitize=address
sanitize: LDFLAGS+=-fno-omit-frame-pointer -fsanitize=address
sanitize: $(SOURCES) $(EXECUTABLE)

release: CPPFLAGS+=-g -O3 -DNDEBUG 
release: $(SOURCES) $(EXECUTABLE)

release2: CPPFLAGS+=-g -O3 -DNDEBUG -DATTAINER
release2: EXE=sedef2
release2: $(SOURCES) $(EXECUTABLE)

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

LIB = libsedef

lib: CPPFLAGS += -O2 -g -fPIC
lib: $(SOURCES) $(LIB)

PYTHON_VERSION = 2.7
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)

$(LIB): $(OBJECTS)
	$(CXX) -I$(PYTHON_INCLUDE) -I. -fPIC -c scratch/sedef.cpp -std=c++14 -o scratch/sedef.po
	$(CXX) -shared -Wl,--export-dynamic scratch/sedef.po $(OBJECTS) -lboost_python -lrt -lz -L/usr/lib/python$(PYTHON_VERSION)/config -lpython$(PYTHON_VERSION) -fopenmp -o $@.so

$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $(EXE)

.cc.o:	
	$(CXX) $(CPPFLAGS) -DGITVER=\"$(GIT_VERSION)\" $< -o $@

-include $(DEP)

clean:
	rm -rf $(EXECUTABLE) $(TESTEXE) $(OBJECTS) $(DEP)

patterns.o:
	ld -r -b binary -o patterns.o patterns.bin