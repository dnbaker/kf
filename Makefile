.PHONY=all tests clean obj update
CXX=g++
CC=gcc

GMATCH=$(findstring g++,$(CXX))

CLHASH_CHECKOUT = "&& git checkout master"
WARNINGS=-Wall -Wextra -Wno-char-subscripts \
		 -Wpointer-arith -Wwrite-strings -Wdisabled-optimization \
		 -Wformat -Wcast-align -Wno-unused-function -Wno-unused-parameter \
		 -pedantic -DUSE_PDQSORT -Wunused-variable \
		-Wduplicated-branches -Wdangling-else  # -Wconversion
ifndef EXTRA
	EXTRA:=
endif
ifndef INCPLUS
	INCPLUS:=
endif
ifndef EXTRA_LD
	EXTRA_LD:=
endif
DBG:=
OS:=$(shell uname)
FLAGS=

OPT_MINUS_OPENMP= -O3 -funroll-loops\
	  -pipe -fno-strict-aliasing -march=native -mpclmul $(FLAGS) $(EXTRA)
OPT=$(OPT_MINUS_OPENMP) -fopenmp
XXFLAGS=-fno-rtti
CXXFLAGS=$(OPT) $(XXFLAGS) -std=c++11 $(WARNINGS)
CCFLAGS=$(OPT) $(CFLAGS) -std=c11 $(WARNINGS)
LIB=-lz
LD=-L. $(EXTRA_LD)

OBJS=$(patsubst %.c,%.o,$(wildcard src/*.c)) $(patsubst %.cpp,%.o,$(wildcard src/*.cpp))
DOBJS=$(patsubst %.c,%.do,$(wildcard src/*.c)) $(patsubst %.cpp,%.do,$(wildcard src/*.cpp))
ZOBJS=$(patsubst %.c,%.zo,$(wildcard src/*.c)) $(patsubst %.cpp,%.zo,$(wildcard src/*.cpp))

TEST_OBJS=$(patsubst %.cpp,%.o,$(wildcard test/*.cpp))
ZTEST_OBJS=$(patsubst %.cpp,%.zo,$(wildcard test/*.cpp))

EXEC_OBJS=$(patsubst %.cpp,%.o,$(wildcard bin/*.cpp))

EX=$(patsubst bin/%.cpp,%,$(wildcard bin/*.cpp))
DEX=$(patsubst %_d,%,$(EX))
ZEX=$(patsubst %_z,%,$(EX))
STATEX=$(patsubst %,%_s,$(EX))
STATEXZ=$(patsubst %,%_sz,$(EX))

INCLUDE=-I. -Iinclude



all: $(OBJS) $(EX) $(DEX)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -DNDEBUG -c $< -o $@ $(LIB)

%: bin/%.cpp
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(OBJS) -DNDEBUG $< -o $@ $(LIB)

%_s: bin/%.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(OBJS) -static-libstdc++ -static-libgcc -DNDEBUG $< -o $@ $(LIB)

%_d: bin/%.cpp $(DOBJS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(DOBJS) $< -o $@ $(LIB)

INCLUDES=-I`python3-config --includes` -Ipybind11/include -I.
SUF=`python3-config --extension-suffix`
OBJS=$(patsubst %.cpp,%$(SUF),$(wildcard *.cpp))
ifeq ($(shell uname),Darwin)
    UNDEFSTR=-undefined dynamic_lookup
else
    UNDEFSTR=
endif

python: py/kf.cpython.so
	python -c "import subprocess;import site; subprocess.check_call('cp "py/*`python3-config --extension-suffix`" %s' % site.getsitepackages()[0], shell=True)"

%.cpython.so: %.cpp
	$(CXX) $(UNDEFSTR) $(INCLUDES) -O3 -Wall $(FLAGS) $(INC) -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` $< -o $*$(SUF) && \
    ln -fs $*$(SUF) $@

tests: clean unit

unit: $(DOBJS) $(TEST_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TEST_OBJS) $(LD) $(DOBJS) -o $@ $(LIB)

zunit: $(ZOBJS) $(ZTEST_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ZTEST_OBJS) $(LD) $(ALL_ZOBJS) -o $@ $(LIB) $(ZCOMPILE_FLAGS)

clean:
	rm -f $(ZOBJS) $(ZTEST_OBJS) $(ZW_OBJS) $(OBJS) $(DEX) $(ZEX) $(EX) $(TEST_OBJS) $(DOBJS) unit lib/*o src/*o \
	&& cd ../hll && make clean && cd ../bonsai

mostlyclean:
	rm -f $(ZOBJS) $(ZTEST_OBJS) $(ZW_OBJS) $(OBJS) $(DEX) $(ZEX) $(EX) $(TEST_OBJS) $(DOBJS) unit lib/*o src/*o
