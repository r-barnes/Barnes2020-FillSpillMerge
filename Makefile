RICHDEM_GIT_HASH="-NA-"
RICHDEM_COMPILE_TIME=`date -u +'%Y-%m-%d %H:%M:%S UTC'`
export RD_CXX_FLAGS=-I../common/richdem/include -DRICHDEM_GIT_HASH="\"$(RICHDEM_GIT_HASH)\"" -DRICHDEM_COMPILE_TIME="\"$(RICHDEM_COMPILE_TIME)\""
export CXXFLAGS=--std=c++17 -g -O3 -march=native -Wall -Wno-unknown-pragmas #-fopenmp #-fsanitize=address
export LIBS=-lnetcdf


a.out: main.cpp DisjointDenseIntSet.hpp dephier.hpp ../common/netcdf.hpp Makefile fill_spill_merge.hpp
	g++-7 $(CXXFLAGS) $(RD_CXX_FLAGS) main.cpp ../common/richdem/include/richdem/richdem.cpp $(LIBS)

clean:
	rm a.out
