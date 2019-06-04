export GDAL_LIBS=`gdal-config --libs` `gdal-config --dep-libs` -DUSEGDAL
export GDAL_CFLAGS=`gdal-config --cflags`
RICHDEM_GIT_HASH="-NA-"
RICHDEM_COMPILE_TIME=`date -u +'%Y-%m-%d %H:%M:%S UTC'`
export RD_CXX_FLAGS=-I../common/richdem/include -DRICHDEM_GIT_HASH="\"$(RICHDEM_GIT_HASH)\"" -DRICHDEM_COMPILE_TIME="\"$(RICHDEM_COMPILE_TIME)\""
export CXXFLAGS=--std=c++17 -g -O3 -march=native -Wall -Wno-unknown-pragmas -pipe -fopenmp #-fsanitize=address
export DS_LIB_FLAGS=-Iparallel-hashmap
#export LIBS=-lnetcdf


a.out: main.cpp DisjointDenseIntSet.hpp dephier.hpp ../common/netcdf.hpp Makefile fill_spill_merge.hpp
	$(CXX) $(GDAL_CFLAGS) $(DS_LIB_FLAGS) $(CXXFLAGS) $(RD_CXX_FLAGS) -o dephier.exe main.cpp ../common/richdem/include/richdem/richdem.cpp $(LIBS) $(GDAL_LIBS)

clean:
	rm a.out
