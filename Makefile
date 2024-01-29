# derived from hal2vg/Makefile
rootPath = ./

all : vg2maf

libbdsgPath = deps/libbdsg-easy
libbdsgLibs = ${libbdsgPath}/lib/libbdsg.a ${libbdsgPath}/lib/libhandlegraph.a ${libbdsgPath}/lib/libsdsl.a ${libbdsgPath}/lib/libdivsufsort.a ${libbdsgPath}/lib/libdivsufsort64.a

# Multithreading with OpenMP.
PARALLEL_FLAGS=-fopenmp -pthread
LIBS=-L$(CWD)

ifeq ($(shell uname -s),Darwin)
    # Our compiler might be clang that lacks -fopenmp support.
    # Sniff that
    ifeq ($(strip $(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l)), 1)
        # The compiler complained about fopenmp instead of its nonsense input file.
        # We need to use the hard way of getting OpenMP not bundled with the compiler.

        # The compiler only needs to do the preprocessing
        PARALLEL_FLAGS = -Xpreprocessor -fopenmp -pthread

        ifeq ($(shell if [ -d /opt/local/lib/libomp ];then echo 1;else echo 0;fi), 1)
            # Use /opt/local/lib/libomp if present, because Macports installs libomp there.
            # Brew is supposed to put it somewhere the compiler can find it by default.
            LIBS += -L/opt/local/lib/libomp
            # And we need to find the includes. Homebrew puts them in the normal place
            # but Macports hides them in "libomp"
            PARALLEL_FLAGS += -I/opt/local/include/libomp
        endif

        # We also need to link it
        LIBS += -lomp

    endif
endif

CXX ?= g++
ifeq (${CXXFLAGS},)
CXXFLAGS = 
endif
#CXXFLAGS += -O0 -fno-inline -fno-omit-frame-pointer -fsanitize=address
CXXFLAGS += -O3
CXXFLAGS += -Werror=return-type -std=c++14 -ggdb -g -MMD -MP $(PARALLEL_FLAGS)

CXXFLAGS += -I deps/taffy/taffy/submodules/sonLib/C/inc/ -I deps/taffy/taffy/inc -I deps/libbdsg-easy/include -I deps/libvgio/include/ -I deps/abPOA/include

# for abpoa
export avx2 = 1
CXXFLAGS += -mavx2

static:
	CFLAGS="$${CFLAGS} -static" \
	CXXFLAGS="$${CXXFLAGS} -static" \
	${MAKE} all

check-static: static
	if [ $(shell ls vg2maf | xargs ldd 2>& 1 | grep "not a dynamic" | wc -l) = $(shell ls hal2vg clip-vg halRemoveDupes halMergeChroms halUnclip filter-paf-deletions count-vg-hap-cov | wc -l) ] ; then\
		echo "ldd verified that all files in bin/ are static";\
	else\
		echo "ldd found dynamic linked binary in bin/";\
		exit 1;\
	fi

cleanFast : 
	rm -f vg2maf *.o

clean :
	rm -f vg2maf *.o deps/htslib/libhts.a deps/htslib/configure deps/taffy/lib/libstTaf.a ${libbdsgPath}/lib/libbdsg.a deps/libvgio/build/libvgio.a
	cd deps/libbdsg-easy && ${MAKE} clean
	cd deps/taffy && ${MAKE} clean
	cd deps/libvgio && ${MAKE} clean
	cd deps/htslib && rm configure && ${MAKE} clean
	cd deps/abPOA && ${MAKE} clean

vg2maf_main.o : vg2maf_main.cpp vg2maf.hpp stream_index.hpp scanner.hpp vg_types.hpp deps/libvgio/build/libvgio.a ${libbdsgPath}/lib/libbdsg.a deps/taffy/lib/libstTaf.a 
	${CXX} ${CXXFLAGS} -I . vg2maf_main.cpp -c

vg2maf.o : vg2maf.cpp vg2maf.hpp stream_index.hpp scanner.hpp vg_types.hpp deps/libvgio/build/libvgio.a ${libbdsgPath}/lib/libbdsg.a deps/taffy/lib/libstTaf.a 
	${CXX} ${CXXFLAGS} -I . vg2maf.cpp -c

insertions.o : insertions.cpp vg2maf.hpp vg_types.hpp deps/libvgio/build/libvgio.a ${libbdsgPath}/lib/libbdsg.a deps/taffy/lib/libstTaf.a 
	${CXX} ${CXXFLAGS} -I . insertions.cpp -c

scanner.o : scanner.cpp scanner.hpp vg_types.hpp scanner.hpp deps/libvgio/build/libvgio.a ${libbdsgPath}/lib/libbdsg.a deps/taffy/lib/libstTaf.a 
	${CXX} ${CXXFLAGS} -I . scanner.cpp -c

stream_index.o : stream_index.cpp stream_index.hpp scanner.hpp vg_types.hpp deps/libvgio/build/libvgio.a ${libbdsgPath}/lib/libbdsg.a deps/taffy/lib/libstTaf.a
	${CXX} ${CXXFLAGS} -I . stream_index.cpp -c

# cant build static without going through configure (todo: maybe can set a env var)
deps/htslib/configure :
	cd deps/htslib && autoreconf --install

deps/htslib/libhts.a : deps/htslib/configure 
	cd deps/htslib && ./configure --disable-libcurl && ${MAKE}

deps/taffy/lib/libstTaf.a:
	cd deps/taffy && ${MAKE}

${libbdsgPath}/lib/libbdsg.a:
	cd deps/libbdsg-easy && ${MAKE}

deps/libvgio/build/libvgio.a:
	cd deps/libvgio && rm -rf build && mkdir build && cd build && cmake .. && ${MAKE} && cp vg.pb.h ../include/vg

deps/abPOA/lib/libabpoa.a:
	cd deps/abPOA && ${MAKE}

vg2maf : vg2maf_main.o vg2maf.o insertions.o scanner.o stream_index.o deps/taffy/lib/libstTaf.a ${libbdsgPath}/lib/libbdsg.a deps/libvgio/build/libvgio.a deps/htslib/libhts.a deps/abPOA/lib/libabpoa.a
	${CXX} ${CXXFLAGS} vg2maf_main.o vg2maf.o insertions.o stream_index.o scanner.o deps/taffy/lib/libstTaf.a deps/taffy/lib/libsonLib.a ${libbdsgLibs} -ljansson deps/libvgio/build/libvgio.a deps/abPOA/lib/libabpoa.a -lprotobuf deps/htslib/libhts.a -lcurl -lm -lz -llzma -lbz2 -ldeflate -fopenmp -pthread -o vg2maf

all : vg2maf

test : vg2maf
	cd test && python3 vg2mafTest.py
