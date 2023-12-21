# derived from hal2vg/Makefile
rootPath = ./

all : vg2maf

libbdsgPath = deps/libbdsg-easy
libbdsgLibs = ${libbdsgPath}/lib/libbdsg.a ${libbdsgPath}/lib/libhandlegraph.a ${libbdsgPath}/lib/libsdsl.a ${libbdsgPath}/lib/libdivsufsort.a ${libbdsgPath}/lib/libdivsufsort64.a

CXX ?= g++
CXXFLAGS += -std=c++14 -I deps/taffy/taffy/inc -I deps/libbdsg-easy/include -I deps/libvgio/include/ 

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
	rm -f vg2maf

clean :
	rm -f vg2maf *.o
	cd deps/libbdsg-easy && make clean
	cd deps/taffy && make clean
	cd deps/libvgio && make clean

vg2maf.o : vg2maf.cpp 
	${CXX} ${CXXFLAGS} -I . vg2maf.cpp -c

deps/taffy/lib/libstTaf.a:
	cd deps/taffy && ${MAKE}

${libbdsgPath}/lib/libbdsg.a:
	cd deps/libbdsg-easy && ${MAKE}

deps/libvgio/build/libvgio.a:
	cd deps/libvgio && rm -rf build && mkdir build && cd build && cmake .. && ${MAKE}

vg2maf : vg2maf.o deps/taffy/lib/libstTaf.a ${libbdsgPath}/lib/libbdsg.a deps/libvgio/build/libvgio.a
	${CXX} ${CXXFLAGS} -lm -lz -llzma -lbz2 -ldeflate -fopenmp -pthread vg2maf.o deps/taffy/lib/libstTaf.a ${libbdsgLibs} -ljansson deps/libvgio/build/libvgio.a -lprotobuf -lhts -o vg2maf

#test :
#	make
#	cd tests && prove -v t
