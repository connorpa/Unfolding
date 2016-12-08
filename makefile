CXX = g++ -g -Wall -std=c++1y
CFLAGS = $(shell root-config --cflags)
LIBS = $(shell root-config --libs)
BOOST = /usr/lib64/
ROOUNFOLDPATH = $(shell echo ${ROOUNFOLD})

unfolding: parser.o resolution.o unfolding.o spectrum.o efficiency.o sigma.o
	${CXX} -o $@ $^ ${LIBS} -L${ROOUNFOLDPATH} -lRooUnfold -L${BOOST} 

parser.o: parser.cpp
	${CXX} -c $^ ${CFLAGS}

resolution.o: resolution.cpp
	${CXX} -c $^ ${CFLAGS}

spectrum.o: spectrum.cpp
	${CXX} -c $^ ${CFLAGS}

sigma.o: sigma.cpp
	${CXX} -c $^ ${CFLAGS}

efficiency.o: efficiency.cpp
	${CXX} -c $^ ${CFLAGS}

unfolding.o: unfolding.cpp
	${CXX} -c $^ ${CFLAGS} -I${ROOUNFOLDPATH}/src

clean:
	rm -rf *.o *.exe *.h.gch *.pcm *.d *_cpp.so

.PHONY: clean
