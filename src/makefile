CXXFLAGS = -std=c++11 -g -O0 -Wpedantic -Wall -Wextra -Winit-self -Winline -Wconversion -Weffc++ -Wctor-dtor-privacy -Woverloaded-virtual -Wconversion -Wsign-promo
TARGET = hmm_test
LDFLAGS = -lm

${TARGET}: constants.o hmm_algorithms.o hmm_test.o
	${CXX} -o $@ $^ ${LDFLAGS}

.PHONY: clean mrproper

clean:	
	rm -f *.o

mrproper: clean
	rm -f ${TARGET}