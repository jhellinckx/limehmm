CC = g++ -std=c++11
CPPFLAGS = -Wpedantic -Wall -Wextra -Winit-self -Winline -Wconversion -Weffc++ -Wctor-dtor-privacy -Woverloaded-virtual -Wconversion -Wsign-promo


all: hmm_test
	./hmm_test

hmm_test: hmm_test.cpp
	${CC} ${CPPFLAGS} $^ -o $@ 

clean: 
	rm -f hmm_test