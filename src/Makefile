CC = clang++
CPPFLAGS = -std=c++11 -Wpedantic -Wall -Wextra -Winit-self -Winline -Wconversion -Weffc++ -Wctor-dtor-privacy -Woverloaded-virtual -Wconversion -Wsign-promo


all: main

main: main.cpp
	${CC} ${CPPFLAGS} $^ -o $@ 

GRAPHS := $(patsubst %.dot,%.eps,$(wildcard *.dot))

%.eps: %.dot
	dot -Teps $< > $@

graph:	${GRAPHS}

clean: 
	rm -f main && rm -f *.eps