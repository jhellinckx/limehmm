ifeq (${CC},)
	CC = clang++
endif

ifeq (${CPPFLAGS},)
	CPPFLAGS = -std=c++11 -Wpedantic -Wall -Wextra -Winit-self -Winline -Wconversion -Weffc++ -Wctor-dtor-privacy -Woverloaded-virtual -Wconversion -Wsign-promo
endif

all: main

main: main.cpp
	${CC} ${CPPFLAGS} $^ -o $@ 

clean: 
	rm -f main