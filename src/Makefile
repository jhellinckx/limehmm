CC = clang++
CPPFLAGS = -std=c++11 -Wpedantic -Wall -Wextra -Winit-self -Winline -Wconversion -Weffc++ -Wctor-dtor-privacy -Woverloaded-virtual -Wconversion -Wsign-promo


all: main

main: main.cpp
	${CC} ${CPPFLAGS} $^ -o $@ 

clean: 
	rm -f main