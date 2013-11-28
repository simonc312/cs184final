CC = g++
RM = /bin/rm -f
CFFLAGS = -I ./Eigen -O2 -fopenmp
LDFLAGS = -freeimage ./FreeImage/libfreeimage.a

all: main

main: main.o scene.o film.o
	$(CC) $(CFFLAGS) main.o scene.o film.o -o main $(LDFLAGS)

main.o: main.cpp
	$(CC) -c main.cpp

scene.o: scene.cpp
	$(CC) -c scene.cpp

film.o: film.cpp
	$(CC) -c film.cpp

clean:
	$(RM) *.o main