CC = g++
ifeq ($(shell sw_vers 2>/dev/null | grep Mac | awk '{ print $$2}'),Mac)
	CFLAGS = -freeimage ./FreeImage/libfreeimage.a -I ./Eigen -g -DGL_GLEXT_PROTOTYPES -I./include/ -I/usr/X11/include -DOSX -fopenmp
	LDFLAGS = -framework GLUT -framework OpenGL \
    	-L"/System/Library/Frameworks/OpenGL.framework/Libraries" \
    	-lGL -lGLU -lm -lstdc++
else
	CFLAGS = -g -DGL_GLEXT_PROTOTYPES -Iglut-3.7.6-bin
	LDFLAGS = -lglut -lGLU  -freeimage ./FreeImage/libfreeimage.a
endif
RM = /bin/rm -f

all: main

main: main.o scene.o film.o marchingCubes.o
	$(CC) $(CFLAGS) main.o scene.o film.o marchingCubes.o -o main $(LDFLAGS)

main.o: main.cpp
	$(CC) -c main.cpp

scene.o: scene.cpp
	$(CC) -c scene.cpp

film.o: film.cpp
	$(CC) -c film.cpp

marchingCubes.o: marchingCubes.cpp
	$(CC) -c marchingCubes.cpp

clean:
	$(RM) *.o main

images:
	$(RM) *.png