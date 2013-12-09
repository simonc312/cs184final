#include <iostream>
#include "FreeImage/FreeImage.h"
#include "Eigen/Eigen"
#include <vector>
#include <string>
#include <cmath>  //may not be needed
#include <climits>
#include <stack>
#include <map>
#include <cstdlib>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <cstring>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include <GLUT/glut.h>
#include <OpenGL/glu.h>

#include <time.h>
#include <math.h>


#define _USE_MATH_DEFINES

// #define MAXINITPART 1
// #define MASS 10
// #define BPP 24
// #define H 20
// #define GRAVITY -10
// #define IDEALDENSITY 275
// #define STIFFNESS 100
// #define VISC 20
// #define DELTAT 0.2
// #define RADIUS 4

#define MAXINITPART 1
#define MASS 0.002
#define BPP 24
#define H 0.25//0.0625
#define GRAVITY -10
#define IDEALDENSITY 3.5
#define STIFFNESS 0.2//0.4
#define VISC 1
#define DELTAT 0.01
#define RADIUS 0.01


#define LEFT -0.5
#define TOP 0.7
#define RIGHT 0.5
#define BOTTOM -0.6
#define WIDTH 800
#define HEIGHT 800
#define EPISILON 3



using namespace Eigen;
using namespace std;

class Particle;
class Film;
class Scene;
class Viewport;


typedef struct{
    double r;
    double g;
    double b;
    double a;
} Color;

typedef struct neighbourAndDist{
    Particle *p;
    double dist;
} neighbourAndDist;

class Scene{
public:
    Scene(int p, double t, double s);
    vector<Particle *> *particles;
    void init();
    void render();


private:
    Film *film;
    int maxParts;
    double timeStep;
    double step;
    double fRand(double fMin, double fMax);
};


class Particle{
public:
    Particle(double m, Vector3f p, Vector3f v);
    double getDensity();
    double getDistance(Particle p);
    double getDistance(Vector3f v);
    double getMass();
    Vector3f getPosition();
    Vector3f getVelocity();
    void setMass(double m);
    void setPosition(Vector3f p);
    void setVelocity(Vector3f v);
    void setDensity(double d);
    double getKernel(double r);
    Vector3f getKernDerive(double r,  Vector3f rij);
    double getKernSecond(double r);
    double calcPressure();
private:
    double mass, density;
    Vector3f position;
    Vector3f velocity;
};


class Film{
public:
    Film(int w, int h);
    void saveImage(vector<vector<Color> > m);
private:
    int width, height, count;
    string name;
};

class Viewport {
  public:
    int w, h; // width and height
};
