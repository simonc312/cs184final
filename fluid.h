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

#define MAXINITPART 1
#define MASS 10
#define BPP 24
#define H 5
#define GRAVITY -5
#define RESTDENSITY 500
#define GASCONSTANT 500
#define VISC 0.7978
#define DELTAT 0.3
#define RADIUS 4

/*#define MAXINITPART 1
#define MASS 0.02
#define BPP 24
#define H 5 //0.0625
#define GRAVITY -1
#define RESTDENSITY 1000
#define GASCONSTANT 3.5
#define VISC 3.5
#define DELTAT 0.6
#define RADIUS 2*/



using namespace Eigen;
using namespace std;

class Particle;
class Film;
class Scene;

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
    double getKernDerive(double r);
    double getKernSecond(double r);
    double getPressure();
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
