#include "fluid.h"

int main(int argc, char* argv[]){
    int maxParticles = atoi(argv[1]);
    double timeStep = atof(argv[2]);
    double stepSize = 0.1;
    Scene *s = new Scene(maxParticles, timeStep, stepSize);
}