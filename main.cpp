#include "fluid.h"

int main(int argc, char* argv[]){
    int maxParticles = 100;
    double timeStep = 100;
    double stepSize = 0.1;
    Scene *s = new Scene(maxParticles, timeStep, stepSize);
}