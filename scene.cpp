#include "fluid.h"

Scene::Scene(int p, double t, double s){
    maxParts = p;
    timeStep = t;
    step = s;
    particles = new vector<Particle *>();
    film = new Film(640, 480);
    init();
    render();
}

void Scene::init(){
    Color colour;
    colour.r = (double)1.0;
    colour.g = (double)1.0;
    colour.b = (double)1.0;
    vector<vector<Color> > m = vector<vector<Color> > (640, vector<Color>(480, colour));
    //map<Vector3f, bool> positions;
    //for(int i = 0; i < timeStep; i++){
        for(int j= 0 ; j < maxParts; j ++){
            int x = rand() % 450 + 50;
            int y = rand() % 200 +50;
            Vector3f pos(x, y, 0);
            //std::pair<std::map<char, int>::iterator, bool> ret;
            //ret = positions.insert(std::pair<Vector3f, bool>(pos, true));
            //if(ret.second == true){
                Particle *p = new Particle(MASS, pos, Vector3f(-10.0, 0, -0.1));
                particles->push_back(p);
                m[x][y].r = 0;
                m[x][y].g = 0;
                m[x][y].b = 1.0;
            //}
        }
        //cout << particles->size() << endl;
        //film->saveImage(m);
    //}
}

void Scene::render(){
    Color colour;
    colour.r = (double)1.0;
    colour.g = (double)1.0;
    colour.b = (double)1.0;

    //Particle *p = new Particle(MASS, Vector3f(0, 0, 0), Vector3f(0, 0, 0));
    //neighbourAndDist * n = new neighbourAndDist();
    //n->p = p;
    //n->dist = 0;
    //map<Vector3f, bool> positions;
    for(int i = 0; i < timeStep; i++){  //for every timestep
        vector<vector<Color> > m = vector<vector<Color> > (640, vector<Color>(480, colour));
        //for(int j= 0 ; j < maxParts; j ++){
        //vector<vector<neighbourAndDist * > > neighbours = vector<vector<neighbourAndDist * > >(100, vector<neighbourAndDist *>(100, n));
        for(int j = 0; j < particles->size(); j++){  //for every particle
            //cout << "j: " << j << endl;
            double density = 0;
            Particle *particle = particles->at(j);
            for(int k = 0; k < particles->size(); k++){  //comparison to all other particles to see if they're close enough to effect the density
                Particle *tempParticle = particles->at(k);
                double dist = particle->getDistance(*tempParticle);

                if (dist < H){  //if the particle is close enough, add its mass * kernel to the density
                    //neighbourAndDist * nAD = new neighbourAndDist();
                    //nAD->p = tempParticle;
                    //nAD->dist = dist;
                    //neighbours[j].push_back(nAD);
                    double kern = particle->getKernel(dist);
                    //cout << "kern: " << kern << endl;
                    density += tempParticle->getMass() * kern;
                }

            }
            //density += particle->getMass() * particle->getKernel(0);
            //cout << "density: " << density << endl;
            particle->setDensity(density);
        }


        for(int j = 0; j < particles->size(); j++){
            Particle *particle = particles->at(j);
            Vector3f gravityForce(0, particle->getDensity() * GRAVITY, 0);
            double pressure = 0;
            Vector3f viscosityForce;
            double pressureJ = particle->getPressure();
            //vector<neighbourAndDist * > currNeighs = neighbours[j];

            for(int k = 0; k < particles->size(); k++){//currNeighs.size(); k++){
                Particle *tempParticle = particles->at(k);// currNeighs[k]->p;
                double tempMass = tempParticle->getMass();
                double tempDens = tempParticle->getDensity();
                double dist = particle->getDistance(*tempParticle);
                if(dist < H){

                    //Pressure
                    double kernDerive = particle->getKernDerive(dist);
                    //cout << kernDerive << endl;
                    double pressureK = tempParticle->getPressure();
                    //cout << tempDens << endl;
                    pressure += (pressureJ + pressureK) / 2 * tempMass / tempDens * kernDerive;
                    //Viscosity
                    double kernSecond = particle->getKernSecond(dist);
                    viscosityForce += (tempParticle->getVelocity() - particle->getVelocity()) * tempMass / tempDens * kernSecond;
                }
            }

            pressure *= -1;
            viscosityForce *= VISC;


            Vector3f pressureForce(pressure, pressure, pressure);
            //cout << "pForce: " << pressureForce << endl;
            //cout << "glForce: " << gravityForce << endl;
            //cout << "vForce: " << viscosityForce << endl;

            Vector3f totalForce = gravityForce + pressureForce + viscosityForce;
            //cout << "totalForce: " << totalForce << endl;
            Vector3f acceleration = totalForce/particle->getDensity();
            //cout << "1. " << particle->getVelocity() << endl;
            Vector3f velocity = particle->getVelocity() + DELTAT * acceleration;
            //cout << "2. " << velocity << endl;

            if((particle->getPosition().x() < 50) || (particle->getPosition().y() < 50) || (particle->getPosition().x() > 600) || (particle->getPosition().y() > 450) ){
                    velocity *= -0.7;
                //cout << "2. " << velocity << endl;
            }


            Vector3f position = particle->getPosition() + DELTAT * velocity;
            //cout << "Original: " << particle->getPosition() << endl;
            //cout << "New: " << position << endl;
            particle->setPosition(position);



            if((position.x() > 0) && (position.y() > 0) && (position.x() < 640) && (position.y() < 480)){
                for(int k = -RADIUS; k < RADIUS; k++){
                    for(int l = -RADIUS; l < RADIUS; l++){
                        if (particle->getDistance(particle->getPosition() + Vector3f(k, 0, 0) + Vector3f(0, l, 0)) < RADIUS) {
                            m[position.x() + k][position.y()+ l].r = 0;
                            m[position.x() + k][position.y() + l].g = 0;
                            m[position.x() + k][position.y() + l].b = 1.0;
                        }
                    }
                }

                particle->setVelocity(velocity);
            }
        }

            //std::pair<std::map<char, int>::iterator, bool> ret;
            //ret = positions.insert(std::pair<Vector3f, bool>(pos, true));
            //if(ret.second == true){
                //Particle *p = new Particle(MASS, pos, Vector3f(0, 0, 0));
                //particles->push_back(p);
                //m[x][y].r = 0;
               // m[x][y].g = 0;
                //m[x][y].b = 1.0;
            //}
        //}
        //if(!particles->empty())
            film->saveImage(m);
    }
}




Particle::Particle(double m, Vector3f p, Vector3f v){
    mass = m;
    position = p;
    velocity = v;
}

double Particle::getDensity(){
    return density;
}

double Particle::getMass(){
    return mass;
}

Vector3f Particle::getPosition(){
    return position;
}

Vector3f Particle::getVelocity(){
    return velocity;
}

void Particle::setMass(double m){
    mass = m;
}

void Particle::setPosition(Vector3f p){
    position = p;
}

void Particle::setVelocity(Vector3f v){
    velocity = v;
}

void Particle::setDensity(double d){
    density = d;
}


double Particle::getDistance(Particle p){
    double x = pow(position.x() - p.getPosition().x(), 2);
    double y = pow(position.y() - p.getPosition().y(), 2);
    double z = pow(position.z() - p.getPosition().z(), 2);
    return sqrt(x + y + z);
}

double Particle::getDistance(Vector3f v){
    double x = pow(position.x() - v.x(), 2);
    double y = pow(position.y() - v.y(), 2);
    double z = pow(position.z() - v.z(), 2);
    return sqrt(x + y + z);
}

double Particle::getPressure(){
    return GASCONSTANT * ( pow(density/ RESTDENSITY, 7) - 1 );
}

double Particle::getKernel(double r){
    //return exp(-4 * pow(r / H, 2));  //Gaussian Kernel form (from intel)
    double  q = r / H;
    return pow((1 - pow(q,2)) , 3);
}

double Particle::getKernDerive(double r){
    double q = r / H;
    //return -8 * q  * exp(-4 * pow(q, 2));
    return -6 * q * pow((1 - pow(q, 2)), 2);
}

double Particle::getKernSecond(double r){
//8 e^(-4 q^2) (-1+8 q^2)
    double q = r / H;
    //return 8 * exp(-4 * pow(q, 2) * (-1 + 8 * pow(q, 2)));
    return -6 * (5 * pow(q, 4) - 6 * pow(q, 2) + 1);
}
