#include "fluid.h"
#define LEFT  50
#define TOP 450
#define RIGHT 600
#define BOTTOM 50
#define WIDTH 480
#define HEIGHT 640
#define EPISILON 3

Scene::Scene(int p, double t, double s){
    maxParts = p;
    timeStep = t;
    step = s;
    particles = new vector<Particle *>();
    film = new Film(HEIGHT, WIDTH);
    init();
    render();
}

void Scene::init(){
    // Color colour;
    // colour.r = (double)1.0;
    // colour.g = (double)1.0;
    // colour.b = (double)1.0;
    // vector<vector<Color> > m = vector<vector<Color> > (HEIGHT, vector<Color>(WIDTH, colour));
    srand(time(NULL));
    //map<Vector3f, bool> positions;
    //for(int i = 0; i < timeStep; i++){
        for(int j= 0 ; j < maxParts; j ++){
            int x = rand() % 50 + HEIGHT / 2 - 25;
            int y = rand() % 50 + WIDTH - 100;
            Vector3f pos(x, y, 0);
            //std::pair<std::map<char, int>::iterator, bool> ret;
            //ret = positions.insert(std::pair<Vector3f, bool>(pos, true));
            //if(ret.second == true){
                Particle *p = new Particle(MASS, pos, Vector3f(0, 0, 0));
                particles->push_back(p);
                // m[x][y].r = 0;
                // m[x][y].g = 0;
                // m[x][y].b = 1.0;
            //}
        }
        //cout << particles->size() << endl;
        //film->saveImage(m);
    //}
}

// This is the hash function taken from the paper by kelager its from email
// page 47
// However the unorder map data structure I use from c++ has its own default
// hash function so I don't know if this will work...somehow make our own custom
// hash function and embed it to the data structure?
//int hashFunction(Vector3f pos){
  //  return (((int)pos.x()*73856093) xor ((int)pos.y()*19349663) xor ((int)pos.z()*83492791)) % getNextPrime(particles//->size());
//}

//int getNextPrime(int n){
// this function should return the closest prime number >= n;
    //it's hard finding code to use on the internet for this function which is
    //ridiculous
  //  return 1;
//}



void Scene::render(){
    Color colour;
    colour.r = (double)1.0;
    colour.g = (double)1.0;
    colour.b = (double)1.0;

    // The spatial hash map to add particles to
    //typedef std::unordered_map<int,Particle> neighbormap;


    Particle *p = new Particle(MASS, Vector3f(0, 0, 0), Vector3f(0, 0, 0));
    //neighbourAndDist * n = new neighbourAndDist();
    //n->p = p;
    //n->dist = 0;
    //map<Vector3f, bool> positions;
    for(int i = 0; i < timeStep; i++){  //for every timestep
        vector<vector<Color> > m = vector<vector<Color> > (HEIGHT, vector<Color>(WIDTH, colour));
        //for(int j= 0 ; j < maxParts; j ++){
        vector<vector<Particle * > > neighbors;
        //= vector<vector<Particle * > >(particles->size());
        for(int j = 0; j < particles->size(); j++){  //for every particle
            //cout << "j: " << j << endl;
            double density = 0;
            Particle *particle = particles->at(j);
            vector<Particle *> findNeighs;
            for(int k = 0; k < particles->size(); k++){  //comparison to all other particles to see if they're close enough to effect the density
                Particle *tempParticle = particles->at(k);
                double dist = particle->getDistance(*tempParticle);

                if (dist < H){  //if the particle is close enough, add its mass * kernel to the density
                    //neighbourAndDist * nAD = new neighbourAndDist();
                    //nAD->p = tempParticle;
                    //nAD->dist = dist;

                    double kern = particle->getKernel(dist);
                    //cout << "kern: " << kern << endl;
                    density += tempParticle->getMass() * kern;
                    findNeighs.push_back(tempParticle);
                }

            }
            //density += particle->getMass() * particle->getKernel(0);
            //cout << "density: " << density << endl;
            particle->setDensity(density);
            neighbors.push_back(findNeighs);
        }


        //second iteration of particles and only their neighbors

        for(int j = 0; j < particles->size(); j++){
            Particle *particle = particles->at(j);

            Vector3f position = particle->getPosition();
            Vector3f velocity = particle->getVelocity();


            for(int k = -RADIUS; k < RADIUS; k++){
               for(int l = -RADIUS; l < RADIUS; l++){
                    if (particle->getDistance(position + Vector3f(k, 0, 0) + Vector3f(0, l, 0)) < RADIUS) {
                        if(abs(velocity.x()) > 20 || abs(velocity.y()) > 20 || abs(velocity.z()) > 20 ){
                            m[position.x() + k][position.y()+ l].r = 1.0;
                            m[position.x() + k][position.y() + l].g = 0;
                            m[position.x() + k][position.y() + l].b = 0;
                            // m[position.x()][position.y()].r = 1.0;
                            // m[position.x()][position.y()].g = 0;
                            // m[position.x()][position.y()].b = 0;
                        }else {
                            m[position.x() + k][position.y()+ l].r = 0;
                            m[position.x() + k][position.y() + l].g = 0;
                            m[position.x() + k][position.y() + l].b = 1.0;
                            // m[position.x()][position.y()].r = 0;
                            // m[position.x()][position.y()].g = 0;
                            // m[position.x()][position.y()].b = 1.0;
                        }
                   }
               }
            }


            Vector3f viscosityForce = Vector3f::Zero();
            double pressure = 0;
            double pressureJ = particle->calcPressure();
            //vector<neighbourAndDist * > currNeighs = neighbours[j];

            vector<Particle * > curNeighs = neighbors[j];


            for(int k = 0; k < curNeighs.size(); k++){//currNeighs.size(); k++){
                Particle *tempParticle = curNeighs[k];// currNeighs[k]->p;
                double tempMass = tempParticle->getMass();
                double tempDens = tempParticle->getDensity();
                Vector3f tempVel = tempParticle->getVelocity();
                double dist = particle->getDistance(*tempParticle);

                //if(dist < H){


                    //Pressure
                    double kernDerive = particle->getKernDerive(dist);
                    //cout << kernDerive << endl;
                    double pressureK = tempParticle->calcPressure();
                    //cout << tempDens << endl;
                    pressure += ((pressureJ + pressureK) / 2 )* tempMass / tempDens * kernDerive;

                    //Viscosity
                    double kernSecond = particle->getKernSecond(dist);

                    viscosityForce += (tempVel - velocity) * tempMass / tempDens * kernSecond;
                //}

            }

            pressure *= -1;

            viscosityForce *= VISC;

            Vector3f gravityForce(0, particle->getDensity() * GRAVITY, 0);
            Vector3f pressureForce(pressure, pressure, 0);//pressure);
            //cout << "pForce: " << pressureForce << endl;
            //cout << "glForce: " << gravityForce << endl;
            //cout << "vForce: " << viscosityForce << endl;

            Vector3f totalForce = /*gravityForce + */pressureForce;/* + viscosityForce;*/
            //cout << "totalForce: " << totalForce << endl;
            Vector3f acceleration = totalForce/particle->getDensity();
            //cout << "1. " << particle->getVelocity() << endl;
            velocity = velocity + DELTAT * acceleration;  //maybe implement some kind of terminal velocity?
            //cout << "2. " << velocity << endl;




            position = position + DELTAT * velocity;
            //cout << "Original: " << particle->getPosition() << endl;
            //cout << "New: " << position << endl;

            //cout << particle->getVelocity().z() << endl;

            if((position.x() <= LEFT) || (position.y() <= BOTTOM) || (position.x() >= RIGHT) || (position.y() >= TOP) ){
                    velocity *= -0.1;
                    while((position.x() <= LEFT) || (position.y() <= BOTTOM) || (position.x() >= RIGHT) || (position.y() >= TOP) ){
                        cout << position << endl;
                        position = position + DELTAT * velocity;
                    }
                    //cout << "NEW: " << position << endl << endl;
                //cout << "2. " << velocity << endl;
            }

            //if((position.x() > 0) && (position.y() > 0) && (position.x() < HEIGHT) && (position.y() < WIDTH)){

            particle->setPosition(position);
            particle->setVelocity(velocity);
           // }
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
            for(int j = LEFT; j < RIGHT; j ++){
                m[j][BOTTOM].r = 0;
                m[j][BOTTOM].g = 0;
                m[j][BOTTOM].b = 0;
                m[j][TOP].r = 0;
                m[j][TOP].g = 0;
                m[j][TOP].b = 0;
            }
            for(int j = BOTTOM; j < TOP; j++){
                m[LEFT][j].r = 0;
                m[LEFT][j].g = 0;
                m[LEFT][j].b = 0;
                m[RIGHT][j].r = 0;
                m[RIGHT][j].g = 0;
                m[RIGHT][j].b = 0;
            }

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

double Particle::calcPressure(){
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
