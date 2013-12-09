//ffmpg

#include "fluid.h"

inline float sqr(float x) { return x*x; }



Scene::Scene(int p, double t, double s){
    maxParts = p;
    timeStep = t;
    step = s;
    particles = new vector<Particle *>();
    film = new Film(WIDTH, HEIGHT);
    init();
    //render();
}

void Scene::init(){
    // Color colour;
    // colour.r = (double)1.0;
    // colour.g = (double)1.0;
    // colour.b = (double)1.0;
    // vector<vector<Color> > m = vector<vector<Color> > (HEIGHT, vector<Color>(WIDTH, colour));
    srand(time(NULL));
    //map<Vector3f, bool> positions;
    //for(int t = 0; t < timeStep; t++){
        for(int i= 0 ; i < maxParts; i ++){
            // int x = rand() % 50 + WIDTH/ 2 - 25;//(RIGHT - LEFT - RADIUS)) + LEFT + RADIUS;
            // int y = rand() % 50 + HEIGHT - 100; //(TOP-BOTTOM - RADIUS)) + BOTTOM + RADIUS ;
            double x = fRand(-0.1, 0.1);
            double y = fRand(0.5, 0.6);
            double z = fRand(-1.1, -1.6);
            // double x = fRand(300.0, 350.0);
            // double y = fRand(200.0, 275.0);
            // double z = fRand(20.0, 30.0);
            Vector3f pos(x, y, z);
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

    // for(int i = 0 ; i < particles->size(); i ++){
    //     Particle *p = particles->at(i);
    //     Vector3f pos = p->getPosition();

    // }
    // Particle *p = new Particle(MASS, Vector3f(0, 0, -0.1), Vector3f(0, 0, 0));
    // Vector3f pos = p->getPosition();
    // cout << "Here" << endl;
    // glPushMatrix();
    //     glTranslated(pos.x(), pos.y(), pos.z());
    //     glutSolidSphere(RADIUS,10, 10);
    // glPopMatrix();

    // double maxD = -109238;
    // double minD = 109238;
    //neighbourAndDist * n = new neighbourAndDist();
    //n->p = p;
    //n->dist = 0;
    //map<Vector3f, bool> positions;

    for(int t = 0; t < timeStep; t++){  //for every timestep
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);               // clear the color buffer
        glMatrixMode(GL_MODELVIEW);                   // indicate we are specifying camera transformations
        glLoadIdentity();
        vector<vector<Color> > m = vector<vector<Color> > (WIDTH, vector<Color>(HEIGHT, colour));
        //for(int i= 0 ; i < maxParts; i ++){
        vector<vector<Particle * > > neighbors;
        //cout << particles->size() << endl;
        //= vector<vector<Particle * > >(particles->size());
        for(int i = 0; i < particles->size(); i++){  //for every particle
            //cout << "i: " << i << endl;
            //double density = MASS;

            Particle *particle = particles->at(i);
            double density = MASS;
            vector<Particle *> findNeighs;
            for(int j = 0; j < particles->size(); j++){  //comparison to all other particles to see if they're close enough to effect the density
                Particle *tempParticle = particles->at(j);
                double dist = particle->getDistance(*tempParticle);

                if (dist <= H && i != j){  //if the particle is close enough, add its mass * kernel to the density
                    //neighbourAndDist * nAD = new neighbourAndDist();
                    //nAD->p = tempParticle;
                    //nAD->dist = dist;

                    double kern = particle->getKernel(dist);
                    // cout << "dist: " << dist << endl;
                    // cout << "kern: " << kern << endl;
                    //cout << "kern: " << kern << endl;
                    density += tempParticle->getMass() * kern;
                    findNeighs.push_back(tempParticle);
                }

            }
            //density += particle->getMass() * particle->getKernel(0);
            // cout << "density: " << density << endl;
            //cout << "1. " << i << ": " << density << endl;
            //particle->setDensity(density); //(*particle).setDensity(density);
            //cout << "2. " << i << ": " << particle->getDensity() << endl;
            neighbors.push_back(findNeighs);
            particle->setDensity(density);
        }

        //second iteration of particles and only their neighbors

        for(int i = 0; i < particles->size(); i++){
            //cout << "3. " << i << ": " << neighbors[i][0]->getDensity() << endl;
            Particle *particle = particles->at(i);
            //cout << "1. " << i << ": " << particle->getDensity() << endl;
            Vector3f position = particle->getPosition();
            //cout << t << ". " << position << endl << endl;
            Vector3f velocity = particle->getVelocity();
            // maxD = max(maxD, particle->getDensity());
            // minD = min(minD, particle->getDensity());
            //cout << i << ": " << position << endl;
            glPushMatrix();
                glTranslated(position.x(), position.y(), position.z());
                glutSolidSphere(RADIUS, 10, 10);
            glPopMatrix();


            // for(int j = -RADIUS; j < RADIUS; j++){
            //     for(int k = -RADIUS; k < RADIUS; k++){
            //         if (particle->getDistance(position + Vector3f(j, 0, 0) + Vector3f(0, k, 0)) < RADIUS) {
            //             if(abs(velocity.x()) > 20 || abs(velocity.y()) > 20 || abs(velocity.z()) > 20 ){
            //                 m[position.x() + j][position.y()+ k].r = 1.0;
            //                 m[position.x() + j][position.y() + k].g = 0;
            //                 m[position.x() + j][position.y() + k].b = 0;
            //             }else {
            //                 m[position.x() + j][position.y()+ k].r = 0;
            //                 m[position.x() + j][position.y() + k].g = 0;
            //                 m[position.x() + j][position.y() + k].b = 1.0;
            //             }
            //             // if(i == 0){
            //             //     m[position.x() + j][position.y()+ k].r = 1.0;
            //             //     m[position.x() + j][position.y() + k].g = 0;
            //             //     m[position.x() + j][position.y() + k].b = 0;
            //             // } else if(i == 1){
            //             //     m[position.x() + j][position.y()+ k].r = 0;
            //             //     m[position.x() + j][position.y() + k].g = 1.0;
            //             //     m[position.x() + j][position.y() + k].b = 0;
            //             // } else if(i == 2){
            //             //     m[position.x() + j][position.y()+ k].r = 0;
            //             //     m[position.x() + j][position.y() + k].g = 0;
            //             //     m[position.x() + j][position.y() + k].b = 1.0;
            //             // }
            //         }
            //     }
            // }


            Vector3f viscosityForce = Vector3f::Zero();
            //double pressure = 0;
            Vector3f pressureForce = Vector3f::Zero();
            double pressureJ = particle->calcPressure();
            //vector<neighbourAndDist * > currNeighs = neighbours[i];

            vector<Particle * > curNeighs = neighbors[i];
            //cout << i << " " << curNeighs.size() << endl;

            for(int j = 0; j < particles->size(); j++){//curNeighs.size(); j++){//currNeighs.size(); j++){
                Particle *tempParticle = particles->at(j);//curNeighs[j];// currNeighs[j]->p;
                //cout << "3. " << i << " " << j << ": " << tempParticle->getDensity() << endl;
                double tempMass = tempParticle->getMass();
                double tempDens = tempParticle->getDensity();
                Vector3f tempVel = tempParticle->getVelocity();
                double dist = particle->getDistance(*tempParticle);

                if(dist <= H && i!=j){

                    //Pressure
                    Vector3f rij = tempParticle->getPosition() - position;
                    Vector3f kernDerive = particle->getKernDerive(dist, rij);
                    //cout << kernDerive << endl;
                    double pressureK = tempParticle->calcPressure();
                    //cout << tempDens << endl;
                    //pressure += ((pressureJ + pressureK) / 2 )* tempMass / tempDens * kernDerive;
                    //pressure += tempMass * (pressureJ + pressureK) / (2 * tempDens) * kernDerive;
                    pressureForce += tempMass * (pressureJ + pressureK) / (2 * tempDens) * kernDerive;

                    //Viscosity
                    double kernSecond = particle->getKernSecond(dist);

                    viscosityForce += (tempVel - velocity) * tempMass / tempDens * kernSecond;
                }

            }

            pressureForce *= -1;

            viscosityForce *= VISC;

            Vector3f gravityForce(0, particle->getDensity() * GRAVITY, 0);
            //Vector3f pressureForce(0, pressure, 0);//pressure);
            cout << i << ": ";
            cout << "pForce: " << pressureForce << endl;
            //cout << "dens: " << particle->getDensity() << endl;
             cout << "glForce: " << gravityForce << endl;
             cout << "vForce: " << viscosityForce << endl;

            Vector3f totalForce = gravityForce + pressureForce + viscosityForce;
            //cout << "totalForce: " << totalForce << endl;
            Vector3f acceleration = totalForce/particle->getDensity();
            //cout << "1. " << particle->getVelocity() << endl;
            velocity = velocity + DELTAT * acceleration;  //maybe implement some kind of terminal velocity?
            //cout << "2. " << velocity << endl;

            position = position + DELTAT * velocity;
            //cout << "Original: " << particle->getPosition() << endl;
            //cout << "New: " << position << endl;

            //cout << particle->getVelocity().z() << endl;
         //cout << "Here" << endl;
            if((position.x() - RADIUS <= LEFT) || (position.y() - RADIUS <= BOTTOM) || (position.x() + RADIUS >= RIGHT) || (position.y() + RADIUS >= TOP)){
                    velocity *= -0.4;
                    while((position.x() - RADIUS <= LEFT) || (position.y() - RADIUS <= BOTTOM) || (position.x() + RADIUS >= RIGHT) || (position.y() + RADIUS >= TOP)){
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
            // for(int i = LEFT; i < RIGHT; i ++){
            //     m[i][BOTTOM].r = 0;
            //     m[i][BOTTOM].g = 0;
            //     m[i][BOTTOM].b = 0;
            //     m[i][TOP].r = 0;
            //     m[i][TOP].g = 0;
            //     m[i][TOP].b = 0;
            // }
            // for(int i = BOTTOM; i < TOP; i++){
            //     m[LEFT][i].r = 0;
            //     m[LEFT][i].g = 0;
            //     m[LEFT][i].b = 0;
            //     m[RIGHT][i].r = 0;
            //     m[RIGHT][i].g = 0;
            //     m[RIGHT][i].b = 0;
            // }
        FreeImage_Initialise();
        BYTE* pixels = new BYTE[ 3 * WIDTH * HEIGHT];

        glReadPixels(0, 0, WIDTH, HEIGHT, GL_BGR, GL_UNSIGNED_BYTE, pixels);
        FIBITMAP* image = FreeImage_ConvertFromRawBits(pixels, WIDTH, HEIGHT, 3 * WIDTH, 24, 0x0000FF, 0xFF0000, 0x00FF00, false);
        // get a string to concat with an int
        stringstream nameStream;
        nameStream << "test" << t << ".png";
        string name = nameStream.str();
        char *a = new char[name.size() + 1];
        a[name.size()] = 0;
        memcpy(a, name.c_str(), name.size());
        if(FreeImage_Save(FIF_BMP, image, a, 0))
            cout << "Image " << t << " successfully saved! " << endl ;
        FreeImage_DeInitialise(); //Cleanup !

            //film->saveImage(m);
        glFlush();
        glutSwapBuffers();
        glPopMatrix();
    }
    // cout << maxD << endl;
    // cout << minD << endl;
}


double Scene::fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


Particle::Particle(double m, Vector3f p, Vector3f v){
    mass = m;
    position = p;
    velocity = v;
    density = 0;
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
    //return STIFFNESS * ( pow(density/ IDEALDENSITY, 7) - 1 );
    return STIFFNESS * (density - IDEALDENSITY);
}

double Particle::getKernel(double r){
    //return exp(-4 * pow(r / H, 2));  //Gaussian Kernel form (from intel)
    //double  q = r / H;
    //return pow((1 - pow(q,2)) , 3);
    double c = 315 / (64 * M_PI * pow((double)H, 9.0));
    return c * pow((H * H - r * r), 3);
}

Vector3f Particle::getKernDerive(double r, Vector3f rij){
    //double q = r / H;
    //return -8 * q  * exp(-4 * pow(q, 2));
    //return -6 * q * pow((1 - pow(q, 2)), 2);
    double c = -45 / (M_PI * pow((double) H, 6.0));
    return c * rij / r * (pow(H - r, 2));
}

double Particle::getKernSecond(double r){
//8 e^(-4 q^2) (-1+8 q^2)
    //double q = r / H;
    //return 8 * exp(-4 * pow(q, 2) * (-1 + 8 * pow(q, 2)));
   // return -6 * (5 * pow(q, 4) - 6 * pow(q, 2) + 1);
    double c = 45 / (M_PI * pow((double)H, 6.0));
    return c * (H - r);
}
