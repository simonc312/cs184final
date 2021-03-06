//ffmpg
#include "fluid.h"

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
    srand(time(NULL));
        for(int i= 0 ; i < maxParts; i ++){
            // double x = fRand(-0.1, 0.1);
            // double y = fRand(0.5, 0.6);
            // double z = fRand(-1.1, -1.6);
            // double x = fRand(WIDTH/2 - 25,WIDTH/2 + 25);
            // double y = fRand(HEIGHT-175, HEIGHT-125);
            // double z = fRand(-30, -35);
            double x = fRand(WIDTH/2 - 25,WIDTH/2 + 25);
            double y = fRand(400, 450);
            double z = fRand(-80, -125);
            Vector3f pos(x, y, z);
            Particle *p = new Particle(MASS, pos, Vector3f(0, 0, 0));
            particles->push_back(p);
        }
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
    for(int t = 0; t < timeStep; t++){  //for every timestep
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);               // clear the color buffer
        glMatrixMode(GL_MODELVIEW);                   // indicate we are specifying camera transformations
        glLoadIdentity();
        gluLookAt(0.0, 1.0, -0.25, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0 );
        // if(t % 5 == 0) {
        //     init();
        // }

        //Draw the boundaries
        drawBoundaries();

        vector<vector<Particle * > > neighbors;

        //density calculations
        #pragma omp for
        for(int i = 0; i < particles->size(); i++){  //for every particle

            Particle *particle = particles->at(i);
            double density = MASS;
            vector<Particle *> findNeighs;

            for(int j = 0; j < particles->size(); j++){  //comparison to all other particles to see if they're close enough to effect the density
                Particle *tempParticle = particles->at(j);
                double dist = particle->getDistance(*tempParticle);

                if (dist <= H && i != j){  //if the particle is close enough, add its mass * kernel to the density
                    double kern = particle->getKernel(dist);
                    density += tempParticle->getMass() * kern;
                    findNeighs.push_back(tempParticle);
                }
            }
            neighbors.push_back(findNeighs);
            particle->setDensity(density);
        }

        //second iteration of particles and only their neighbors
        #pragma omp for
        for(int i = 0; i < particles->size(); i++){
            Particle *particle = particles->at(i);
            Vector3f position = particle->getPosition();
            Vector3f velocity = particle->getVelocity();

            //http://stackoverflow.com/questions/17565664/gluproject-and-2d-display
            //Render particle
            GLdouble posX, posY, posZ;//3D point
            posX=convert(position.x(), WIDTH);
            posY=convert(position.y(), HEIGHT);
            posZ=convert(position.z(), LENGTH);

            glPushMatrix();
                glTranslated(posX, posY, posZ);
                glutSolidSphere(SRADIUS, 10, 10);
            glPopMatrix();

            //Force calculations
            Vector3f viscosityForce = Vector3f::Zero();
            Vector3f pressureForce = Vector3f::Zero();
            Vector3f surfaceNormal = Vector3f::Zero();
            double colorField = 0;
            double pressureJ = particle->calcPressure();

            vector<Particle * > curNeighs = neighbors[i];

            for(int j = 0; j < curNeighs.size(); j++){//currNeighs.size(); j++){

                Particle *tempParticle = curNeighs[j];// currNeighs[j]->p;
                double tempMass = tempParticle->getMass();
                double tempDens = tempParticle->getDensity();
                Vector3f tempVel = tempParticle->getVelocity();
                double dist = particle->getDistance(*tempParticle);
                double kern = particle->getKernel(dist);
                colorField += tempMass / tempDens * kern;

                //Pressure and surfaceNormal
                Vector3f rij = tempParticle->getPosition() - position;
                Vector3f kernDerive = particle->getKernDerive(dist, rij);
                double pressureK = tempParticle->calcPressure();
                pressureForce += tempMass * (pressureJ + pressureK) / (2 * tempDens) * kernDerive;

                surfaceNormal += tempMass / tempDens * kernDerive;

                //Viscosity
                double kernSecond = particle->getKernSecond(dist);
                viscosityForce += (tempVel - velocity) * tempMass / tempDens * kernSecond;
            }

            pressureForce *= -1;
            viscosityForce *= VISC;

            Vector3f surfaceTension = Vector3f::Zero();
            Vector3f gravityForce(0, particle->getDensity() * GRAVITY, 0);
            // cout << "pForce: " << pressureForce << endl;
            // //cout << "dens: " << particle->getDensity() << endl;
            //  cout << "glForce: " << gravityForce << endl;
            //  cout << "vForce: " << viscosityForce << endl;

            //Update next position
            Vector3f totalForce = gravityForce + pressureForce + viscosityForce;
            //cout << "totalForce: " << totalForce << endl;
            Vector3f acceleration = totalForce/particle->getDensity();
            //cout << "1. " << particle->getVelocity() << endl;
            velocity = velocity + DELTAT * acceleration;  //maybe implement some kind of terminal velocity?
            //cout << "2. " << velocity << endl;
            Vector3f newPosition = position + DELTAT * velocity;

            //Boundary check next position
            bool bounce = false;
             while((newPosition.x() - RADIUS <= LEFT) || (newPosition.y() - RADIUS <= BOTTOM) || (newPosition.x() + RADIUS >= RIGHT) || (newPosition.y() + RADIUS >= TOP) || (newPosition.z() + RADIUS <= BACK) || (newPosition.z() - RADIUS >= FRONT)){
                bounce = true;
                velocity *= 0.9;
                newPosition = position + DELTAT * velocity;
             }
             if(bounce) velocity *= -1;

            // if(newPosition.y() - RADIUS <= BOTTOM){
            //     int boundTime = (BOTTOM - position.y()) / velocity;
            //     Vector3f collision = position + boundTime * velocity;
            //     Vector3f collNorm(collision.x(), 1, collision.z()).normalized();
            //     double penDist = newPosition.dist(collision);
            //     newPosition = newPosition + penDist * collNorm;
            //     velocity = velocity - 0.3 * (velocity.dot(collNorm)) * collNorm;
            // }

            particle->setPosition(newPosition);
            particle->setVelocity(velocity);
        }
        saveImage(t);

        glFlush();
        glutSwapBuffers();
        glPopMatrix();
    }
}

void Scene::drawBoundaries(){
        glDisable(GL_LIGHTING);

        glBegin(GL_QUADS); //bottom
            glColor3f(1.0, 0.0, 0.0);
            glVertex3d(convert(LEFT, WIDTH), convert(BOTTOM, HEIGHT), convert(FRONT, LENGTH));
            glVertex3d(convert(RIGHT, WIDTH), convert(BOTTOM, HEIGHT), convert(FRONT, LENGTH));
            glVertex3d(convert(RIGHT, WIDTH), convert(BOTTOM, HEIGHT), convert(BACK, LENGTH));
            glVertex3d(convert(LEFT, WIDTH), convert(BOTTOM, HEIGHT), convert(BACK, LENGTH));
        glEnd();
        glEnable(GL_LIGHTING);
        // // glBegin(GL_QUADS); //top
        // //     glVertex3d(convert(LEFT, WIDTH), convert(TOP, HEIGHT), convert(FRONT, LENGTH));
        // //     glVertex3d(convert(RIGHT, WIDTH), convert(TOP HEIGHT), convert(FRONT, LENGTH));
        // //     glVertex3d(convert(RIGHT, WIDTH), convert(TOP, HEIGHT), convert(BACK, LENGTH));
        // //     glVertex3d(convert(LEFT, WIDTH), convert(TOP, HEIGHT), convert(BACK, LENGTH));
        // // glEnd();

        // glBegin(GL_QUADS); //left
        //         // glColor3f(0.0, 1.0, 0.0);
        //     glVertex3d(convert(LEFT, WIDTH), convert(BOTTOM, HEIGHT), convert(FRONT, LENGTH));
        //     glVertex3d(convert(LEFT, WIDTH), convert(BOTTOM, HEIGHT), convert(BACK, LENGTH));
        //     glVertex3d(convert(LEFT, WIDTH), convert(TOP, HEIGHT), convert(BACK, LENGTH));
        //     glVertex3d(convert(LEFT, WIDTH), convert(TOP, HEIGHT), convert(FRONT, LENGTH));
        // glEnd();
        // glBegin(GL_QUADS); //right
        //     glVertex3d(convert(RIGHT, WIDTH), convert(BOTTOM, HEIGHT), convert(FRONT, LENGTH));
        //     glVertex3d(convert(RIGHT, WIDTH), convert(TOP, HEIGHT), convert(FRONT, LENGTH));
        //     glVertex3d(convert(RIGHT, WIDTH), convert(TOP, HEIGHT), convert(BACK, LENGTH));
        //     glVertex3d(convert(RIGHT, WIDTH), convert(BOTTOM, HEIGHT), convert(BACK, LENGTH));
        // glEnd();
        // glBegin(GL_QUADS); //back
        //     glVertex3d(convert(LEFT, WIDTH), convert(BOTTOM, HEIGHT), convert(BACK, LENGTH));
        //     glVertex3d(convert(RIGHT, WIDTH), convert(BOTTOM, HEIGHT), convert(BACK, LENGTH));
        //     glVertex3d(convert(RIGHT, WIDTH), convert(TOP, HEIGHT), convert(BACK, LENGTH));
        //     glVertex3d(convert(LEFT, WIDTH), convert(TOP, HEIGHT), convert(BACK, LENGTH));
        // glEnd();
        // glBegin(GL_QUADS); //front
        //     glVertex3d(convert(LEFT, WIDTH), convert(BOTTOM, HEIGHT), convert(FRONT, LENGTH));
        //     glVertex3d(convert(RIGHT, WIDTH), convert(BOTTOM, HEIGHT), convert(FRONT, LENGTH));
        //     glVertex3d(convert(RIGHT, WIDTH), convert(TOP, HEIGHT), convert(FRONT, LENGTH));
        //     glVertex3d(convert(LEFT, WIDTH), convert(TOP, HEIGHT), convert(FRONT, LENGTH));
        // glEnd();
}


void Scene::saveImage(int t){
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
}



double Scene::fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

double Scene::convert(double point, double comp){
    if(point >= comp/2){
        point = (point - comp/2)/(comp/2);
    } else{
        point = (point - comp/2)/(comp/2);
    }
    return point;
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
    return STIFFNESS * (density - IDEALDENSITY);
}

double Particle::getKernel(double r){
    double c = 315 / (64 * M_PI * pow((double)H, 9.0));
    return c * pow((H * H - r * r), 3);
}

Vector3f Particle::getKernDerive(double r, Vector3f rij){
    double c = -45 / (M_PI * pow((double) H, 6.0));
    return c * rij / r * (pow(H - r, 2));
}

double Particle::getKernSecond(double r){
    double c = 45 / (M_PI * pow((double)H, 6.0));
    return c * (H - r);
}
