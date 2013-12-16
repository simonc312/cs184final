//ffmpg
#include "fluid.h"

Scene::Scene(int p, double t, double s, int m){
    maxParts = p;
    timeStep = t;
    step = s;
    march = m;
    particles = new vector<Particle *>();
    film = new Film(WIDTH, HEIGHT);
    cubes = new Cubes();
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
            // double x = fRand(385, 415);
            // double y = fRand(425, 465);
            // double z = fRand(-385, -415);
            double x = fRand(350, 450);
            double y = fRand(450, 475);
            double z = fRand(-350, -450);
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
        gluLookAt(0.0, 0.15, -1.20, 0.0, 0.15, -10.0, 0.0, 1.0, 0.0);
        //gluLookAt(0.0, 0.6, -1.25, .0, 0.0, -2.1, 0.0, 1.0, 0.0);
        // gluLookAt()
        // if(t % 5 == 0) {
        //     init();
        // }

        //Draw the boundaries
       // if (t < 100)
        drawBoundaries();

        vector<vector<Particle * > > neighbors;
        double totalDens = 0;

        //density calculations
        #pragma omp parallel for
        for(int i = 0; i < particles->size(); i++){  //for every particle

            Particle *particle = particles->at(i);
            double density = MASS;
            vector<Particle *> findNeighs;

            for(int j = 0; j < particles->size(); j++){  //comparison to all other particles to see if they're close enough to effect the density
                Particle *tempParticle = particles->at(j);
                double dist = particle->getDistance(*tempParticle);

                if (dist <= H && i != j){  //if the particle is close enough, add its mass * kernel to the density
                    double kern = particle->getKernel(dist);
                    //cout << dist << endl;
                    density += tempParticle->getMass() * kern;
                    findNeighs.push_back(tempParticle);
                }
            }
            neighbors.push_back(findNeighs);
            particle->setDensity(density);
            totalDens += density;
        }
        cout << totalDens / particles->size() << endl;

        if(march == 1){

            int x = (RIGHT - LEFT) / GRID + 1;
            int y = (TOP - BOTTOM)/GRID + 1;
            int z = (abs(BACK - FRONT))/GRID + 1;
            cout << x << " " << y << " " << z << endl;

            Particle * grid[x][y][z];

            for(int i = 0; i < x; i++){
                for(int j = 0; j < y; j++){
                    for(int k = 0; k < z; k++){
                        Particle* particle = new Particle(MASS, Vector3f(i  * GRID + LEFT, j * GRID + BOTTOM, k * GRID + BACK), Vector3f(0, 0, 0));
                        double density = MASS;
                         for(int l = 0; l < particles->size(); l++){  //naive
                                Particle *tempParticle = particles->at(l);
                                double dist = particle->getDistance(*tempParticle);

                                if (dist <= H){  //if the particle is close enough, add its mass * kernel to the density
                                    double kern = particle->getKernel(dist);

                                    density += tempParticle->getMass() * kern;
                                }
                           }
                           particle->setDensity(density);
                           grid[i][j][k] = particle;
                           //cout << density << endl;
                    }
                }
            }


            for(int i = LEFT; i < RIGHT; i += GRID){
                for(int j = BOTTOM; j < TOP; j += GRID){
                    for(int k = BACK; k < FRONT; k +=GRID){
                        GRIDCELL g;
                        int x1 = (i - LEFT)/GRID;
                        int x2 = (i + GRID - LEFT)/GRID;
                        int y1 = (j - BOTTOM)/GRID;
                        int y2 = (j + GRID - BOTTOM)/GRID;
                        int z1 = (k - BACK)/GRID;
                        int z2 = (k + GRID - BACK)/GRID;

                        Particle *p = grid[x1][y1][z1];
                        if(p->getDistance(Vector3f(i, j, k)) != 0 ){
                            cerr << "Values are not the same" << endl;
                            cout << "Particle: " << p->getPosition() << endl;
                            cout << "Vector: " << Vector3f(i, j, k) << endl;
                        }
                        g.p[0] = Vector3f(p->getPosition());
                        g.val[0] = p->getDensity();

                        p = grid[x2][y1][z1];
                        if(p->getDistance(Vector3f(i + GRID, j, k)) != 0 ){
                            cerr << "Values are not the same" << endl;
                            cout << "Particle: " << p->getPosition() << endl;
                            cout << "Vector: " << Vector3f(i, j, k) << endl;
                        }
                        g.p[1] = Vector3f(p->getPosition());
                        g.val[1] = p->getDensity();

                        p = grid[x2][y1][z2];
                        if(p->getDistance(Vector3f(i + GRID, j, k + GRID)) != 0 ){
                            cerr << "Values are not the same" << endl;
                            cout << "Particle: " << p->getPosition() << endl;
                            cout << "Vector: " << Vector3f(i, j, k) << endl;
                        }
                        g.p[2] = Vector3f(p->getPosition());
                        g.val[2] = p->getDensity();

                        p = grid[x1][y1][z2];
                        if(p->getDistance(Vector3f(i, j, k + GRID)) != 0 ){
                            cerr << "Values are not the same" << endl;
                            cout << "Particle: " << p->getPosition() << endl;
                            cout << "Vector: " << Vector3f(i, j, k) << endl;
                        }
                        g.p[3] = Vector3f(p->getPosition());
                        g.val[3] = p->getDensity();

                        p = grid[x1][y2][z1];
                        if(p->getDistance(Vector3f(i, j + GRID, k)) != 0 ){
                            cerr << "Values are not the same" << endl;
                            cout << "Particle: " << p->getPosition() << endl;
                            cout << "Vector: " << Vector3f(i, j, k) << endl;
                        }
                        g.p[4] = Vector3f(p->getPosition());
                        g.val[4] = p->getDensity();

                        p = grid[x2][y2][z1];
                        if(p->getDistance(Vector3f(i + GRID, j + GRID, k)) != 0 ){
                            cerr << "Values are not the same" << endl;
                            cout << "Particle: " << p->getPosition() << endl;
                            cout << "Vector: " << Vector3f(i, j, k) << endl;
                        }
                        g.p[5] = Vector3f(p->getPosition());
                        g.val[5] = p->getDensity();

                        p = grid[x2][y2][z2];
                        if(p->getDistance(Vector3f(i+ GRID, j + GRID, k + GRID)) != 0 ){
                            cerr << "Values are not the same" << endl;
                            cout << "Particle: " << p->getPosition() << endl;
                            cout << "Vector: " << Vector3f(i, j, k) << endl;
                        }
                        g.p[6] = Vector3f(p->getPosition());
                        g.val[6] = p->getDensity();

                        p = grid[x1][y2][z2];
                        if(p->getDistance(Vector3f(i, j + GRID, k+ GRID)) != 0 ){
                            cerr << "Values are not the same" << endl;
                            cout << "Particle: " << p->getPosition() << endl;
                            cout << "Vector: " << Vector3f(i, j, k) << endl;
                        }
                        g.p[7] = Vector3f(p->getPosition());
                        g.val[7] = p->getDensity();



                        // glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
                        // glBegin(GL_QUADS); //bottom
                        //     glColor4f(0.1, 0.8, 0.6, 0.8);//235/255.0, 244/255.0, 250/255.0, 1.0);
                        //     glVertex3d(convert(g.p[0].x(), WIDTH), convert(g.p[0].y(), HEIGHT), convert(g.p[0].z(), LENGTH));
                        //     glVertex3d(convert(g.p[3].x(), WIDTH), convert(g.p[3].y(), HEIGHT), convert(g.p[3].z(), LENGTH));
                        //     glVertex3d(convert(g.p[7].x(), WIDTH), convert(g.p[7].y(), HEIGHT), convert(g.p[7].z(), LENGTH));
                        //     glVertex3d(convert(g.p[4].x(), WIDTH), convert(g.p[4].y(), HEIGHT), convert(g.p[4].z(), LENGTH));
                        // glEnd();

                        // glBegin(GL_QUADS); //bottom
                        //     glVertex3d(convert(g.p[4].x(), WIDTH), convert(g.p[4].y(), HEIGHT), convert(g.p[4].z(), LENGTH));
                        //     glVertex3d(convert(g.p[5].x(), WIDTH), convert(g.p[5].y(), HEIGHT), convert(g.p[5].z(), LENGTH));
                        //     glVertex3d(convert(g.p[6].x(), WIDTH), convert(g.p[6].y(), HEIGHT), convert(g.p[6].z(), LENGTH));
                        //     glVertex3d(convert(g.p[7].x(), WIDTH), convert(g.p[7].y(), HEIGHT), convert(g.p[7].z(), LENGTH));
                        // glEnd();
                        // glPolygonMode( GL_FRONT_AND_BACK, GL_FILL);

                        cubes->polygonise(g, totalDens/particles->size());



                    }
                }
            }
        }
        //second iteration of particles and only their neighbors
        #pragma omp parallel for
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

                //Pressure
                Vector3f rij = tempParticle->getPosition() - position;
                Vector3f kernDerive = particle->getKernDerive(dist, rij);
                double pressureK = tempParticle->calcPressure();
                pressureForce += tempMass * (pressureJ + pressureK) / (2 * tempDens) * kernDerive;

                //Viscosity
                double kernSecond = particle->getKernSecond(dist);
                viscosityForce += (tempVel - velocity) * tempMass / tempDens * kernSecond;

                //Surface Tension
                surfaceNormal += tempMass / tempDens * kernDerive;
                colorField += tempMass / tempDens * kernSecond;
            }
            // cout << surfaceNormal << endl;
            pressureForce *= -1;
            viscosityForce *= VISC;

            Vector3f surfaceTension = -TENSION * colorField * surfaceNormal / surfaceNormal.norm();
            Vector3f gravityForce(0, particle->getDensity() * GRAVITY, 0);
            // cout << "pForce: " << pressureForce << endl;
            // //cout << "dens: " << particle->getDensity() << endl;
            //  cout << "glForce: " << gravityForce << endl;
            //  cout << "vForce: " << viscosityForce << endl;

            //Update next position
            Vector3f totalForce = gravityForce + pressureForce + viscosityForce;
            if (surfaceNormal.norm() >= 0.00001)
                totalForce += surfaceTension;
            //cout << "totalForce: " << totalForce << endl;
            Vector3f acceleration = totalForce/particle->getDensity();
            //cout << "1. " << particle->getVelocity() << endl;
            Vector3f newVelocity = velocity + DELTAT * acceleration;  //maybe implement some kind of terminal velocity?
            //cout << "2. " << velocity << endl;
            Vector3f newPosition = position + DELTAT * newVelocity;

            //Boundary check next position
            // int c = 0.1;
            // double dist = abs(position.x() - LEFT - EPSILON);
            // if(dist < EPSILON ){
            //     Vector3f normal(1, 0, 0);
            //     double incidentAngle = velocity.dot(normal);
            //     if(incidentAngle < 0){
            //         //double reflectAngle = incidentAngle * 1.1;
            //         newVelocity = newVelocity * -0.3;
            //         newVelocity = newVelocity + c * (2 * EPSILON - dist) * normal;
            //         newPosition = position + DELTAT * newVelocity;
            //     }
            // }
            // dist = abs(position.x() - RIGHT + EPSILON);
            // if(dist < EPSILON ){
            //     Vector3f normal(-1, 0, 0);
            //     double incidentAngle = velocity.dot(normal);
            //     if(incidentAngle < 0){
            //         //double reflectAngle = incidentAngle * 1.1;
            //         newVelocity = velocity * -0.3;
            //         newVelocity = newVelocity + c * (2 * EPSILON - dist) * normal;
            //         newPosition = position + DELTAT * newVelocity;
            //     }
            // }
            // dist = abs(position.y() - BOTTOM - EPSILON);
            // if(dist < EPSILON ){
            //     Vector3f normal(0, 1, 0);
            //     double incidentAngle = velocity.dot(normal);
            //     if(incidentAngle < 0){
            //         //double reflectAngle = incidentAngle * 1.1;
            //         newVelocity = velocity * -0.3;
            //         newVelocity = newVelocity + c * (2 * EPSILON - dist) * normal;
            //         newPosition = position + DELTAT * newVelocity;
            //     }
            // }
            // dist = abs(position.y() - TOP + EPSILON);
            // if(dist < EPSILON ){
            //     Vector3f normal(0, -1, 0);
            //     double incidentAngle = velocity.dot(normal);
            //     if(incidentAngle < 0){
            //         //double reflectAngle = incidentAngle * 1.1;
            //         newVelocity = velocity * -0.3;
            //         newVelocity = newVelocity + c * (2 * EPSILON - dist) * normal;
            //         newPosition = position + DELTAT * newVelocity;
            //     }
            // }
            // dist = abs(position.z() - FRONT + EPSILON);
            // if(dist < EPSILON ){
            //     Vector3f normal(0, 0, -1);
            //     double incidentAngle = velocity.dot(normal);
            //     if(incidentAngle < 0){
            //         //double reflectAngle = incidentAngle * 1.1;
            //         newVelocity = velocity * -0.3;
            //         newVelocity = newVelocity + c * (2 * EPSILON - dist) * normal;
            //         newPosition = position + DELTAT * newVelocity;
            //     }
            // }
            // dist = abs(position.z() - BACK - EPSILON);
            // if(dist < EPSILON ){
            //     Vector3f normal(0, 0, 1);
            //     double incidentAngle = velocity.dot(normal);
            //     if(incidentAngle < 0){
            //         //double reflectAngle = incidentAngle * 1.1;
            //         newVelocity = velocity * -0.3;
            //         newVelocity = newVelocity + c * (2 * EPSILON - dist) * normal;
            //         newPosition = position + DELTAT * newVelocity;
            //     }
            // }


            // bool bounce = false;
            //  while((newPosition.x() - RADIUS <= LEFT) || (newPosition.y() - RADIUS <= BOTTOM) || (newPosition.x() + RADIUS >= RIGHT) || (newPosition.y() + RADIUS >= TOP) || (newPosition.z() + RADIUS <= BACK) || (newPosition.z() - RADIUS >= FRONT)){
            //     bounce = true;
                // newVelocity *= 0.9;
            //     newPosition = position + DELTAT * newVelocity;
            //  }
             // if(bounce) velocity *= -1;

            // bool bound = false;
            double c = -0.3;
            // Vector3f collNorm = Vector3f::Zero();
            //if(t < 100 ){
                if(newPosition.y() - EPSILON < BOTTOM){
                    // double boundTime = (position.y() - (BOTTOM)) / velocity.norm();
                    // cout << boundTime << endl;  //maybe boundTime = min(boundTime, DELTA);
                    // Vector3f collision = position + boundTime * velocity;
                    // // cout << "collision: " << collision << endl;
                    // collNorm << 0, 1, 0;
                    // double penDist = (newPosition - collision).norm();
                    // //cout << penDist << endl;
                    // // cout << "1. " << newPosition << endl;
                    // newPosition = newPosition + penDist * collNorm;
                    // bound = true;
                    // velocity = velocity - ( 1 + cr) * (velocity.dot(collNorm) * collNorm);
                    // cout << "2. " << newPosition << endl;
                    // bound = true;
                    // velocity *= -0.3;
                    // while(newPosition.y() - EPSILON < BOTTOM){
                    //     newPosition = newPosition + velocity * DELTAT;
                    // }
                    double d = abs(newPosition.y() - BOTTOM);
                    newPosition = Vector3f(newPosition.x(), BOTTOM + EPSILON, newPosition.z());
                    newVelocity = Vector3f(newVelocity.x(), newVelocity.y() * c, newVelocity.z());
                }
                if(newPosition.y() + EPSILON > TOP){
                    // double boundTime = (-position.y() + (TOP)) / velocity.norm();
                    // Vector3f collision = position + boundTime * velocity;
                    // collNorm << 0, -1, 0;
                    // double penDist = (newPosition - collision).norm();
                    // newPosition = newPosition + penDist * collNorm;
                    // // velocity = velocity - ( 1 + cr) * (velocity.dot(collNorm) * collNorm);
                    // bound = true;
                    // velocity *= -0.3;
                    // while(newPosition.y() + EPSILON > TOP){
                    //     newPosition = newPosition + velocity * DELTAT;
                    // }
                    double d = abs(newPosition.y() - TOP);
                    newPosition = Vector3f(newPosition.x(), newPosition.y() - d, newPosition.z());
                    newVelocity = Vector3f(newVelocity.x(), newVelocity.y() * c, newVelocity.z());
                }
                if(newPosition.x() - EPSILON < LEFT){
                    // double boundTime = (position.x() - (LEFT)) / velocity.norm();
                    // Vector3f collision = position + boundTime * velocity;
                    // collNorm << 1, 0, 0;
                    // double penDist = (newPosition - collision).norm();
                    // newPosition = newPosition + penDist * collNorm;
                    // // velocity = velocity - ( 1 + cr) * (velocity.dot(collNorm) * collNorm);
                    // bound = true;
                    // velocity *= -0.3;
                    // while(newPosition.x() + EPSILON < LEFT){
                    //     newPosition = newPosition + velocity * DELTAT;
                    // }
                    double d = abs(newPosition.x() - LEFT);
                    newPosition = Vector3f(newPosition.x() + d, newPosition.y(), newPosition.z());
                    newVelocity = Vector3f(newVelocity.x() * c, newVelocity.y(), newVelocity.z());
                }
                if(newPosition.x() + EPSILON > RIGHT){
                    // double boundTime = (-position.x() + (RIGHT)) / velocity.norm();
                    // Vector3f collision = position + boundTime * velocity;
                    // collNorm << -1, 0, 0;
                    // double penDist = (newPosition - collision).norm();
                    // newPosition = newPosition + penDist * collNorm;
                    // // velocity = velocity - ( 1 + cr) * (velocity.dot(collNorm) * collNorm);
                    // bound = true;
                    // velocity *= -0.3;
                    // while(newPosition.x() - EPSILON > RIGHT){
                    //     newPosition = newPosition + velocity * DELTAT;
                    // }
                    double d = abs(newPosition.x() - RIGHT);
                    newPosition = Vector3f(newPosition.x() - d, newPosition.y(), newPosition.z());
                    newVelocity = Vector3f(newVelocity.x() * c, newVelocity.y(), newVelocity.z());
                }
                if(newPosition.z() - EPSILON < BACK){
                    // double boundTime = (position.z() - (BACK)) / velocity.norm();
                    // Vector3f collision = position + boundTime * velocity;
                    // collNorm << 0, 0, 1;
                    // double penDist = (newPosition - collision).norm();
                    // newPosition = newPosition + penDist * collNorm;
                    // // velocity = velocity - ( 1 + cr) * (velocity.dot(collNorm) * collNorm);
                    // bound = true;
                    // velocity *= -0.3;
                    // while(newPosition.z() - EPSILON < BACK){
                    //     newPosition = newPosition + velocity * DELTAT;
                    // }
                    double d = abs(newPosition.z() - BACK);
                    newPosition = Vector3f(newPosition.x(), newPosition.y(), newPosition.z() + d);
                    newVelocity = Vector3f(newVelocity.x(), newVelocity.y(), newVelocity.z() * c);
                }
                if(newPosition.z() + EPSILON > FRONT){
                    // double boundTime = (-position.z() + (FRONT)) / velocity.norm();
                    // Vector3f collision = position + boundTime * velocity;
                    // collNorm << 0, 0, -1;
                    // double penDist = (newPosition - collision).norm();
                    // newPosition = newPosition + penDist * collNorm;
                    // // velocity = velocity - ( 1 + cr) * (velocity.dot(collNorm) * collNorm);
                    // bound = true;
                    // velocity *= -0.3;
                    // while(newPosition.z() + EPSILON > FRONT){
                    //     newPosition = newPosition + velocity * DELTAT;
                    // }
                    double d = abs(newPosition.z() - FRONT);
                    newPosition = Vector3f(newPosition.x(), newPosition.y(), newPosition.z() - d);
                    newVelocity = Vector3f(newVelocity.x(), newVelocity.y(), newVelocity.z() * c);
                }
                // if(bound) {
                //     newVelocity *= -0.3;
                // }
            //}
            particle->setPosition(newPosition);
            particle->setVelocity(newVelocity);
        }
        saveImage(t);

        glFlush();
        glutSwapBuffers();
        glPopMatrix();
    }
}


void Scene::drawBoundaries(){
        glDisable(GL_LIGHTING);
        glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glBegin(GL_QUADS); //bottom
            glColor4f(0.1, 0.8, 0.6, 0.8);//235/255.0, 244/255.0, 250/255.0, 1.0);
            glVertex3d(convert(LEFT, WIDTH), convert(BOTTOM, HEIGHT), convert(FRONT, LENGTH));
            glVertex3d(convert(RIGHT, WIDTH), convert(BOTTOM, HEIGHT), convert(FRONT, LENGTH));
            glVertex3d(convert(RIGHT, WIDTH), convert(BOTTOM, HEIGHT), convert(BACK, LENGTH));
            glVertex3d(convert(LEFT, WIDTH), convert(BOTTOM, HEIGHT), convert(BACK, LENGTH));
        glEnd();
        glBegin(GL_QUADS); //top
            glVertex3d(convert(LEFT, WIDTH), convert(TOP, HEIGHT), convert(FRONT, LENGTH));
            glVertex3d(convert(RIGHT, WIDTH), convert(TOP, HEIGHT), convert(FRONT, LENGTH));
            glVertex3d(convert(RIGHT, WIDTH), convert(TOP, HEIGHT), convert(BACK, LENGTH));
            glVertex3d(convert(LEFT, WIDTH), convert(TOP, HEIGHT), convert(BACK, LENGTH));
        glEnd();
        glBegin(GL_QUADS); //left
                // glColor3f(0.0, 1.0, 0.0);
            glVertex3d(convert(LEFT, WIDTH), convert(BOTTOM, HEIGHT), convert(FRONT, LENGTH));
            glVertex3d(convert(LEFT, WIDTH), convert(BOTTOM, HEIGHT), convert(BACK, LENGTH));
            glVertex3d(convert(LEFT, WIDTH), convert(TOP, HEIGHT), convert(BACK, LENGTH));
            glVertex3d(convert(LEFT, WIDTH), convert(TOP, HEIGHT), convert(FRONT, LENGTH));
        glEnd();
        glBegin(GL_QUADS); //right
            glVertex3d(convert(RIGHT, WIDTH), convert(BOTTOM, HEIGHT), convert(FRONT, LENGTH));
            glVertex3d(convert(RIGHT, WIDTH), convert(TOP, HEIGHT), convert(FRONT, LENGTH));
            glVertex3d(convert(RIGHT, WIDTH), convert(TOP, HEIGHT), convert(BACK, LENGTH));
            glVertex3d(convert(RIGHT, WIDTH), convert(BOTTOM, HEIGHT), convert(BACK, LENGTH));
        glEnd();
        glBegin(GL_QUADS); //back
            glVertex3d(convert(LEFT, WIDTH), convert(BOTTOM, HEIGHT), convert(BACK, LENGTH));
            glVertex3d(convert(RIGHT, WIDTH), convert(BOTTOM, HEIGHT), convert(BACK, LENGTH));
            glVertex3d(convert(RIGHT, WIDTH), convert(TOP, HEIGHT), convert(BACK, LENGTH));
            glVertex3d(convert(LEFT, WIDTH), convert(TOP, HEIGHT), convert(BACK, LENGTH));
        glEnd();
        glBegin(GL_QUADS); //front
            glVertex3d(convert(LEFT, WIDTH), convert(BOTTOM, HEIGHT), convert(FRONT, LENGTH));
            glVertex3d(convert(RIGHT, WIDTH), convert(BOTTOM, HEIGHT), convert(FRONT, LENGTH));
            glVertex3d(convert(RIGHT, WIDTH), convert(TOP, HEIGHT), convert(FRONT, LENGTH));
            glVertex3d(convert(LEFT, WIDTH), convert(TOP, HEIGHT), convert(FRONT, LENGTH));
        glEnd();
        glEnable(GL_LIGHTING);
        glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
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

Particle::Particle(){}

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
