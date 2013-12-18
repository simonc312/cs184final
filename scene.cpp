//ffmpg
#include "fluid.h"

//PROBLEMS: particles collect in the back, maybe something to do with using bottom left FRONT vertex as reference; those particles lose getting referred to
//improve marching cubes


Scene::Scene(int p, double t, double s, int m){
    maxParts = p;
    timeStep = t;
    step = s;
    march = m;
    particles = new vector<Particle *>();
    film = new Film(WIDTH, HEIGHT);
    cubes = new Cubes();

    int x = GRIDX + 1;
    int y = GRIDY + 1;
    int z = GRIDZ + 1;
    grids = new GRIDCELL**[x];
    for(int i=0;i<x;i++){
        grids[i] = new GRIDCELL *[y];
    }
    for(int i = 0; i < x; i ++){
        for(int j = 0; j < y; j++){
            grids[i][j] = new GRIDCELL[z];
        }
    }
    // cout << "HERE-2" << endl;

    init();
}

void Scene::init(){
    srand(time(NULL));

    // for(int i = 0; i < GRIDX + 1; i++){
    //     for(int j = 0; j < GRIDY + 1; j++){
    //         for(int k = 0; k < GRIDZ + 1; k++){
    //             Particle * particle = new Particle(0, Vector3f(i * GRID + LEFT, j * GRID + BOTTOM, k * GRID + BACK), Vector3f(0, 0, 0));
    //             int gridX = floor((particle->getPosition().x() - LEFT)/GRID);
    //             int gridY = floor((particle->getPosition().y() - BOTTOM)/GRID);
    //             int gridZ = abs(floor((particle->getPosition().z() - FRONT)/GRID));
    //             // grids[i][j][k].vParts
    //             grids[gridX][gridY][gridZ].particles.push_back(particle);
    //             particles->push_back(particle);
    //         }
    //     }
    // }
    // cout << "HERE-1.5" << endl;
    for(int i= 0 ; i < maxParts; i ++){
            double x = fRand(375, 415);
            double y = fRand(475, 625);
            double z = fRand(-375, -415);
            Vector3f pos(x, y, z);
            // cout << "HERE-1" << endl;
            int gridX = floor((x - LEFT)/GRID);
            int gridY = floor((y - BOTTOM)/GRID);
            int gridZ = floor(abs((z - FRONT)/GRID));
            Particle *p = new Particle(MASS, pos, Vector3f(0, 0, 0));
            p->setGridPosition(Vector3f(gridX, gridY, gridZ));
            grids[gridX][gridY][gridZ].particles.push_back(p);
            particles->push_back(p);
    }
    // cout << particles->size() << endl;
}



void Scene::render(){
    for(int t = 0; t < timeStep; t++){  //for every timestep
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);               // clear the color buffer
        glMatrixMode(GL_MODELVIEW);                   // indicate we are specifying camera transformations
        glLoadIdentity();
        // gluLookAt(0.0, 0.15, -1.0, 0.0, 0.15, -10.0, 0.0, 1.0, 0.0);
        // gluLookAt(0.0, 0.0, 0.65, 0.0, 0.0, -10.0, 0.0, 1.0, 0.0);
        // gluLookAt(0.0, 0.6, -1.25, .0, 0.0, -2.1, 0.0, 1.0, 0.0);
        // gluLookAt(0.0, 0.0, -1.0,0.0, 0.0, -10.0, 0.0, 1.0, 0.0);
        gluLookAt(0.0, 0.6, -0.8, 0.0, -1.0, -10.0, 0.0, 1.0, 0.0);

        //Draw the boundaries
        drawBoundaries();

        //density calculations
        vector<vector<Particle * > > neighbors;
        // double totalDens = 0;
        // cout << "HERE0" << endl;
        for(int i = 0; i < particles->size(); i++){  //for every particle
            //cout << omp_get_num_threads() << endl;
            Particle *particle = particles->at(i);
            double density = MASS;
            vector<Particle *> findNeighs;
            Vector3f gridPos = particle->getGridPosition();
            Vector3f BBmin = gridPos - Vector3f(1, 1, 1);
            Vector3f BBmax = gridPos + Vector3f(1, 1, 1);

            for(int y = BBmin.y(); y <=BBmax.y(); y++){
                for(int x = BBmin.x(); x <=BBmax.x(); x++){
                    for(int z = BBmin.z(); z <= BBmax.z(); z++){

                        if(x >= 0 && y >= 0 && z >= 0 && x <= GRIDX && y <= GRIDY && z <= GRIDZ){
                            // cout << "HERE " << endl;
                            vector<Particle *> gridParticles = grids[x][y][z].particles;

                            for(int j = 0; j < gridParticles.size(); j++){
                                Particle *tempParticle = gridParticles.at(j);
                                double dist = particle->getDistance(*tempParticle);

                                if(dist <= H && dist != 0){
                                    double kern = particle->getKernel(dist);
                                    density += tempParticle->getMass() * kern;
                                    findNeighs.push_back(tempParticle);
                                }
                            }
                        }
                    }
                }
            }
            particle->setDensity(density);
            #pragma omp critical
            {
            neighbors.push_back(findNeighs);
            // totalDens += density;
            }
        }
                    // cout << "HERE1" << endl;
        // cout << totalDens/particles->size() << endl;

        //Marching Cubes
        if(march == 1){

            int x = GRIDX + 1;
            int y = GRIDY + 1;
            int z = GRIDZ + 1;
            // cout << "HERE" << endl;
            // cout << x << " " << y << " " << z << endl;

            Particle **** grid;
            grid = new Particle***[x];
            for(int i=0;i<x;i++){
                grid[i] = new Particle **[y];
            }
            for(int i = 0; i < x; i ++){
                for(int j = 0; j < y; j++){
                    grid[i][j] = new Particle*[z];
                }
            }
            // cout << "HERE2" << endl;

            double totalDens = 0;
            int count = 0;
            for(int i = 0; i < x; i++){
                for(int j = 0; j < y; j++){
                    for(int k = 0; k < z; k++){
                        Particle* particle = new Particle(MASS, Vector3f(i  * GRID + LEFT, j * GRID + BOTTOM, k * GRID + BACK), Vector3f(0, 0, 0));
                        count ++;
                        double density = MASS;
                        Vector3f gridPos((particle->getPosition().x() - LEFT)/GRID, (particle->getPosition().y() - BOTTOM)/GRID, abs((particle->getPosition().z() - FRONT)/GRID));
                        // cout << gridPos << endl << endl;
                        Vector3f BBmin = gridPos - Vector3f(1, 1, 1);
                        Vector3f BBmax = gridPos + Vector3f(1, 1, 1);

                        for(int b = BBmin.y(); b <=BBmax.y(); b++){
                            for(int a = BBmin.x(); a <=BBmax.x(); a++){
                                for(int c = BBmin.z(); c <= BBmax.z(); c++){
                                    // cout << x  << " " << y << " " << z << endl;
                                    if(b >= 0 && a >= 0 && c >= 0 && a <= GRIDX && b <= GRIDY && c <= GRIDZ){
                                        // cout << "HERE " << endl;
                                        vector<Particle *> gridParticles = grids[a][b][c].particles;

                                        for(int l = 0; l < gridParticles.size(); l++){
                                            Particle *tempParticle = gridParticles.at(l);
                                            double dist = particle->getDistance(*tempParticle);

                                            if(dist <= H &&  dist != 0){
                                                double kern = particle->getKernel(dist);
                                                                                                // cout << kern << endl;
                                                density += tempParticle->getMass() * kern;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        particle->setDensity(density);
                        totalDens += density;
                        // cout << totalDens/count << endl;
                        grid[i][j][k] = particle;
                    }
                }
            }



            for(int i = LEFT; i < RIGHT; i += GRID){
                for(int j = BOTTOM; j < TOP; j += GRID){
                    for(int k = BACK; k < FRONT; k +=GRID){
                        // cout << "HERE" << endl;
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

                        cubes->polygonise(g, totalDens/count);



                    }
                }
            }
        }


        //second iteration of particles and only their neighbors
        //#pragma omp parallel for
        for(int i = 0; i < particles->size(); i++){
            Particle *particle = particles->at(i);

            // if(particle->getMass() == 0){ //Particle mass is 0 if it is a vertex particle; we don't want to calculate forces for it
            //     continue;
            // }

            Vector3f position = particle->getPosition();
            Vector3f velocity = particle->getVelocity();
            // cout << position << endl << endl;

            //http://stackoverflow.com/questions/17565664/gluproject-and-2d-display
            //Render particle
            if(march == 0){
                GLdouble posX, posY, posZ;//3D point
                posX=convert(position.x(), WIDTH);
                posY=convert(position.y(), HEIGHT);
                posZ=convert(position.z(), LENGTH);

                glPushMatrix();
                    glTranslated(posX, posY, posZ);
                    glutSolidSphere(SRADIUS, 10, 10);
                glPopMatrix();
            }
            //Force calculations
            Vector3f viscosityForce = Vector3f::Zero();
            Vector3f pressureForce = Vector3f::Zero();
            Vector3f surfaceNormal = Vector3f::Zero();
            double colorField = 0;
            double pressureJ = particle->calcPressure();

            vector<Particle * > curNeighs = neighbors[i];
            // cout << curNeighs.size() << endl;
            // cout << "particle " << i << " num neighbors" << curNeighs.size() << endl;
            for(int j = 0; j < curNeighs.size(); j++){//currNeighs.size(); j++){

                Particle *tempParticle = curNeighs[j];// currNeighs[j]->p;
                double tempMass = tempParticle->getMass();
                double tempDens = tempParticle->getDensity();
                Vector3f tempVel = tempParticle->getVelocity();
                double dist = particle->getDistance(*tempParticle);
                double kern = particle->getKernel(dist);
                // cout << "dist: " << dist << "\npos: " << position << "\ntempPos: " << tempParticle->getPosition() << endl;

                if(dist != 0){
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
            }
            // cout << surfaceNormal << endl;
            pressureForce *= -1;
            viscosityForce *= VISC;

            Vector3f surfaceTension = -TENSION * colorField * surfaceNormal / surfaceNormal.norm();
            Vector3f gravityForce(0, particle->getDensity() * GRAVITY, 0);

            //Update next position
            Vector3f totalForce = gravityForce + pressureForce + viscosityForce;
            //cout << "vForce: " << viscosityForce << endl;
            if (surfaceNormal.norm() >= 0.00001)
                totalForce += surfaceTension;
            //cout << "totalForce: " << totalForce << endl;
            Vector3f acceleration = totalForce/particle->getDensity();
            //cout << "1. " << particle->getVelocity() << endl;
            Vector3f newVelocity = velocity + DELTAT * acceleration;  //maybe implement some kind of terminal velocity?
            //cout << "2. " << velocity << endl;
            Vector3f newPosition = position + DELTAT * newVelocity;

            // cout << "HERE2" << endl;



            //Boundary checks
            double c = -0.3;
            if(newPosition.y() - EPSILON < BOTTOM){
                double d = abs(newPosition.y() - BOTTOM);
                newPosition = Vector3f(newPosition.x(), BOTTOM + EPSILON, newPosition.z());
                newVelocity = Vector3f(newVelocity.x(), newVelocity.y() * c, newVelocity.z());
            }
            if(newPosition.y() + EPSILON > TOP){
                double d = abs(newPosition.y() - TOP);
                newPosition = Vector3f(newPosition.x(), newPosition.y() - d, newPosition.z());
                newVelocity = Vector3f(newVelocity.x(), newVelocity.y() * c, newVelocity.z());
            }
            if(newPosition.x() - EPSILON < LEFT){
                double d = abs(newPosition.x() - LEFT);
                newPosition = Vector3f(newPosition.x() + d, newPosition.y(), newPosition.z());
                newVelocity = Vector3f(newVelocity.x() * c, newVelocity.y(), newVelocity.z());
            }
            if(newPosition.x() + EPSILON > RIGHT){
                double d = abs(newPosition.x() - RIGHT);
                newPosition = Vector3f(newPosition.x() - d, newPosition.y(), newPosition.z());
                newVelocity = Vector3f(newVelocity.x() * c, newVelocity.y(), newVelocity.z());
            }
            if(newPosition.z() - EPSILON < BACK){
                double d = abs(newPosition.z() - BACK);
                newPosition = Vector3f(newPosition.x(), newPosition.y(), newPosition.z() + d);
                newVelocity = Vector3f(newVelocity.x(), newVelocity.y(), newVelocity.z() * c);
            }
            if(newPosition.z() + EPSILON > FRONT){
                double d = abs(newPosition.z() - FRONT);
                newPosition = Vector3f(newPosition.x(), newPosition.y(), newPosition.z() - d);
                newVelocity = Vector3f(newVelocity.x(), newVelocity.y(), newVelocity.z() * c);
            }
                        Vector3f gridPos((newPosition.x() - LEFT)/GRID, (newPosition.y() - BOTTOM)/GRID, abs((newPosition.z() - FRONT)/GRID));

            if(!(gridPos.x() >= 0 && gridPos.y() >= 0 && gridPos.z() >= 0 && gridPos.x() <= GRIDX && gridPos.y() <= GRIDY && gridPos.z() <= GRIDZ)) {
                cout << "gridPos: " << gridPos << endl;
                cout << "surface tension: " << surfaceTension << endl;
                cout << "pressure: " << pressureForce << endl;
                exit(0);
            }

            particle->setPosition(newPosition);
            particle->setVelocity(newVelocity);
        }


        saveImage(t);

        glFlush();
        glutSwapBuffers();
        glPopMatrix();

        // cout << "HERE3" << endl;

        //Clear particle in current grids
        for(int i = 0; i < GRIDX; i++){
            for(int j = 0; j < GRIDY; j++){
                for(int k = 0; k < GRIDZ; k++){
                    grids[i][j][k].particles.clear();
                    // vector<Particle * > particles = grids[i][j][k].particles;
                    // int l = 0;
                    // for(vector<Particle *>::iterator it = particles.begin(); it != particles.end();) {   //Delete only particles with mass from GRIDCELL; particles without mass are vertex particles
                    //     if(particles[l]->getMass() != 0){
                    //         it = particles.erase(it);
                    //     } else{
                    //         ++it;
                    //         ++l;
                    //     }
                    // }
                }
            }
        }
        //cout << "HERE1" << endl;
        // cout << "HERE4" << endl;

        //Update particles in each grid
        for(int i = 0; i < particles->size(); i++){
            // cout << "i: " << i << " size: " << particles->size() << endl;
            // cout << "HERE1" << endl;
            Particle *particle = particles->at(i);
            // if(particle->getMass() != 0){
            Vector3f position = particle->getPosition();
            // cout << position << endl << endl;
            int gridX = floor((position.x() - LEFT)/GRID);
            int gridY = floor((position.y() - BOTTOM)/GRID);
            int gridZ = abs(floor((position.z() - FRONT)/GRID));
            //cout << "HERE2" << endl;
            particle->setGridPosition(Vector3f(gridX, gridY, gridZ));
            // cout << "HERE3" << endl;
            // cout << gridX << " " << gridY << " " << gridZ << endl;
            grids[gridX][gridY][gridZ].particles.push_back(particle);
            //cout << "HERE4" << endl;
            // }
        }
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

void Particle::setGridPosition(Vector3f p){
    gridPos = p;
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
    if (r == 0){
        r = 0.001;
    }
    return c * (rij / r) * (pow(H - r, 2));
}

double Particle::getKernSecond(double r){
    double c = 45 / (M_PI * pow((double)H, 6.0));
    return c * (H - r);
}

Vector3f Particle::getGridPosition(){
    return gridPos;
}
