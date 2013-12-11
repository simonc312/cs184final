#include "fluid.h"



Viewport    viewport;
Scene *s;

//****************************************************
// reshape viewport if the window is resized
//****************************************************
void myReshape(int w, int h) {
  viewport.w = w;
  viewport.h = h;

  glViewport (0,0,viewport.w,viewport.h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective (60, (float)w/(float)h, 0.01f, (float)LENGTH);
  // glTranslatef(0.0, 100.0, 0.0);
  // glRotatef(90.0, 0.0, 1.0, 0.0);

  //glOrtho(0, WIDTH, 0, HEIGHT, -1.0, LENGTH);

  //glMatrixMode(GL_MODELVIEW);
  //glLoadIdentity();
          // gluLookAt (100.0, 0.0, -3.0, 0.0, 0.0, -5.0 , 1.0, 0.0, 0.0);

}

//****************************************************
// Simple init function
//****************************************************
void initScene(){
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f); // Clear to black, fully transparent
    myReshape(viewport.w,viewport.h);
}

//****************************************************
// function that does the actual drawing of stuff
//***************************************************
void myDisplay() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);               // clear the color buffer


  glMatrixMode(GL_MODELVIEW);                   // indicate we are specifying camera transformations
  glLoadIdentity();                     // make sure transformation is "zero'd"
  // Start drawing
   //gluLookAt (0.0, 0.0, 0.0, 0.0, 0.0, -1.0 , 1.0, 0.0, 0.0);

    GLfloat diffuse0[]={0.2, 0.2, 0.2, 1.0};
    GLfloat ambient0[]={0.1, 0.1, 0.1, 1.0};
    GLfloat specular0[]={0.6, 0.6, 0.6, 1.0};
    GLfloat light0_pos[]={1.0, 1.0, -0.5, 0.0};

    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0, GL_POSITION, light0_pos);
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient0);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse0);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specular0);

    glEnable(GL_LIGHT1);
    GLfloat light1_pos[]={-3.0, 0.0, 4.0, 1.0};

    glLightfv(GL_LIGHT1, GL_POSITION, light1_pos);
    glLightfv(GL_LIGHT1, GL_AMBIENT, ambient0);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse0);
    glLightfv(GL_LIGHT1, GL_SPECULAR, specular0);
  glPushMatrix(); // put current matrix on stack

    GLfloat ambient[] = {0.2, 0.2, 0.2, 0.6};
    GLfloat diffuse[] = {0.2, 0.5, 1.0, 0.6};
    GLfloat specular[] = {1.0, 1.0, 1.0, 0.6};
    GLfloat shine = 100.0;

    glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT, GL_SHININESS, shine);

    s->render();
    glFlush();
    glutSwapBuffers();
    glPopMatrix();

}

void myKeyboard(unsigned char key, int x, int y){
    if(key == ' ')
        exit(0);
}


int main(int argc, char* argv[]){
    int maxParticles = atoi(argv[1]);
    double timeStep = atof(argv[2]);
    double stepSize = 0.1;
    s = new Scene(maxParticles, timeStep, stepSize);
    glutInit(&argc, argv);

    //This tells glut to use a double-buffered window with red, green, and blue channels
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    // Initalize theviewport size
    viewport.w = WIDTH;
    viewport.h = HEIGHT;

    //The size and position of the window

    glutInitWindowSize(viewport.w, viewport.h);
    glutInitWindowPosition(0,0);
    glutCreateWindow("Fluids");

    initScene();                          // quick function to set up scene

    //default settings: fill, shading
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
    glShadeModel(GL_SMOOTH);

    GLfloat global_ambient[] = { 0.2f, 0.5f, 1.0f, 1.0f };
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHTING);

    glutDisplayFunc(myDisplay);               // function to run when its time to draw something
    glutReshapeFunc(myReshape);               // function to run when the window gets resized
    glutKeyboardFunc(myKeyboard);

    glutMainLoop();                           // infinite loop that will keep drawing and resizing
    // and whatever else

    return 0;
}