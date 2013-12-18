#include "fluid.h"



Viewport    viewport;
Scene *s;
int rotAngleY = 0;
int rotAngleX = 0;
double transX = 0.0;
double transY = 0.0;
double transZ = -4.5;

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
}

//****************************************************
// Simple init function
//****************************************************
void initScene(){
    glClearColor(0.9f, 0.9f, 0.9f, 0.0f); // Clear to black, fully transparent
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



    // GLfloat diffuse0[]={0.2, 0.2, 0.2, 1.0};
    // GLfloat ambient0[]={0.1, 0.1, 0.1, 1.0};
    // GLfloat specular0[]={0.6, 0.6, 0.6, 1.0};
    glEnable(GL_LIGHT1);
    GLfloat light1_pos[]={0.0, 1.0, 0.0, 1.0};

    glLightfv(GL_LIGHT1, GL_POSITION, light1_pos);
    glLightfv(GL_LIGHT1, GL_AMBIENT, ambient0);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse0);
    glLightfv(GL_LIGHT1, GL_SPECULAR, specular0);
  glPushMatrix(); // put current matrix on stack

    GLfloat ambient[] = {0.2, 0.2, 0.2, 1.0};
    GLfloat diffuse[] = {93/255.0, 170/255.0, 231/255.0, 0.9};
    GLfloat specular[] = {1.0, 1.0, 1.0, 1.0};
    GLfloat shine = 100.0;

    glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT, GL_SHININESS, shine);
    glPushMatrix(); // put current matrix on stack

      glTranslatef(transX, 0.0f, transZ);
      glTranslatef(0.0f, transY, transZ);

      glRotatef(rotAngleY, 0.0f, 1.0f, 0.0f);
      glRotatef(rotAngleX, 1.0f, 0.0f, 0.0f);

      s->render();
    glFlush();
    glutSwapBuffers();
    glPopMatrix();

}

void myKeyboard(unsigned char key, int x, int y){
  if(key == ' ')
    exit(0);
  if(key == '+')
    transZ += 1.5;
  if(key == '-')
    transZ -= 1.5;
  myDisplay();
}

void mySpecialKey(int key, int x, int y){
  int mod = glutGetModifiers();
  double transAmount = 0.1;
  if(mod == GLUT_ACTIVE_SHIFT){
    if(key == GLUT_KEY_LEFT) transX -= transAmount;
    if(key == GLUT_KEY_RIGHT) transX += transAmount;
    if(key == GLUT_KEY_UP) transY += transAmount;
    if(key == GLUT_KEY_DOWN) transY -= transAmount;
  }
  else{ //weird issues with up down transform where it'll go back on itself, might be openGL
    double rotateAmount = 10.0;
    if(key == GLUT_KEY_LEFT) rotAngleY -= rotateAmount;
    if(key == GLUT_KEY_RIGHT) rotAngleY+= rotateAmount;
    if(key == GLUT_KEY_UP) rotAngleX-=rotateAmount;
    if(key == GLUT_KEY_DOWN) rotAngleX+=rotateAmount;
  }
  myDisplay();
}

int main(int argc, char* argv[]){
    int maxParticles = atoi(argv[1]);
    double timeStep = atof(argv[2]);
    double stepSize = 0.1;
    int march = 0;
    march = atoi(argv[3]);
    s = new Scene(maxParticles, timeStep, stepSize, march);
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
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    GLfloat global_ambient[] = { 0.0f, 0.5f, 1.0f, 1.0f };
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHTING);

    glutDisplayFunc(myDisplay);               // function to run when its time to draw something
    glutReshapeFunc(myReshape);               // function to run when the window gets resized
    glutKeyboardFunc(myKeyboard);
      glutSpecialFunc(mySpecialKey);


    glutMainLoop();                           // infinite loop that will keep drawing and resizing
    // and whatever else

    return 0;
}