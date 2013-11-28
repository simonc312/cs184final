//Code from http://inst.eecs.berkeley.edu/~cs184/fa09/resources/sec_UsingFreeImage.pdf
#include "fluid.h"

using namespace std;
using namespace Eigen;

Film::Film(int w, int h){
    width = w;
    height = h;
    count = 0;
}

void Film::saveImage( vector<vector<Color> > m) {
    FreeImage_Initialise();
    FIBITMAP* bitmap = FreeImage_Allocate(width, height, BPP);
    RGBQUAD color;
    if (! bitmap )
        exit(1);
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            color.rgbRed = (double)m[i][j].r * 255.0;
            color.rgbGreen = (double)m[i][j].g * 255.0;
            color.rgbBlue = (double)m[i][j].b * 255.0;
            FreeImage_SetPixelColor(bitmap, i, j, &color);

        }
    }

    // get a string to concat with an int
    stringstream nameStream;
    nameStream << "test" << count << ".png";
    string name = nameStream.str();
    char *a = new char[name.size() + 1];
    a[name.size()] = 0;
    memcpy(a, name.c_str(), name.size());
    if(FreeImage_Save(FIF_PNG, bitmap, a, 0))
        cout << "Image successfully saved! " << endl ;
    FreeImage_DeInitialise(); //Cleanup !
    count ++;
}