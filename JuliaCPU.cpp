/*#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <chrono>
#include "lodepng.h"

namespace lodepng {
unsigned encode(const std::string& filename,
const unsigned char* image,
unsigned w, unsigned h,
LodePNGColorType colortype,
unsigned bitdepth);
}

void Julia(std::vector<double>X, std::vector<double>Y, std::vector<std::vector<double>>image, unsigned char * ArImg, std::complex< double >c_,int N){
    for(double j=0;j<N;++j){
        for(double i=0;i<N;++i){
            int iters = 0;                     //iteration count 
            int count = j*N+i;                 //ArImg index
            std::complex<double> z(X[i],Y[j]); //definign the complex number
            while(iters<250){
                z = std::pow(z,2) + c_;
                if(std::norm(z) <100)
                {
                     image[j][i]=10*iters%255;
                    //  pixelData = 3*iters%255;
                }
                ++iters;
                ArImg[count*3]=image[j][i];
            }
        }
    }
}

int main(){

    //domain size and mapping
    double domain      = 4.0;               //betweev -2 and +2
    double pixel_width = 2048;                 //2048
    double dx          = (domain/pixel_width)*1;
    int N              = (domain/dx); 
    int Arraay_size    = N*N*3;               //this is for png encoding
    std::complex< double >c_ (-0.8,0.2);      //defining const complex number for Julia
    double pixelData;
    std::vector<double>X(N,1);
    std::vector<double>Y(N,1);
    std::vector<std::vector<double>>image(N,X);     //this we can uncomment 
                                                    //if we want all the values as a 2d matrix 
    unsigned char * ArImg;
    ArImg = new unsigned char [Arraay_size];
    //initialization vectors
    for(int i=0;i<N;++i){
        X[i] = -2+i*dx;
        Y[i] = -2+i*dx;//2-i*dx
    }
    auto START = std::chrono::system_clock::now();
    Julia(X,Y,image,ArImg,c_,N);
    auto END   = std::chrono::system_clock::now();
    std::chrono::duration<double>elepsed_time = END-START;
    // For Encoding png
    const char* filename = "min0802i.png";
    unsigned width       = 2048, height = 2048;
    using lodepng::encode;
    encode(filename, ArImg, width, height,LCT_RGB,8);
    printf("RUN_TIME = %f sec. \n",elepsed_time.count());
    delete[]ArImg;
    return 0;
}*/

#include <iostream>
#include <cmath>
#include <complex>
#include <chrono>
#include "lodepng.h"
#include "Complex.h"


namespace lodepng {
unsigned encode(const std::string& filename,
const unsigned char* image,
unsigned w, unsigned h,
LodePNGColorType colortype,
unsigned bitdepth);
}

// void Julia(double* X, double* Y, double* image, unsigned char* ArImg, std::complex< double >c_,int N){
void Julia(double* X, double* Y, double* image, unsigned char* ArImg, Complex c_,int N){
    for(int i=0;i<N;++i){
        for(int j=0;j<N;++j){
            int iters = 0;                     //iteration count 
            int count = i*N+j;                 //ArImg index
            Complex z(X[j],Y[i]); //definign the complex number
            while(iters<250){
                z = z*z + c_;
                if(z.absolute_val(z) < 100)
                {
                    //  image[j][i]=10*iters%255;
                    image[i*N+j] = 10*iters%255; 
                    //  pixelData = 3*iters%255;
                }
                ++iters;
                ArImg[count*3]=image[i*N+j];
            }
        }
    }
}

int main(){
    
    //domain size and mapping
    double domain      = 4.0;               //betweev -2 and +2
    double pixel_width = 2048;                 //2048
    double dx          = (domain/pixel_width)*1;
    int N              = (domain/dx); 
    int Arraay_size    = N*N*3;               //this is for png encoding
    Complex c_(-0.8,0.2);      //defining const complex number for Julia
    // double pixelData;
    int nBytes         = N*sizeof(double);
    int n2DBytes       = N*N*sizeof(double);

    double* X=nullptr; double* Y=nullptr; double* image=nullptr; unsigned char* ArImg=nullptr;
    X     = (double*)malloc(nBytes);
    Y     = (double*)malloc(nBytes);
    image = (double*)malloc(n2DBytes);
    ArImg = (unsigned char*)malloc(Arraay_size);

    /*std::vector<double>X(N,1);
    std::vector<double>Y(N,1);
    std::vector<std::vector<double>>image(N,X); */    //this we can uncomment 
                                                    //if we want all the values as a 2d matrix 
    //
    // ArImg = new unsigned char [Arraay_size];
    //initialization vectors
    for(int i=0;i<N;++i){
        X[i] = -2+i*dx;
        Y[i] = -2+i*dx;//2-i*dx
    }
    auto START = std::chrono::system_clock::now();
    Julia(X,Y,image,ArImg,c_,N);
    auto END   = std::chrono::system_clock::now();
    std::chrono::duration<double>elepsed_time = END-START;
    // For Encoding png
    const char* filename = "min0802i.png";
    unsigned width       = 2048, height = 2048;
    using lodepng::encode;
    encode(filename, ArImg, width, height,LCT_RGB,8);
    printf("RUN_TIME = %f sec. \n",elepsed_time.count());
    free (X);
    free (Y);
    free (image);
    free (ArImg);
    return 0;
}



