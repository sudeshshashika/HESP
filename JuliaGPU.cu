#include <iostream>
#include <cmath>
#include <complex>
#include <chrono>
#include <cuda_runtime.h>
#include "lodepng.h"
#include "Complex.h"

namespace lodepng {
unsigned encode(const std::string& filename,
const unsigned char* image,
unsigned w, unsigned h,
LodePNGColorType colortype,
unsigned bitdepth);  
}


__global__ void Julia(double* X, double* Y, double* image, unsigned char* ArImg, Complex c_,int N){
    
    long long idx = blockIdx.x * blockDim.x + threadIdx.x;
    long long idy = blockIdx.y * blockDim.y + threadIdx.y;
    if((idx<N)&&(idy<N)){
            int iters = 0;                     //iteration count                 
            Complex z(X[idy],Y[idx]);          //definign the complex number
            while(iters<250){
                z = z * z + c_;
                if(z.absolute_val(z) <100)
                {
                    image[idx*N+idy] = 10*iters%255; 
                }
                ++iters;
                
                ArImg[(idx*N+idy)*3]=image[idx*N+idy];
            }
    }
    __syncthreads();
}

int main(){

    //domain size and mapping in HOST
    double domain      = 4.0;                   //betweev -2 and +2
    double pixel_width = 2048;                 //2048
    double dx          = (domain/pixel_width)*1;
    int N              = (domain/dx); 
    int Arraay_size    = N*N*3;               //this is for png encoding
    Complex c_ (-0.8,0.2);                      //defining const complex number for Julia
    int nBytes         = N*sizeof(double);
    int n2DBytes       = N*N*sizeof(double);

    double* X=nullptr; double* Y=nullptr; double* image=nullptr; unsigned char* ArImg=nullptr;
    X     = (double*)malloc(nBytes);
    Y     = (double*)malloc(nBytes);
    image = (double*)malloc(n2DBytes);
    ArImg = (unsigned char*)malloc(Arraay_size);

    //initialization vectors
    for(int i=0;i<N;++i){
        X[i] = -2+i*dx;
        Y[i] = -2+i*dx;                 //2-i*dx
    }

    double* d_X;
    double* d_Y;
    double* d_image;
    unsigned char* d_ArImg;

    cudaMalloc(&d_X, nBytes);
    cudaMalloc(&d_Y, nBytes);
    cudaMalloc(&d_image, n2DBytes);
    cudaMalloc(&d_ArImg, Arraay_size);


    cudaMemcpy(d_X, X, nBytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_Y, Y, nBytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_image, image, n2DBytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_ArImg, ArImg, Arraay_size, cudaMemcpyHostToDevice);
    
    int threads_per_block = 32.0;
    int number_of_blocks  = (int)ceil((double)N/threads_per_block);
    auto START = std::chrono::system_clock::now();
    Julia<<< number_of_blocks, threads_per_block >>>(d_X, d_Y, d_image, d_ArImg, c_, N);
    cudaDeviceSynchronize();
    auto END   = std::chrono::system_clock::now();
    std::chrono::duration<double>elepsed_time = END-START;

    cudaMemcpy(ArImg, d_ArImg, Arraay_size, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    // For Encoding png
    const char* filename = "min0802i.png";
    unsigned width       = 2048, height = 2048;
    using lodepng::encode;
    encode(filename, ArImg, width, height,LCT_RGB,8);
    printf("RUN_TIME = %f sec. \n",elepsed_time.count());


    cudaFree(d_X);
    cudaFree(d_Y);
    cudaFree(d_image);
    cudaFree(d_ArImg);

    free (X);
    free (Y);
    free (image);
    free (ArImg);
    return 0;
}