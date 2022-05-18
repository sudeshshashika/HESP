#pragma once

#include <cuda_runtime.h>
#include <iostream>

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif 


class Complex{
    CUDA_CALLABLE_MEMBER friend std::ostream &operator<<(std::ostream& os, Complex& obj){
        os<<"("<<obj.r<<","<<obj.i<<")";
        return os;
    }
    CUDA_CALLABLE_MEMBER friend Complex operator*(const Complex& a, const Complex& b){
        double t_r = a.r * b.r - a.i * b.i;
        double t_i = a.r * b.i + b.r * a.i;
        Complex temp(t_r,t_i);
        return temp;     
    }
public:
    double r;
    double i;
    CUDA_CALLABLE_MEMBER Complex(){r=0.0;i=0.0;}
    CUDA_CALLABLE_MEMBER Complex(double x, double y):r{x},i{y}{}
    CUDA_CALLABLE_MEMBER Complex(const Complex& obj):r{obj.r},i{obj.i}{}
    CUDA_CALLABLE_MEMBER Complex &operator=(const Complex& obj){
        if(this==&obj)
            return *this;
        r = obj.r;
        i = obj.i;
        return *this;
    }
    CUDA_CALLABLE_MEMBER Complex operator+(const Complex& b){
        
        double t_r = r+b.r;
        double t_i = i+b.i;
        Complex temp(t_r,t_i);
        return temp;
    }
    CUDA_CALLABLE_MEMBER Complex operator-(const Complex& b){

        double t_r = r-b.r;
        double t_i = i-b.i;
        Complex temp(t_r,t_i);
        return temp;
    }
    CUDA_CALLABLE_MEMBER double absolute_val(const Complex& a){return a.r*a.r + a.i*a.i;}
    CUDA_CALLABLE_MEMBER ~Complex(){}
};