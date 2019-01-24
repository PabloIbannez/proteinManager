/*Pablo Ibáñez Freire 2018, pablo.ibannez@uam.es*/

/*
 * Different vector type data structures are defined 
 * as well as different operations such as addition, subtraction ... .
 */


#ifndef VECTOR_DATA_CUH
#define VECTOR_DATA_CUH

#include <iostream>

namespace vectorData{

#ifdef DOUBLE_PRECISION
typedef double real;
#else
typedef float real;
#endif

#ifdef __CUDACC__
#define CUDA_TOKENS __host__ __device__
#endif

struct real2{
    real x;
    real y;
};

struct real3{
    real x;
    real y;
    real z;

    friend std::ostream& operator<<(std::ostream& os, const real3& r);
};

CUDA_TOKENS real3 operator +(const real3 &a, const real3 &b) {
    real3 r;
    r.x = a.x + b.x;
    r.y = a.y + b.y;
    r.z = a.z + b.z;
    return r;
}

CUDA_TOKENS real3 operator +(const real3 &a, const real &b) {
    real3 r;
    r.x = a.x + b;
    r.y = a.y + b;
    r.z = a.z + b;
    return r;
}

CUDA_TOKENS real3 operator -(const real3 &a, const real3 &b) {
    real3 r;
    r.x = a.x - b.x;
    r.y = a.y - b.y;
    r.z = a.z - b.z;
    return r;
}

CUDA_TOKENS real3 operator -(const real3 &a, const real &b) {
    real3 r;
    r.x = a.x - b;
    r.y = a.y - b;
    r.z = a.z - b;
    return r;
}

CUDA_TOKENS real3 operator /(const real3 &a, const real &b) {
    real3 r;
    r.x = a.x / b;
    r.y = a.y / b;
    r.z = a.z / b;
    return r;
}

CUDA_TOKENS real3 operator *(const real3 &a, const real &b) {
    real3 r;
    r.x = a.x * b;
    r.y = a.y * b;
    r.z = a.z * b;
    return r;
}

CUDA_TOKENS void operator +=( real3 &a, const real3 &b) {
    a.x = a.x + b.x;
    a.y = a.y + b.y;
    a.z = a.z + b.z;
}

CUDA_TOKENS void operator /=( real3 &a, const real &b) {
    a.x = a.x / b;
    a.y = a.y / b;
    a.z = a.z / b;
}

CUDA_TOKENS real dot(const real3 &a,const real3 &b){
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

std::ostream& operator<<(std::ostream& os, const real3& r) {
    os << r.x << " " << r.y << " " << r.z;
    return os;
}

struct real4{
    real x;
    real y;
    real z;
    real w;
};

CUDA_TOKENS real3 make_real3(real4 a){
    return {a.x,a.y,a.z};
}

}
#endif
