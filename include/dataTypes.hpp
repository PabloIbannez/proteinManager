#ifndef DATA_TYPES_HPP
#define DATA_TYPES_HPP

#include <iostream>

namespace proteinManager {

    enum DATA_FORMAT {PDB,PDBQ,PQR,PDRS,SP,SPQ,GRO,XYZ};
    
    //#define DOUBLE_PRECISION
    
    #ifdef DOUBLE_PRECISION
    typedef double real;
    #else
    typedef float real;
    #endif
    
    #ifdef __CUDACC__
        #define CUDA_TOKENS __host__ __device__
    #else
        #define CUDA_TOKENS
    #endif
    
    struct int3 {
        int x;
        int y;
        int z;
    };
    
    struct real2 {
        real x;
        real y;
    };
    
    struct real3 {
        real x;
        real y;
        real z;
    
        friend std::ostream& operator<<(std::ostream& os, const real3& r);
    };
    
    CUDA_TOKENS real3 operator +(const real3 &a, const real3 &b);
    CUDA_TOKENS real3 operator +(const real3 &a, const real &b);
    CUDA_TOKENS real3 operator -(const real3 &a, const real3 &b);
    CUDA_TOKENS real3 operator -(const real3 &a, const real  &b);
    CUDA_TOKENS real3 operator /(const real3 &a, const real &b);
    CUDA_TOKENS real3 operator *(const real3 &a, const real &b);
    CUDA_TOKENS void operator +=( real3 &a, const real3 &b);
    CUDA_TOKENS void operator /=( real3 &a, const real &b);
    
    CUDA_TOKENS real dot(const real3 &a,const real3 &b);
    
    std::ostream& operator<<(std::ostream& os, const real3& r);
    
    struct real4 {
        real x;
        real y;
        real z;
        real w;
    };
    
    std::ostream& operator<<(std::ostream& os, const real4& r);

}
#endif
