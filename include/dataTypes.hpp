#ifndef DATA_TYPES_HPP
#define DATA_TYPES_HPP
/* Pablo Ibanez Freire, pablo.ibannez@uam.es*/

/* Here the compilation precission is set. The in/out file formats are listed and all
 * the types used by proteinManager are defined as well as their operators.
 */

#include <iostream>

namespace proteinManager {
    
    //Uncomment the next line and recompile to use the double precision.
    //#define DOUBLE_PRECISION
    
    #ifdef DOUBLE_PRECISION
    typedef double real;
    #else
    typedef float real;
    #endif
    
    //in/out admitted file formats
    enum DATA_FORMAT {PDB,PDBQ,PQR,PDRS,SP,SPQ,GRO,XYZ};
    
    //Flag to enable data CUDA compatibility if nvcc is used as compiler
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
    CUDA_TOKENS real3 operator *(const real  &a, const real3 &b);
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
