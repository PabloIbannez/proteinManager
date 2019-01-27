#ifndef DATA_TYPES_HPP
#define DATA_TYPES_HPP

#include <iostream>

namespace proteinManager {

    enum DATA_FORMAT {PDB,PDBQ,PQR,PDRS,SP,SPQ,GRO};
    
    //#define DOUBLE_PRECISION
    
    #ifdef DOUBLE_PRECISION
    typedef double real;
    #else
    typedef float real;
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
    
    real3 operator +(const real3 &a, const real3 &b);
    real3 operator +(const real3 &a, const real &b);
    real3 operator -(const real3 &a, const real3 &b);
    real3 operator -(const real3 &a, const real  &b);
    real3 operator /(const real3 &a, const real &b);
    real3 operator *(const real3 &a, const real &b);
    void operator +=( real3 &a, const real3 &b);
    void operator /=( real3 &a, const real &b);
    
    real dot(const real3 &a,const real3 &b);
    
    std::ostream& operator<<(std::ostream& os, const real3& r);
    
    struct real4 {
        real x;
        real y;
        real z;
        real w;
    };

}
#endif