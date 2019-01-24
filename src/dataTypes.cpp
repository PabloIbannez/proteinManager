#include "dataTypes.hpp"

namespace proteinManager {

    real3 operator +(const real3 &a, const real3 &b) {
        real3 r;
        r.x = a.x + b.x;
        r.y = a.y + b.y;
        r.z = a.z + b.z;
        return r;
    }
    
    real3 operator +(const real3 &a, const real &b) {
        real3 r;
        r.x = a.x + b;
        r.y = a.y + b;
        r.z = a.z + b;
        return r;
    }
    
    real3 operator -(const real3 &a, const real3 &b) {
        real3 r;
        r.x = a.x - b.x;
        r.y = a.y - b.y;
        r.z = a.z - b.z;
        return r;
    }
    
    real3 operator -(const real3 &a, const real &b) {
        real3 r;
        r.x = a.x - b;
        r.y = a.y - b;
        r.z = a.z - b;
        return r;
    }
    
    real3 operator /(const real3 &a, const real &b) {
        real3 r;
        r.x = a.x / b;
        r.y = a.y / b;
        r.z = a.z / b;
        return r;
    }
    
    real3 operator *(const real3 &a, const real &b) {
        real3 r;
        r.x = a.x * b;
        r.y = a.y * b;
        r.z = a.z * b;
        return r;
    }
    
    void operator +=( real3 &a, const real3 &b) {
        a.x = a.x + b.x;
        a.y = a.y + b.y;
        a.z = a.z + b.z;
    }
    
    void operator /=( real3 &a, const real &b) {
        a.x = a.x / b;
        a.y = a.y / b;
        a.z = a.z / b;
    }
    
    real dot(const real3 &a,const real3 &b) {
        return a.x*b.x + a.y*b.y + a.z*b.z;
    }
    
    std::ostream& operator<<(std::ostream& os, const real3& r) {
        os << r.x << " " << r.y << " " << r.z;
        return os;
    }
}
