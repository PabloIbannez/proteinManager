#include "dataTypes.hpp"

namespace proteinManager {

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
    
    CUDA_TOKENS real3 operator *(const real &a, const real3 &b) {
        real3 r;
        r.x = a * b.x;
        r.y = a * b.y;
        r.z = a * b.z;
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
    
    CUDA_TOKENS real dot(const real3 &a,const real3 &b) {
        return a.x*b.x + a.y*b.y + a.z*b.z;
    }
    
    CUDA_TOKENS real3 cross(const real3 &a,const real3 &b) {
        return {a.y*b.z - a.z*b.y, (-a.x*b.z + a.z*b.x), a.x*b.y - a.y*b.x};
    }
    
    std::ostream& operator<<(std::ostream& os, const real3& r) {
        os << r.x << " " << r.y << " " << r.z;
        return os;
    }
    
    std::ostream& operator<<(std::ostream& os, const real4& r) {
        os << r.x << " " << r.y << " " << r.z << " " << r.w;
        return os;
    }
}
