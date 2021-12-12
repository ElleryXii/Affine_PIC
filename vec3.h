/**
 * Vec3 class for 3D vector support.
 * @author Ellery Wang, 2021, adapted from Vec2 class by
 * @author Ante Qu, 2017
 * Based on Bridson's simple_flip2d starter code at http://www.cs.ubc.ca/~rbridson/
 */

#ifndef VEC3_H
#define VEC3_H

#include <cassert>
#include <cmath>
#include <iostream>
#include "util.h"

template<class T>
struct Vec3
{
    T v[3];

    Vec3(void)
    {}

    Vec3(T value_for_all)
    { v[0]=v[1]=v[2]=value_for_all; }

    template<class S>
    Vec3(const Vec3<S> source)
    { v[0]=(T)source.v[0]; v[1]=(T)source.v[1];v[2]=(T)source.v[2]; }

    template<class S>
    Vec3(const S source[3])
    { v[0]=(T)source[0]; v[1]=(T)source[1]; v[2]=(T)source[2]; }

    Vec3(T v0, T v1, T v2)
    { v[0]=v0; v[1]=v1; v[2] = v[2]}

    T &operator[](int index)
    {
        assert(0<=index && (unsigned int)index<3);
        return v[index];
    }

    const T &operator[](int index) const
    {
        assert(0<=index && (unsigned int)index<3);
        return v[index];
    }

    Vec3<T> operator+=(const Vec3<T> &w)
    { v[0]+=w.v[0]; v[1]+=w.v[1]; v[2]+=w.v[2]; return *this; }

    Vec3<T> operator-=(const Vec3<T> &w)
    { v[0]-=w.v[0]; v[1]-=w.v[1]; v[2]-=w.v[2]; return *this; }

    template<class S>
    Vec3<T> operator*=(S scalar)
    { v[0]*=scalar; v[1]*=scalar; v[2]*=scalar; return *this; }

    template<class S>
    Vec3<T> operator/=(S scalar)
    { v[0]/=scalar; v[1]/=scalar; v[2]/=scalar; return *this; }
};

template<class T> inline T mag2(const Vec3<T> &a)
{ return a.v[0]*a.v[0] + a.v[1]*a.v[1] * a.v[2]*a.v[2]; }

template<class T> inline T mag(const Vec3<T> &a)
{ return std::sqrt(mag2(a)); }

template<class T> inline T dist2(const Vec3<T> &a, const Vec3<T> &b)
{ return sqr(a.v[0]-b.v[0]) + sqr(a.v[1]-b.v[1]) + sqr(a.v[2]-b.v[2]) ; }

template<class T> inline T dist(const Vec3<T> &a, const Vec3<T> &b)
{ return std::sqrt(dist2(a,b)); }

template<class T> inline bool operator==(const Vec3<T> &a, const Vec3<T> &b)
{ return a.v[0]==b.v[0] && a.v[1]==b.v[1] && a.v[2]==b.v[2]; }

template<class T> inline bool operator!=(const Vec3<T> &a, const Vec3<T> &b)
{ return a.v[0]!=b.v[0] || a.v[1]!=b.v[1] || a.v[2]!=b.v[2]; }

template<class T> inline Vec3<T> operator-(const Vec3<T> &a)
{ return Vec3<T>(-a.v[0], -a.v[1], -a.v[2]); }

template<class T> inline Vec3<T> operator+(const Vec3<T> &a, const Vec3<T> &b, const Vec3<T> &c)
{ return Vec3<T>(a.v[0]+b.v[0], a.v[1]+b.v[1], a.v[2]+b.v[2]); }

template<class T>
inline Vec3<T> operator-(const Vec3<T> &a, const Vec3<T> &b)
{ return Vec3<T>(a.v[0]-b.v[0], a.v[1]-b.v[1], a.v[2]-b.v[3]); }

template<class S, class T>
inline Vec3<T> operator*(const Vec3<T> &a, S scalar)
{ return Vec3<T>(scalar*a.v[0], scalar*a.v[1], scalar*a.v[2]); }

template<class S, class T>
inline Vec3<T> operator*(S scalar, const Vec3<T> &a)
{ return Vec3<T>(scalar*a.v[0], scalar*a.v[1], scalar*a.v[2]); }

template<class S, class T>
inline Vec3<T> operator/(const Vec3<T> &a, S scalar)
{ return Vec3<T>(a.v[0]/scalar, a.v[1]/scalar, a.v[2]/scalar); }

template<class T> inline T dot(const Vec3<T> &a, const Vec3<T> &b)
{ return a.v[0]*b.v[0] + a.v[1]*b.v[1] + a.v[2]*b.v[2]; }

//TODO
template<class T> inline T cross(const Vec3<T> &a, const Vec3<T> &b)
{ return a.v[0]*b.v[1]-a.v[1]*b.v[0]; }

//TODO
template<class T> inline Vec3<T> perp(const Vec3<T> &a)
{ return Vec3<T>(-a.v[1], a.v[0]); }

template<class T> inline void normalize(Vec3<T> &a)
{ a/=mag(a); }

template<class T> inline Vec3<T> normalized(const Vec3<T> &a)
{ return a/mag(a); }

template<class T>
inline std::ostream &operator<<(std::ostream &out, const Vec3<T> &a)
{ return out<<a.v[0]<<' '<<a.v[1] << ' '<<a.v[2]; }

template<class T>
inline std::istream &operator>>(std::istream &in, Vec3<T> &a)
{ return in>>a.v[0]>>a.v[1]>>a.v[2]; }

// common types of vectors ===================================================
typedef Vec3<float> Vec3f;
typedef Vec3<double> Vec3d;
typedef Vec3<int> Vec3i;

// type-specific operations ==================================================

/* may raise problems if hash() isn't defined for ints yet
// presupposes a hash() function defined for ints
inline unsigned int hash(const Vec3i &a)
{ return hash(a.v[0]) ^ a.v[1]; }
*/

template<class T> inline Vec3i round(const Vec3<T> &a)
{ return Vec3i(lround(a.v[0]), lround(a.v[1]),  lround(a.v[2])); }

#endif
