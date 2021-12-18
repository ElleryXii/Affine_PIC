/**
 * The Array3 class.
 * @author Ellery Wang, 2021, adapted from Array2 class by
 * @author Ante Qu, 2017
 * Based on Bridson's simple_flip2d starter code at http://www.cs.ubc.ca/~rbridson/
 */

//TODO: NOT IMPLEMENTED!

#ifndef ARRAY3
#define ARRAY3

#include <cstdio>
#include <cmath>
#include <cstring>

template<class T>
struct Array3{
    int nx, ny, nz;
    int size;
    T *data;

    Array3()
            :nx(0), ny(0), nz(0), size(0), data(0)
    {}

    Array3(int nx_, int ny_, int nz_)
            :nx(0), ny(0), nz(0), size(0), data(0)
    { init(nx_, ny_, nz_); }

    void init(int nx_, int ny_, int nz_)
    {
        delete_memory();
        nx=nx_;
        ny=ny_;
        nz=nz_
        size=nx*ny*nz;
        data=new T[size];
        zero();
    }

    ~Array3()
    { delete_memory(); }

    void delete_memory()
    {
        delete[] data; data=0;
        nx=ny=nz=size=0;
    }


    const T &operator() (int i, int j, int z) const
    { return data[i+nx*j]; }


    T &operator() (int i, int j, int z)
    { return data[i+nx*j]; }

    T bilerp(int i, int j, T fx, T fy)
    { return (1-fx)*((1-fy)*(*this)(i,j)+fy*(*this)(i,j+1))+fx*((1-fy)*(*this)(i+1,j)+fy*(*this)(i+1,j+1)); }

    void copy_to(Array3 &a) const
    { std::memcpy(a.data, data, size*sizeof(T)); }


    T infnorm() const
    {
        T r=0;
        for(int i=0; i<size; ++i)
            if(!(std::fabs(data[i])<=r)) r=std::fabs(data[i]);
        return r;
    }

    void zero()
    { std::memset(data, 0, size*sizeof(T)); }


    double dot(const Array3 &a) const
    {
        double r=0;
        for(int i=0; i<size; ++i)
            r+=data[i]*a.data[i];
        return r;
    }

    void increment(double scale, const Array3 &a)
    { for(int i=0; i<size; ++i) data[i]+=scale*a.data[i]; }

    void scale_and_increment(double scale, const Array3 &a)
    { for(int i=0; i<size; ++i) data[i]=scale*data[i]+a.data[i]; }

    void write_matlab(FILE *fp, const char *variable_name)
    {
        fprintf(fp, "%s=[", variable_name);
        for(int i=0; i<nx; ++i){
            for(int j=0; j<ny; ++j)
                fprintf(fp, "%lg ", (double)(*this)(i,j));
            fprintf(fp, ";\n");
        }
        fprintf(fp,"];\n");
    }
};

typedef Array3<float> Array3f;
typedef Array3<double> Array3d;
typedef Array3<char> Array3c;

template<class T>
struct Array3x4{
    int nx, ny;
    int size;
    T *data;

    Array3x3()
            :nx(0), ny(0), size(0), data(0)
    {}

    Array3x3(int nx_, int ny_)
            :nx(0), ny(0), size(0), data(0)
    { init(nx_, ny_); }

    void init(int nx_, int ny_)
    {
        delete_memory();
        nx=nx_;
        ny=ny_;
        size=4*nx*ny*nz;
        data=new T[size];
        zero();
    }

    ~Array3x3()
    { delete_memory(); }

    void delete_memory()
    {
        delete[] data; data=0;
        nx=ny=size=0;
    }

    const T &operator() (int i, int j, int k) const
    { return data[(i+nx*j)*4+k]; }

    T &operator() (int i, int j, int k)
    { return data[(i+nx*j)*4+k]; }

    void zero()
    { std::memset(data, 0, size*sizeof(T)); }

    T infnorm() const
    {
        T r=0;
        for(int i=0; i<size; ++i)
            if(!(std::fabs(data[i])<=r)) r=std::fabs(data[i]);
        return r;
    }

    void write_matlab(FILE *fp, const char *variable_name)
    {
        fprintf(fp, "%s=reshape([", variable_name);
        for(int i=0; i<size; ++i)
            fprintf(fp, "%lg ", (double)data[i]);
        fprintf(fp,"],[3 %d %d]);\n", nx, ny);
    }
};

typedef Array3x4<float> Array3x4f;
typedef Array3x4<double> Array3x4d;

#endif
