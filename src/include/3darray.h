#ifndef __3DARRAY_H__
#define __3DARRAY_H__

#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <fstream>

template<class T> struct v3d_elem { T u, v, w; };

template<class T>
class s3darray {

protected:
  const size_t _imax;
  const size_t _jmax;
  const size_t _kmax;
  const size_t _size;

  T* _data;

  inline size_t index(const size_t i, const size_t j, const size_t k) const
  {
    return k * _jmax * _imax + j * _imax + i;
  }

public:
  s3darray(const size_t imax, const size_t jmax, const size_t kmax)
    : _imax(imax), _jmax(jmax), _kmax(kmax), _size(imax * jmax * kmax)
  {
    _data = new T[_size];
    memset(_data, 0, _size * sizeof(T));
  }

  virtual ~s3darray()
  {
    delete[] _data;
  }

  void set(const T& val)
  {
    for (size_t i = 0; i < _size; i++) { _data[i] = val; }
  }

  inline T* data() { return _data; }
  inline const T* data() const { return _data; }

  inline T& operator() (const size_t i, const size_t j, const size_t k)
  {
    return _data[index(i,j,k)];
  }

  inline const T& operator() (const size_t i, const size_t j, const size_t k) const
  {
    return _data[index(i,j,k)];
  }

  virtual inline operator       T*()       { return _data; }
  virtual inline operator const T*() const { return _data; }

  inline size_t isize() const { return _imax; }
  inline size_t jsize() const { return _jmax; }
  inline size_t ksize() const { return _kmax; }
  inline size_t size()  const { return _size; }

  inline void copy(const s3darray<T>& s) { std::memcpy(_data, s._data,  _size * sizeof(T)); }
  inline void copy(const s3darray<T>* s) { std::memcpy(_data, s->_data, _size * sizeof(T)); }

  virtual void outsph(const char* fname,
                      const REAL dx, const REAL dy, const REAL dz,
                      const int step, const REAL time)
  {
    using namespace std;

    FILE *fd;
    fd = fopen(fname, "wb");

    if (!fd) return;

    const int imax = _imax;
    const int jmax = _jmax;
    const int kmax = _kmax;

    int _size  = 8;
    int svtype = 1; // scalar
    int dtype  = sizeof(T) / sizeof(float);
    fwrite(&_size,  sizeof(int), 1, fd);
    fwrite(&svtype, sizeof(int), 1, fd);
    fwrite(&dtype,  sizeof(int), 1, fd);
    fwrite(&_size,  sizeof(int), 1, fd);

    _size = dtype * 12;
    fwrite(&_size, sizeof(int), 1, fd);
    fwrite(&imax,  sizeof(int), 1, fd);
    fwrite(&jmax,  sizeof(int), 1, fd);
    fwrite(&kmax,  sizeof(int), 1, fd);
    fwrite(&_size, sizeof(int), 1, fd);

    // origin
    _size = dtype * 12;
    REAL org = 0.0;
    fwrite(&_size, sizeof(int),  1, fd);
    fwrite(&org,   sizeof(REAL), 1, fd);
    fwrite(&org,   sizeof(REAL), 1, fd);
    fwrite(&org,   sizeof(REAL), 1, fd);
    fwrite(&_size, sizeof(int),  1, fd);

    // pitch
    _size = dtype * 12;
    REAL ddx = 1.0;
    REAL ddy = 1.0;
    REAL ddz = 1.0;
    fwrite(&_size, sizeof(int),  1, fd);
    fwrite(&ddx,   sizeof(REAL), 1, fd);
    fwrite(&ddy,   sizeof(REAL), 1, fd);
    fwrite(&ddz,   sizeof(REAL), 1, fd);
    fwrite(&_size, sizeof(int),  1, fd);

    // pitch
    _size = dtype * 8;
    fwrite(&_size, sizeof(int),  1, fd);
    fwrite(&step,  sizeof(int),  1, fd);
    fwrite(&time,  sizeof(REAL), 1, fd);
    fwrite(&_size, sizeof(int),  1, fd);

    // data
    _size = 4;
    fwrite(&_size, sizeof(int), 1, fd);
    fwrite(_data,  sizeof(REAL), _size, fd);
    fwrite(&_size, sizeof(int), 1, fd);

    fclose(fd);
  }

  void dump() const
  {
    for (size_t s = 0; s < _size; s++) {
      std::cout << ((REAL*)_data)[s] << ' ';
    }
    std::cout << std::endl;
  }
};


template<class T>
class v3darray : public s3darray<v3d_elem<T> > {

  typedef s3darray<v3d_elem<T> > base;

public:
  v3darray(const size_t imax, const size_t jmax, const size_t kmax)
    : base(imax,jmax,kmax)
  {
  }

  void set(const T& val)
  {
    for (size_t i = 0; i < base::_size; i++) {
      base::_data[i].u = val;
      base::_data[i].v = val;
      base::_data[i].w = val;
    }
  }

  virtual inline operator       T*()       { return reinterpret_cast<T*>(base::_data); }
  virtual inline operator const T*() const { return reinterpret_cast<const T*>(base::_data); }

  void setu(const T& val)
  {
    for (size_t i = 0; i < base::_size; i++) { base::_data[i].u = val; }
  }

  void setv(const T& val)
  {
    for (size_t i = 0; i < base::_size; i++) { base::_data[i].v = val; }
  }

  void setw(const T& val)
  {
    for (size_t i = 0; i < base::_size; i++) { base::_data[i].w = val; }
  }

  void reduction(const s3darray<T>& x, const s3darray<T>& y, const s3darray<T>& z)
  {
    for (size_t k = 0; k < base::_kmax - 1; k++) {
    for (size_t j = 0; j < base::_jmax - 1; j++) {
    for (size_t i = 0; i < base::_imax - 1; i++) {
      v3d_elem<T>& e = base::_data[base::index(i,j,k)];
      e.u = (x(i,j,k) + x(i+1,j,k)) * 0.5;
      e.v = (y(i,j,k) + y(i,j+1,k)) * 0.5;
      e.w = (z(i,j,k) + z(i,j,k+1)) * 0.5;
    }}}

    // i = imax - 1
    for (size_t k = 0; k < base::_kmax; k++) {
    for (size_t j = 0; j < base::_jmax; j++) {
      size_t i = base::_imax - 1;
      v3d_elem<T>& e = base::_data[base::index(i,j,k)];
      e.u = x(i,j,k);
      e.v = y(i,j,k);
      e.w = z(i,j,k);
    }}
    // j = jmax - 1
    for (size_t k = 0; k < base::_kmax; k++) {
    for (size_t i = 0; i < base::_imax; i++) {
      size_t j = base::_jmax - 1;
      v3d_elem<T>& e = base::_data[base::index(i,j,k)];
      e.u = x(i,j,k);
      e.v = y(i,j,k);
      e.w = z(i,j,k);
    }}
    // k = kmax - 1
    for (size_t j = 0; j < base::_jmax; j++) {
    for (size_t i = 0; i < base::_imax; i++) {
      size_t k = base::_kmax - 1;
      v3d_elem<T>& e = base::_data[base::index(i,j,k)];
      e.u = x(i,j,k);
      e.v = y(i,j,k);
      e.w = z(i,j,k);
    }}
  }

  virtual void outsph(const char* fname,
            const REAL dx, const REAL dy, const REAL dz,
            const int step, const REAL time, const size_t vc = 0)
  {
    const int imax = base::_imax - vc * 2;
    const int jmax = base::_jmax - vc * 2;
    const int kmax = base::_kmax - vc * 2;

    FILE *fd;
    fd = fopen(fname, "wb");

    if (!fd) return;

    int _size  = 8;
    int svtype = 2; // vector
    int dtype  = sizeof(T) / sizeof(float);
    fwrite(&_size,  sizeof(int), 1, fd);
    fwrite(&svtype, sizeof(int), 1, fd);
    fwrite(&dtype,  sizeof(int), 1, fd);
    fwrite(&_size,  sizeof(int), 1, fd);

    _size = dtype * 12;
    fwrite(&_size, sizeof(int), 1, fd);
    fwrite(&imax,  sizeof(int), 1, fd);
    fwrite(&jmax,  sizeof(int), 1, fd);
    fwrite(&kmax,  sizeof(int), 1, fd);
    fwrite(&_size, sizeof(int), 1, fd);

    // origin
    _size = dtype * 12;
    REAL org = 0.0;
    fwrite(&_size, sizeof(int),  1, fd);
    fwrite(&org,   sizeof(REAL), 1, fd);
    fwrite(&org,   sizeof(REAL), 1, fd);
    fwrite(&org,   sizeof(REAL), 1, fd);
    fwrite(&_size, sizeof(int),  1, fd);

    // pitch
    _size = dtype * 12;
    REAL ddx = 1.0;
    REAL ddy = 1.0;
    REAL ddz = 1.0;
    fwrite(&_size, sizeof(int),  1, fd);
    fwrite(&ddx,   sizeof(REAL), 1, fd);
    fwrite(&ddy,   sizeof(REAL), 1, fd);
    fwrite(&ddz,   sizeof(REAL), 1, fd);
    fwrite(&_size, sizeof(int),  1, fd);

    // time
    _size = dtype * 8;
    fwrite(&_size, 1, sizeof(int),  fd);
    fwrite(&step,  1, sizeof(int),  fd);
    fwrite(&time,  1, sizeof(REAL), fd);
    fwrite(&_size, 1, sizeof(int),  fd);

    // data
    _size = imax * jmax * kmax * (4 * dtype) * 3;
    fwrite(&_size,      sizeof(int), 1, fd);
    for (size_t k = 0; k < kmax; k++) {
    for (size_t j = 0; j < jmax; j++) {
    for (size_t i = 0; i < imax; i++) {
      const size_t idx = base::index(i,j,k);
      fwrite(&(base::_data[idx]), sizeof(v3d_elem<T>), 1, fd);
    }}}
    fwrite(&_size,      sizeof(int), 1, fd);

    fclose(fd);
  }

  virtual void outgplot(const char* fname,
            const REAL dx, const REAL dy, const REAL dz,
            const REAL ox, const REAL oy, const REAL oz,
            const REAL sx, const REAL sy, const REAL sz)
  {
    std::ofstream ofs(fname);
    for (size_t k = 0; k < base::_kmax; k++) {
    for (size_t j = 0; j < base::_jmax; j++) {
    for (size_t i = 0; i < base::_imax; i++) {
      ofs << "set arrow from "
        << ox + i * dx << ','
        << oy + j * dy << ','
        << oz + k * dz << " to "
        << ox + i * dx + sx * base::_data[base::index(i,j,k)].u << ','
        << oy + j * dy + sy * base::_data[base::index(i,j,k)].v << ','
        << oz + k * dz + sz * base::_data[base::index(i,j,k)].w << " # "
        << base::_data[base::index(i,j,k)].u << ','
        << base::_data[base::index(i,j,k)].v << ','
        << base::_data[base::index(i,j,k)].w << std::endl;
    }}}
    ofs.close();
  }

  virtual void dump() const
  {
    for (size_t k = 0; k < base::_kmax; k++) {
    for (size_t j = 0; j < base::_jmax; j++) {
    for (size_t i = 0; i < base::_imax; i++) {
      char idx[256];
      sprintf(idx, "%d, %d, %d : ",i,j,k);
      std::cout << idx
        << base::_data[base::index(i,j,k)].u << ' '
        << base::_data[base::index(i,j,k)].v << ' '
        << base::_data[base::index(i,j,k)].w << std::endl;
    }}}
  }

};

#endif // __3DARRAY_H__

