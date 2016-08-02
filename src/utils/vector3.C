#include <stdio.h>
#include <stdarg.h>
using namespace std;
#include "vector3.h"
#include "global_job_parameter.h"
#include "error.h"

Vector3::Vector3(int x0, int y0, int z0)
: x(x0),
  y(y0),
  z(z0)
{
}

Vector3::Vector3(int w)
: x(w),
  y(w),
  z(w)
{
}

Vector3::Vector3()
: x(0),
  y(0),
  z(0)
{
}

Vector3::~Vector3()
{
}

void Vector3::SetX(int x0) { x=x0; }

void Vector3::SetY(int y0) { y=y0; }

void Vector3::SetZ(int z0) { z=z0; }

void Vector3::Set(int x0, int y0, int z0) { x=x0; y=y0; z=z0; }

int Vector3::GetX() const { return x; }

int Vector3::GetY() const { return y; }

int Vector3::GetZ() const { return z; }

Vector3& Vector3::operator=(int rhs) { x = rhs; y = rhs; z = rhs; }

Vector3& Vector3::operator+=(const Vector3& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; }

Vector3& Vector3::operator-=(const Vector3& rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; }

Vector3& Vector3::operator*=(const Vector3& rhs) { x *= rhs.x; y *= rhs.y; z *= rhs.z; }

Vector3& Vector3::operator/=(const Vector3& rhs) { x /= rhs.x; y /= rhs.y; z /= rhs.z; }

Vector3& Vector3::operator%=(const Vector3& rhs) { x %= rhs.x; y %= rhs.y; z %= rhs.z; }

Vector3& Vector3::operator*=(int rhs) { x *= rhs; y *= rhs; z *= rhs; }

Vector3& Vector3::operator/=(int rhs) { x /= rhs; y /= rhs; z /= rhs; }

Vector3& Vector3::operator%=(int rhs) { x %= rhs; y %= rhs; z %= rhs; }

int Vector3::L1NormSquared() { return abs(x)+abs(y)+abs(z); }

int Vector3::L2NormSquared() { return x*x+y*y+z*z; }

const Vector3 operator+(const Vector3& lhs, const Vector3& rhs) {
  Vector3 result(lhs);
  result += rhs;
  return result;
  //return Vector3(lhs) += rhs; // This is better, but seg faults for some reason...
}

const Vector3 operator-(const Vector3& lhs, const Vector3& rhs) {
  Vector3 result(lhs);
  result -= rhs;
  return result;
  //return Vector3(lhs) -= rhs; // This is better, but seg faults for some reason...
}

const Vector3 operator*(const Vector3& lhs, const Vector3& rhs) {
  Vector3 result(lhs);
  result *= rhs;
  return result;
  //return Vector3(lhs) *= rhs; // This is better, but seg faults for some reason...
}

const Vector3 operator/(const Vector3& lhs, const Vector3& rhs) {
  Vector3 result(lhs);
  result /= rhs;
  return result;
  //return Vector3(lhs) /= rhs; // This is better, but seg faults for some reason...
}

const Vector3 operator%(const Vector3& lhs, const Vector3& rhs) {
  Vector3 result(lhs);
  result %= rhs;
  return result;
  //return Vector3(lhs) %= rhs; // This is better, but seg faults for some reason...
}

const Vector3 operator*(int lhs, const Vector3& rhs) {
  Vector3 result(rhs);
  result *= lhs;
  return result;
  //return Vector3(rhs) *= lhs; // This is better, but seg faults for some reason...
}

const Vector3 operator*(const Vector3& lhs, int rhs) {
  Vector3 result(lhs);
  result *= rhs;
  return result;
  //return Vector3(lhs) *= rhs; // This is better, but seg faults for some reason...
}

const Vector3 operator/(const Vector3& lhs, int rhs) {
  Vector3 result(lhs);
  result /= rhs;
  return result;
  //return Vector3(lhs) /= rhs; // This is better, but seg faults for some reason...
}

const Vector3 operator%(const Vector3& lhs, int rhs) {
  Vector3 result(lhs);
  result %= rhs;
  return result;
  //return Vector3(lhs) %= rhs; // This is better, but seg faults for some reason...
}

const Vector3 operator/(int lhs, const Vector3& rhs) {
  return Vector3( lhs/rhs.GetX(), lhs/rhs.GetY(), lhs/rhs.GetZ() );
}

const Vector3 operator%(int lhs, const Vector3& rhs) {
  return Vector3( lhs%rhs.GetX(), lhs%rhs.GetY(), lhs%rhs.GetZ() );
}


bool operator==(const Vector3& lhs, const Vector3& rhs) {
  return ( lhs.GetX() == rhs.GetX() ) && ( lhs.GetY() == rhs.GetY() ) && ( lhs.GetZ() == rhs.GetZ() );
}


const Vector3 abs(const Vector3& vec)
{
  return Vector3( abs(vec.GetX()), abs(vec.GetY()), abs(vec.GetZ()) );
}

int Vector3ToIndex(const Vector3& vec)
{
  int x_sites = GJP.Xsites();
  int y_sites = GJP.Ysites();
  int z_sites = GJP.Zsites();

  int x = vec.GetX();
  if (x<0) { x += x_sites*(1 - x/x_sites); }
  x %= x_sites;

  int y = vec.GetY();
  if (y<0) { y += y_sites*(1 - y/y_sites); }
  y %= y_sites;

  int z = vec.GetZ();
  if (z<0) { z += z_sites*(1 - z/z_sites); }
  z %= z_sites;

  int index(z);
  index += y*z_sites;
  index += x*y_sites*z_sites;
  return index;
}

const Vector3 IndexToVector3(int index)
{
  const char* fname = "Vector3 Vector3ToIndex(int)";
  if ( (index<0) || (index>=GJP.Vol()) ) { ERR.General(fname, "Index out of bounds."); }
  int x_sites = GJP.Xsites();
  int y_sites = GJP.Ysites();
  int z_sites = GJP.Zsites();
  int z = index%z_sites;
  index /= z_sites;
  int y = index%y_sites;
  int x = index/y_sites;
  return Vector3(x,y,z);
}


int Vector3ToIndex(const Vector3& vec, const Vector3& dims)
{
  const char* fname = "Vector3 Vector3ToIndex(const Vector3&, const Vector3 &)";

  int x_sites = dims.GetX();
  int y_sites = dims.GetY();
  int z_sites = dims.GetZ();

  int x = vec.GetX();
  if (x<0) { x += x_sites*(1 - x/x_sites); }
  x %= x_sites;

  int y = vec.GetY();
  if (y<0) { y += y_sites*(1 - y/y_sites); }
  y %= y_sites;

  int z = vec.GetZ();
  if (z<0) { z += z_sites*(1 - z/z_sites); }
  z %= z_sites;

  int index(z);
  index += y*z_sites;
  index += x*y_sites*z_sites;

  return index;
}

const Vector3 IndexToVector3(int index, const Vector3& dims)
{
  const char* fname = "Vector3 IndexToVector3(int, const Vector3 &)";

  int x_sites = dims.GetX();
  int y_sites = dims.GetY();
  int z_sites = dims.GetZ();

  if ( index < 0 ) { ERR.General(fname, "Index out of bounds."); }
  if ( index >= x_sites*y_sites*z_sites ) { ERR.General(fname, "Index out of bounds."); }

  int z = index%z_sites;
  index /= z_sites;
  int y = index%y_sites;
  int x = index/y_sites;

  return Vector3(x,y,z);
}


int KroneckerDelta(const Vector3& lhs)
{
  if (lhs==Vector3(0)) { return 1; } else { return 0; }
}



int KroneckerDelta(const Vector3& lhs, const Vector3& rhs)
{
  if (lhs==rhs) { return 1; } else { return 0; }
}

int Dot(const Vector3& lhs, const Vector3& rhs)
{
  return  lhs.GetX() * rhs.GetX() +  lhs.GetY() * rhs.GetY() + lhs.GetZ() * rhs.GetZ(); 
}
