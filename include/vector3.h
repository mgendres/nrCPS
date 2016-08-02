#include "config.h"
#include "enum.h"
#include <vector>


#ifndef INCLUDED_VECTOR3
#define INCLUDED_VECTOR3
class Vector3
{
  private:

    int x;
    int y;
    int z;

  public:
    explicit Vector3(int, int, int);
    explicit Vector3(int);
    explicit Vector3();
    ~Vector3();

    void SetX(int);
    void SetY(int);
    void SetZ(int);
    void Set(int, int, int);

    int GetX() const;
    int GetY() const;
    int GetZ() const;

    Vector3& operator=(int);
    Vector3& operator+=(const Vector3&);
    Vector3& operator-=(const Vector3&);
    Vector3& operator*=(const Vector3&);
    Vector3& operator*=(int);
    Vector3& operator/=(const Vector3&);
    Vector3& operator/=(int);
    Vector3& operator%=(const Vector3&);
    Vector3& operator%=(int);

    int L1NormSquared();
    int L2NormSquared();

};

// Helper functions
const Vector3 operator+(const Vector3&, const Vector3&);
const Vector3 operator-(const Vector3&, const Vector3&);
const Vector3 operator*(const Vector3&, const Vector3&);
const Vector3 operator*(int, const Vector3&);
const Vector3 operator*(const Vector3&, int);
const Vector3 operator/(const Vector3&, const Vector3&);
const Vector3 operator/(const Vector3&, int);
const Vector3 operator/(int, const Vector3&);
const Vector3 operator%(const Vector3&, const Vector3&);
const Vector3 operator%(const Vector3&, int);
const Vector3 operator%(int, const Vector3&);
bool operator==(const Vector3&, const Vector3&);
const Vector3 abs(const Vector3&);


int Vector3ToIndex(const Vector3&); // This should be obsolete
const Vector3 IndexToVector3(int); // This should be obsolete

int Vector3ToIndex(const Vector3&, const Vector3&);
const Vector3 IndexToVector3(int, const Vector3&);

int KroneckerDelta(const Vector3&);
int KroneckerDelta(const Vector3&, const Vector3&);
int Dot(const Vector3&, const Vector3&);

#endif
