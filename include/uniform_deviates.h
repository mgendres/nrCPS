#include "arg.h"
#include "enum.h"

#ifndef INCLUDED_UNIFORM_DEVIATES
#define INCLUDED_UNIFORM_DEVIATES
class UniformDeviate
// Abstract random number generator base class
{

  private:
    UniformDeviate& operator=(const UniformDeviate&);
    UniformDeviate(const UniformDeviate&);

  public:
    explicit UniformDeviate() {};
    virtual ~UniformDeviate() {};
    virtual void Init(long)=0; // Initialize the generator
    virtual Float Run()=0; // Generate a random number
    virtual int StateSize()=0; // Determine the size of the generator state
    virtual void GetState(long*)=0; // Save the state in integer array
    virtual void SetState(long*)=0; // Set the state from integer array
    virtual RandomType TypeQ()=0; // Returns the generator type
};


class Ran0 : public UniformDeviate
// Random number generator class derived from Random
{
  private:
    long idum;

    Ran0& operator=(const Ran0&);
    Ran0(const Ran0&);

  public:
    explicit Ran0();
    ~Ran0();
    void Init(long);
    Float Run();
    int StateSize();
    void GetState(long*);
    void SetState(long*);
    RandomType TypeQ();
};

class Ran1 : public UniformDeviate
// Random number generator class derived from Random
{
  private:
    long idum;
    long iy;
    long* iv;

    Ran1& operator=(const Ran1&);
    Ran1(const Ran1&);

  public:
    explicit Ran1();
    ~Ran1();
    void Init(long);
    Float Run();
    int StateSize();
    void GetState(long*);
    void SetState(long*);
    RandomType TypeQ();
};



class Ran2 : public UniformDeviate
// Random number generator class derived from Random
{
  private:
    long idum;
    long idum2;
    long iy;
    long* iv;

    Ran2& operator=(const Ran2&);
    Ran2(const Ran2&);

  public:
    explicit Ran2();
    ~Ran2();
    void Init(long);
    Float Run();
    int StateSize();
    void GetState(long*);
    void SetState(long*);
    RandomType TypeQ();
};

class Ran3 : public UniformDeviate
// Random number generator class derived from Random
{
  private:
    int inext;
    int inextp;
    long* ma;

    Ran3& operator=(const Ran3&);
    Ran3(const Ran3&);

  public:
    explicit Ran3();
    ~Ran3();
    void Init(long);
    Float Run();
    int StateSize();
    void GetState(long*);
    void SetState(long*);
    RandomType TypeQ();
};

//---- Ranlxs random generator declarations
typedef struct
{
   int c1,c2,c3,c4;
} vec_t;

typedef struct
{
   vec_t c1,c2;
} dble_vec_t;

class Ranlxs : public UniformDeviate
// Random number generator class derived from Random
{
  private:
//    int init;
    int pr;
    int prm;
    int ir;
    int jr;
    int is;
    int is_old;
    int next[96];
    float one_bit;
    vec_t carry;
    union { dble_vec_t vec[12]; int num[96]; } x;
 
    void Update(void);
    void DefineConstants(void);

    void LuxLevel(int);
    int lux_level;

    Ranlxs& operator=(const Ranlxs&);
    Ranlxs(const Ranlxs&);

  public:
    explicit Ranlxs(int);
    ~Ranlxs();
    void Init(long);
    Float Run();
    int StateSize();
    void GetState(long*);
    void SetState(long*);
    RandomType TypeQ();
};

class Ranlxd : public UniformDeviate
// Random number generator class derived from Random
{
  private:
//    int init;
    int pr;
    int prm;
    int ir;
    int jr;
    int is;
    int is_old;
    int next[96];
    float one_bit;
    vec_t carry;
    union { dble_vec_t vec[12]; int num[96]; } x;
 
    void Update(void);
    void DefineConstants(void);

    void LuxLevel(int);
    int lux_level;

    Ranlxd& operator=(const Ranlxd&);
    Ranlxd(const Ranlxd&);

  public:
    explicit Ranlxd(int);
    ~Ranlxd();
    void Init(long);
    Float Run();
    int StateSize();
    void GetState(long*);
    void SetState(long*);
    RandomType TypeQ();
};

#endif
