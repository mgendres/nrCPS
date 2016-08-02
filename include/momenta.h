#include <vector>
#include "enum.h"
#include "arg.h"

#ifndef INCLUDED_MOMENTA
#define INCLUDED_MOMENTA

class Momenta
{
  private:
    Float* omega;
    vector<int>* momentum;
    vector<int> shell_count;
    int* parity_list;
    int momenta_count;

    Momenta& operator=(const Momenta&);
    Momenta(const Momenta&);

  public:
    explicit Momenta(const MomentaArg &);
    ~Momenta();
    int MomentaCount();
    vector<int> GetMomentum(int);
    int OppositeParityIndex(int);
    int ShellCount(int);
    int AccShellCount(int);
    int NumShells();
};

#endif
