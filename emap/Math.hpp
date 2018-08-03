#include "Common.hpp"

class FFT
{
public:
  cvec in;

  FFT(bool reverse = false, unsigned prepare = 1);

  cvec & execute(cvec const & copyee);
  cvec & execute();
  
  ~FFT();

private:
  cvec out;
  int direction, flags;
  void * plan;

  cmplx * inp;
  int n;
};
