#include "Math.hpp"

#include <fftw3.h>

// TODO: we can support changing the cvec type here
// but need to handle fftw function prefixes changing for different types

FFT::FFT(bool reverse, unsigned prepare)
: direction(reverse ? FFTW_BACKWARD : FFTW_FORWARD),
  flags(prepare == 0 ? FFTW_ESTIMATE : (prepare == 1 ? FFTW_MEASURE : FFTW_PATIENT)),
  plan(nullptr), inp(nullptr), n(0)
{ }


cvec & FFT::execute(cvec const & copyee)
{
  in = copyee;
  return execute();
}

cvec & FFT::execute()
{
  if (inp != &in[0] || n != in.size()) {
    inp = &in[0];
    n = in.size();
    out.resize(n);
    if (plan != nullptr) {
      fftw_destroy_plan(reinterpret_cast<fftw_plan>(plan));
    }
    plan = fftw_plan_dft_1d(n, reinterpret_cast<fftw_complex*>(inp), reinterpret_cast<fftw_complex*>(&out[0]), direction, flags);
  }
  fftw_execute(reinterpret_cast<fftw_plan>(plan));
  return out;
}

FFT::~FFT()
{
  fftw_destroy_plan(reinterpret_cast<fftw_plan>(plan));
}
