/* -*- c++ -*- */

#define OPENEMISSIONS_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "openemissions_swig_doc.i"

%{
#include "openemissions/pigpio_sink.h"
%}


%include "openemissions/pigpio_sink.h"
GR_SWIG_BLOCK_MAGIC2(openemissions, pigpio_sink);
