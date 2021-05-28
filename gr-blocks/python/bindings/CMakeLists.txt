# Copyright 2020 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

########################################################################
# Check if there is C++ code at all
########################################################################
if(NOT openemissions_sources)
    MESSAGE(STATUS "No C++ sources... skipping python bindings")
    return()
endif(NOT openemissions_sources)

########################################################################
# Check for pygccxml
########################################################################
GR_PYTHON_CHECK_MODULE_RAW(
    "pygccxml"
    "import pygccxml"
    PYGCCXML_FOUND
    )

include(GrPybind)

########################################################################
# Python Bindings
########################################################################

list(APPEND openemissions_python_files
    accumulate_python.cc
    pad_crop_python.cc
    pigpio_sink_python.cc
    histogram_python.cc
    tagged_stream_histogram_python.cc
    histogram_solve_python.cc
    python_bindings.cc
)

GR_PYBIND_MAKE_OOT(openemissions 
   ../..
   gr::openemissions
   "${openemissions_python_files}")

install(TARGETS openemissions_python DESTINATION ${GR_PYTHON_DIR}/openemissions COMPONENT pythonapi)
