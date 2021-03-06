/*
 * Copyright 2021 Free Software Foundation, Inc.
 *
 * This file is part of GNU Radio
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 */

/***********************************************************************************/
/* This file is automatically generated using bindtool and can be manually edited  */
/* The following lines can be configured to regenerate this file during cmake      */
/* If manual edits are made, the following tags should be modified accordingly.    */
/* BINDTOOL_GEN_AUTOMATIC(0)                                                       */
/* BINDTOOL_USE_PYGCCXML(0)                                                        */
/* BINDTOOL_HEADER_FILE(pigpio_sink.h)                                             */
/* BINDTOOL_HEADER_FILE_HASH(14d89953a737c96b5e78f0c1a1c7c55a)                     */
/***********************************************************************************/

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <openemissions/pigpio_sink.h>
// pydoc.h is automatically generated in the build directory
#include <pigpio_sink_pydoc.h>

template <typename T>
void bind_pigpio_sink_template(py::module& m, const char* classname)
{

    using pigpio_sink    = gr::openemissions::pigpio_sink<T>;


    py::class_<pigpio_sink,
               gr::sync_block,
               gr::block,
               gr::basic_block,
               std::shared_ptr<pigpio_sink>>(m, classname)

        .def(py::init(&pigpio_sink::make),
           py::arg("samp_rate"),
           py::arg("pin") =  4,
           py::arg("level") =  0,
           py::arg("address") =  "127.0.0.1:8888",
           py::arg("wave_buffer_percent") =  50,
           py::arg("hardware_clock_frequency") =  30000000,
           py::arg("pad_milliamps") = 0,
           D(pigpio_sink,make)
        ) ;
}

void bind_pigpio_sink(py::module& m)
{
    bind_pigpio_sink_template<float>(m, "pigpio_sink_float");
    bind_pigpio_sink_template<int>(m, "pigpio_sink_int");
    bind_pigpio_sink_template<short>(m, "pigpio_sink_short");
    bind_pigpio_sink_template<signed char>(m, "pigpio_sink_byte");
    bind_pigpio_sink_template<bool>(m, "pigpio_sink_bit");
}
