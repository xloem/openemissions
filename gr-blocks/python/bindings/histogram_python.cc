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
/* BINDTOOL_HEADER_FILE(histogram.h)                                               */
/* BINDTOOL_HEADER_FILE_HASH(b0567cc375f7e5e6c3e41189247ee9af)                     */
/***********************************************************************************/

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <openemissions/histogram.h>
// pydoc.h is automatically generated in the build directory
#include <histogram_pydoc.h>

template <typename input_type, typename freq_type>
void bind_histogram_template(py::module& m, const char* classname)
{

    using histogram    = gr::openemissions::histogram<input_type, freq_type>;


    py::class_<histogram, gr::sync_block, gr::block, gr::basic_block,
        std::shared_ptr<histogram>>(m, classname)

        .def(py::init(&histogram::make),
           py::arg("min"),
           py::arg("max"),
           py::arg("nbuckets") = 1024,
           py::arg("vinlen") = 1,
           py::arg("prop_tag_keys") = std::vector<std::string>{},
           py::arg("filename") = "",
           D(histogram,make)
        ) ;

}

void bind_histogram(py::module& m)
{
    bind_histogram_template<double,   double>(m, "histogram_f64_f64");
    bind_histogram_template<float,    double>(m, "histogram_f32_f64");
    bind_histogram_template<int64_t,  double>(m, "histogram_s64_f64");
    bind_histogram_template<int32_t,  double>(m, "histogram_s32_f64");
    bind_histogram_template<int16_t,  double>(m, "histogram_s16_f64");
    bind_histogram_template<int8_t,   double>(m, "histogram_s8_f64");
    bind_histogram_template<uint64_t, double>(m, "histogram_u64_f64");
    bind_histogram_template<uint8_t,  double>(m, "histogram_u8_f64");

    bind_histogram_template<double,   float>(m, "histogram_f64_f32");
    bind_histogram_template<float,    float>(m, "histogram_f32_f32");
    bind_histogram_template<int64_t,  float>(m, "histogram_s64_f32");
    bind_histogram_template<int32_t,  float>(m, "histogram_s32_f32");
    bind_histogram_template<int16_t,  float>(m, "histogram_s16_f32");
    bind_histogram_template<int8_t,   float>(m, "histogram_s8_f32");
    bind_histogram_template<uint64_t, float>(m, "histogram_u64_f32");
    bind_histogram_template<uint8_t,  float>(m, "histogram_u8_f32");

    bind_histogram_template<double,   uint64_t>(m, "histogram_f64_u64");
    bind_histogram_template<float,    uint64_t>(m, "histogram_f32_u64");
    bind_histogram_template<int64_t,  uint64_t>(m, "histogram_s64_u64");
    bind_histogram_template<int32_t,  uint64_t>(m, "histogram_s32_u64");
    bind_histogram_template<int16_t,  uint64_t>(m, "histogram_s16_u64");
    bind_histogram_template<int8_t,   uint64_t>(m, "histogram_s8_u64");
    bind_histogram_template<uint64_t, uint64_t>(m, "histogram_u64_u64");
    bind_histogram_template<uint8_t,  uint64_t>(m, "histogram_u8_u64");
}
