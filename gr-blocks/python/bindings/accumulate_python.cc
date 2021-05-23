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
/* BINDTOOL_HEADER_FILE(accumulate.h)                                              */
/* BINDTOOL_HEADER_FILE_HASH(cbb329ca45bd60803c9b1e29342bd832)                     */
/***********************************************************************************/

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <openemissions/accumulate.h>
// pydoc.h is automatically generated in the build directory
#include <accumulate_pydoc.h>

template <typename T>
void bind_accumulate_template(py::module& m, const char* classname)
{

    using accumulate    = gr::openemissions::accumulate<T>;


    py::class_<accumulate, gr::tagged_stream_block, gr::block, gr::basic_block,
        std::shared_ptr<accumulate>>(m, classname)

        .def(py::init(&accumulate::make),
           py::arg("average") = true,
           py::arg("len_tag_key") = "packet_len",
           py::arg("reset_tag_key") = "reset_sum",
           D(accumulate,make)
        ) ;
}

void bind_accumulate(py::module& m)
{
    bind_accumulate_template<std::complex<double>>(m, "accumulate_c64");
    bind_accumulate_template<std::complex<float>>(m, "accumulate_c32");
    bind_accumulate_template<double>(m, "accumulate_f64");
    bind_accumulate_template<float>(m, "accumulate_f32");
    bind_accumulate_template<int64_t>(m, "accumulate_s64");
    bind_accumulate_template<int>(m, "accumulate_int");
}