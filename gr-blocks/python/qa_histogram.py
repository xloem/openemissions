#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2021 Free Software Foundation, Inc..
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

from gnuradio import gr, gr_unittest
from gnuradio import blocks
import numpy as np
try:
    from openemissions import *
except ImportError:
    import os
    import sys
    dirname, filename = os.path.split(os.path.abspath(__file__))
    sys.path.append(os.path.join(dirname, "bindings"))
    from openemissions import *

def tag(offset, key, value = None):
    if np.issubdtype(type(value), np.float):
        value = gr.pmt.from_double(value)
    elif np.issubdtype(type(value), np.integer):
        value = gr.pmt.from_long(value)
    elif type(value) is str:
        value = gr.pmt.string_to_symbol(value)
    if type(key) is str:
        key = gr.pmt.string_to_symbol(key)
    return gr.tag_utils.python_to_tag((offset, key, value))

class qa_histogram(gr_unittest.TestCase):

    def setUp(self):
        self.tb = gr.top_block()
        self.id_key = 'id'

    def tearDown(self):
        self.tb = None

    def test_instance(self):
        f64_u64 = histogram_f64_f64(-1, 1)
        f32_u64 = histogram_f32_f64(-1, 1)
        s64_u64 = histogram_s64_f64(-1024, 1024)
        s32_u64 = histogram_s32_f64(-1024, 1024)
        s16_u64 = histogram_s16_f64(-1024, 1024)
        s8_u64  = histogram_s8_f64(-127, 127)
        u64_u64 = histogram_u64_f64(0, 1024)
        u8_u64  = histogram_u8_f64(0, 255)

        f64_u64 = histogram_f64_f32(-1, 1)
        f32_u64 = histogram_f32_f32(-1, 1)
        s64_u64 = histogram_s64_f32(-1024, 1024)
        s32_u64 = histogram_s32_f32(-1024, 1024)
        s16_u64 = histogram_s16_f32(-1024, 1024)
        s8_u64  = histogram_s8_f32(-127, 127)
        u64_u64 = histogram_u64_f32(0, 1024)
        u8_u64  = histogram_u8_f32(0, 255)

        f64_u64 = histogram_f64_u64(-1, 1)
        f32_u64 = histogram_f32_u64(-1, 1)
        s64_u64 = histogram_s64_u64(-1024, 1024)
        s32_u64 = histogram_s32_u64(-1024, 1024)
        s16_u64 = histogram_s16_u64(-1024, 1024)
        s8_u64  = histogram_s8_u64(-127, 127)
        u64_u64 = histogram_u64_u64(0, 1024)
        u8_u64  = histogram_u8_u64(0, 255)

    def test_vinlen1(self):
        hist_f32_f32 = histogram_f32_f32(-1, 1, nbuckets=4, prop_tag_keys=[self.id_key])
        src_f32_f32 = blocks.vector_source_f(
            data = [
                0.75, -0.25, 0.75, -0.75, 0.75, -0.25
            ],
            tags = [
                tag(2, self.id_key, 'id1'),
                tag(4, self.id_key, None),
            ]
        )
        sink_f32_f32 = blocks.vector_sink_f(4)

        self.tb.connect(src_f32_f32, hist_f32_f32, sink_f32_f32)

        self.tb.run()

        np.testing.assert_almost_equal([
            0.0, 0.0, 0.0, 1.0,
            0.0, 1/2, 0.0, 1/2,
            0.0, 0.0, 0.0, 1.0,
            1/2, 0.0, 0.0, 1/2,
            0.0, 1/3, 0.0, 2/3,
            0.0, 2/4, 0.0, 2/4
        ], sink_f32_f32.data(), 7)

    def test_vinlen2(self):
        hist_f32_f32 = histogram_f32_f32(-1, 1, nbuckets=4, vinlen=2, prop_tag_keys=[self.id_key])
        src_f32_f32 = blocks.vector_source_f(
            data = [
                0.75, -0.25,
                0.75, -0.75,
                0.75, -0.25
            ],
            tags = [
                tag(1, self.id_key, 'id1'),
                tag(2, self.id_key, None),
            ],
            vlen = 2
        )
        sink_f32_f32 = blocks.vector_sink_f(4)

        self.tb.connect(src_f32_f32, hist_f32_f32, sink_f32_f32)

        self.tb.run()

        np.testing.assert_almost_equal([
            0.0, 1/2, 0.0, 1/2,
            1/2, 0.0, 0.0, 1/2,
            0.0, 2/4, 0.0, 2/4
        ], sink_f32_f32.data(), 7)


if __name__ == '__main__':
    gr_unittest.run(qa_histogram)
