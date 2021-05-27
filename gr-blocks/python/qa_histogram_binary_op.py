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

class qa_histogram_binary_op(gr_unittest.TestCase):

    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_instance(self):
        mul_f32 = histogram_binary_op_f32(0, 2, lambda x, y: x * y)
        mul_f64 = histogram_binary_op_f64(0, 2, lambda x, y: x * y)
        mul_u64 = histogram_binary_op_u64(0, 2, lambda x, y: x * y)

    def test_mul_f64(self):
        src1_f32 = blocks.vector_source_f(
            data = [
                #-3   -2   -1     0    1     2     3     4
                0.00, 0.00, 0.25, 0.0, 0.25, 0.50, 0.00, 0.00
            ],
            vlen = 8
        )
        src2_f32 = blocks.vector_source_f(
            data = [
                #-3   -2     -1    0     1     2    3     4
                0.00, 0.00, 0.66, 0.15, 0.00, 0.19, 0.00, 0.00
            ],
            vlen = 8
        )
        mul_f32 = histogram_binary_op_f32(-3.5, 4.5, lambda x, y: x * y, nbuckets=8)
        sink_f32 = blocks.vector_sink_f(8)

        self.tb.connect(src1_f32, mul_f32, sink_f32)
        self.tb.connect(src2_f32, (mul_f32, 1))

        self.tb.run()

        self.assertAlmostEqual(1.0, np.sum(sink_f32.data()), 6)

        np.testing.assert_almost_equal([
            # [-3,-2)
            0.00, 
            # [-2,-1)
                # -1 * 2
            0.25 * 0.19 +
                # 2 * -1
            0.50 * 0.66 +
            0.00, 
            # [-1, 0)
                # 1 * -1
            0.25 * 0.66 +
            0.00, 
            # [ 0, 1)
                # -1 * 0
            0.25 * 0.15 +
                # 1 * 0
            0.25 * 0.15 +
                # 2 * 0
            0.50 * 0.15 + 
            0.00, 
            # [ 1, 2)
                # -1 * -1
            0.25 * 0.66 + 
            0.00, 
            # [ 2, 3)
                # 1 * 2
            0.25 * 0.19 +
            0.00, 
            # [ 3, 4)
            0.00, 
            # [ 4, 5)
                # 2 * 2
            0.50 * 0.19 +
            0.00,
        ], sink_f32.data(), 7)

if __name__ == '__main__':
    gr_unittest.run(qa_histogram_binary_op)
