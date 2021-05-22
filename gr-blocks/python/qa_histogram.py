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

class qa_histogram(gr_unittest.TestCase):

    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_instance(self):
        # FIXME: Test will fail until you pass sensible arguments to the constructor
        instance = histogram_f32_f32(-1, 1)

    def test_001_descriptive_test_name(self):
        f32_f32 = histogram_f32_f32(-1, 1, 4)
        src = blocks.vector_source_f(
            data = [
                0.75, -0.25, 0.75, -0.75, 0.75, -0.25
            ]
        )
        sink_f32_f32 = blocks.vector_sink_f(4)

        self.tb.connect(src, f32_f32, sink_f32_f32)

        self.tb.run()

        np.testing.assert_almost_equal([
            0.0, 0.0, 0.0, 1.0,
            0.0, 1/2, 0.0, 1/2,
            0.0, 1/3, 0.0, 2/3,
            1/4, 1/4, 0.0, 2/4,
            1/5, 1/5, 0.0, 3/5,
            1/6, 2/6, 0.0, 3/6
        ], sink_f32_f32.data(), 7)


if __name__ == '__main__':
    gr_unittest.run(qa_histogram)
