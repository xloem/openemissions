#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2021 Free Software Foundation, Inc..
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

from gnuradio import gr, gr_unittest
# from gnuradio import blocks
try:
    from openemissions import histogram_binary_op
except ImportError:
    import os
    import sys
    dirname, filename = os.path.split(os.path.abspath(__file__))
    sys.path.append(os.path.join(dirname, "bindings"))
    from openemissions import histogram_binary_op

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
        src_f64 = blocks.vector_source_f(
            data = [
                
            ]
        )
        self.assertEqual("implemented qa", "did not implement qa yet")
        mul_f64 = histogram_binary_op_f64(0, 2, lambda x, y: x * y)
        # set up fg
        self.tb.run()
        # check data


if __name__ == '__main__':
    gr_unittest.run(qa_histogram_binary_op)
