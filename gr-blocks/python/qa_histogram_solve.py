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

class qa_histogram_solve(gr_unittest.TestCase):

    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_instance(self):
        f64_1 = histogram_solve_f64_1(-1, 1, lambda a: a, nbuckets=4)
        f32_1 = histogram_solve_f32_1(-1, 1, lambda a: a, nbuckets=4)
        u64_1 = histogram_solve_u64_1(-1, 1, lambda a: a, nbuckets=4)
        f64_2 = histogram_solve_f64_2(-1, 1, lambda a, b: a, nbuckets=4)
        f32_2 = histogram_solve_f32_2(-1, 1, lambda a, b: a, nbuckets=4)
        u64_2 = histogram_solve_u64_2(-1, 1, lambda a, b: a, nbuckets=4)
        f64_3 = histogram_solve_f64_3(-1, 1, lambda a, b, c: a, nbuckets=4)
        f32_3 = histogram_solve_f32_3(-1, 1, lambda a, b, c: a, nbuckets=4)
        u64_3 = histogram_solve_u64_3(-1, 1, lambda a, b, c: a, nbuckets=4)

    def test_scalar_reverse_add_f32(self):
        src_sum_f32 = blocks.vector_source_f(
            data = [
                #-3   -2   -1     0    1     2     3     4
                0,     0,   0,    1,   2,    1,    2,    3
            ],
            vlen = 8
        )
        src_addend_f32 = blocks.vector_source_f(
            data = [
                #-3   -2   -1     0    1     2     3     4
                1,     2,   1,    2,   3,    0,    0,    0
            ],
            vlen = 8
        )
        reverse_add_f32 = histogram_solve_f32_2(-3.5, 4.5, lambda x, y: x + y, output_idx = 1, nbuckets = 8)
        sink_addend_f32 = blocks.vector_sink_f(8)

        self.tb.connect(src_sum_f32, reverse_add_f32, sink_addend_f32)
        self.tb.connect(src_addend_f32, (reverse_add_f32, 1))

        self.tb.run()

        print(sink_addend_f32.data())

    def test_scalar_reverse_add_f32_random_precise(self):
        from random import randint
        sum_data = [1, randint(2, 40), randint(2, 40), randint(2, 40), randint(2, 40), 1, 0, 0]
        addend_data = [0, 0, 1, randint(2, 40), randint(2, 40), randint(2, 40), randint(2, 40), 1]
        src_sum_f32 = blocks.vector_source_f(
            data = sum_data,
            vlen = 8
        )
        src_addend_f32 = blocks.vector_source_f(
            data = addend_data,
            vlen = 8
        )
        reverse_add_f32 = histogram_solve_f32_2(-3.5, 4.5, lambda x, y: x + y, output_idx = 1, nbuckets = 8)
        sink_addend_f32 = blocks.vector_sink_f(8)

        self.tb.connect(src_sum_f32, reverse_add_f32, sink_addend_f32)
        self.tb.connect(src_addend_f32, (reverse_add_f32, 1))

        self.tb.run()

        print(addend_data, '+ x =', sum_data)
        print(sink_addend_f32.data())

    def test_scalar_reverse_add_f32_sum_random(self):
        from random import choices
        sum_dist = [0, 1/4, 1/4, 1/4, 1/4, 0, 0, 0]
        addend_dist = [0, 0, *sum_dist[:-2]]
        sum_data = [0] * len(sum_dist)
        addend_data = [0] * len(addend_dist)
        for item in choices(range(len(sum_dist)), sum_dist, k=100000):
            sum_data[item] += 1
        src_sum_f32 = blocks.vector_source_f(
            data = sum_data,
            #data = [0, 20/4, 20/4, 20/4, 20/4, 0, 0, 0],
            #data = [0, 1/4, 1/4, 1/4, 1/4, 0, 0, 0],
            vlen = 8
        )
        #break1 = randint(0, 40)
        #break2 = randint(0, 40) + break1
        #break3 = randint(0, 40) + break2
        #total = randint(0, 40) + break3
        for item in choices(range(len(sum_dist)), sum_dist, k=100000):
            addend_data[item] += 1
        src_addend_f32 = blocks.vector_source_f(
            #data = [0, 0, 0, randint(0, 40), randint(0, 40), randint(0, 40), randint(0, 40), 0],
            #data = [0, 0, 0, break1 / total, (break2 - break1) / total, (break3 - break2) / total, (total - break3) / total, 0],
            #data = [0, 0, 0, 300, 501, 150, 503, 0],
            data = addend_data,
            vlen = 8
        )
        reverse_add_f32 = histogram_solve_f32_2(-3.5, 4.5, lambda x, y: x + y, output_idx = 1, nbuckets = 8)
        sink_addend_f32 = blocks.vector_sink_f(8)

        self.tb.connect(src_sum_f32, reverse_add_f32, sink_addend_f32)
        self.tb.connect(src_addend_f32, (reverse_add_f32, 1))

        self.tb.run()

        print(sink_addend_f32.data())

    def test_solve_identity_mul_f32(self):
        from random import choices, randint
        dist = [randint(0, 100) for x in range(8)]
        #data1 = 
        #for item in choices(range(len(sum_dist)), sum_dist, k=100000):
        sum_dist = [0, 1/4, 1/4, 1/4, 1/4, 0, 0, 0]
        addend_dist = [0, 0, *sum_dist[:-2]]
        sum_data = [0] * len(sum_dist)
        addend_data = [0] * len(addend_dist)
        for item in choices(range(len(sum_dist)), sum_dist, k=100000):
            sum_data[item] += 1
        src_sum_f32 = blocks.vector_source_f(
            data = sum_data,
            #data = [0, 20/4, 20/4, 20/4, 20/4, 0, 0, 0],
            #data = [0, 1/4, 1/4, 1/4, 1/4, 0, 0, 0],
            vlen = 8
        )
        #break1 = randint(0, 40)
        #break2 = randint(0, 40) + break1
        #break3 = randint(0, 40) + break2
        #total = randint(0, 40) + break3
        for item in choices(range(len(sum_dist)), sum_dist, k=100000):
            addend_data[item] += 1
        src_addend_f32 = blocks.vector_source_f(
            #data = [0, 0, 0, randint(0, 40), randint(0, 40), randint(0, 40), randint(0, 40), 0],
            #data = [0, 0, 0, break1 / total, (break2 - break1) / total, (break3 - break2) / total, (total - break3) / total, 0],
            #data = [0, 0, 0, 300, 501, 150, 503, 0],
            data = addend_data,
            vlen = 8
        )
        reverse_add_f32 = histogram_solve_f32_2(-3.5, 4.5, lambda x, y: x + y, output_idx = 1, nbuckets = 8)
        sink_addend_f32 = blocks.vector_sink_f(8)

        self.tb.connect(src_sum_f32, reverse_add_f32, sink_addend_f32)
        self.tb.connect(src_addend_f32, (reverse_add_f32, 1))

        self.tb.run()

        print(sink_addend_f32.data())

    def test_2hist_forward_mul_f32(self):
        src1_f32 = blocks.vector_source_f(
            data = [
                #-3   -2   -1     0    1     2     3     4
                0,     0,   1,    0,   1,    2,    0,    0
            ],
            vlen = 8
        )
        src2_f32 = blocks.vector_source_f(
            data = [
                #-3   -2   -1     0    1     2     3     4
                0,     0,   4,    1,   0,    1,    0,    0 
            ],
            vlen = 8
        )
        mul_f32 = histogram_solve_f32_2(-3.5, 4.5, lambda x, y: x * y, nbuckets=8)
        sink_f32 = blocks.vector_sink_f(8)

        self.tb.connect(src1_f32, mul_f32, sink_f32)
        self.tb.connect(src2_f32, (mul_f32, 1))

        self.tb.run()

        np.testing.assert_almost_equal([
            # [-3,-2)
            0.00, 
            # [-2,-1)
                # -1 * 2
            1 * 1 +
                # 2 * -1
            2 * 4 +
            0.00, 
            # [-1, 0)
                # 1 * -1
            1 * 4 +
            0.00, 
            # [ 0, 1)
                # -1 * 0
            1 * 1 +
                # 1 * 0
            1 * 1 +
                # 2 * 0
            2 * 1 + 
            0.00, 
            # [ 1, 2)
                # -1 * -1
            1 * 4 + 
            0.00, 
            # [ 2, 3)
                # 1 * 2
            1 * 1 +
            0.00, 
            # [ 3, 4)
            0.00, 
            # [ 4, 5)
                # 2 * 2
            2 * 1 +
            0.00,
        ], sink_f32.data(), 7)


if __name__ == '__main__':
    gr_unittest.run(qa_histogram_solve)
