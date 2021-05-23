#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2021 Free Software Foundation, Inc..
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

from gnuradio import gr, gr_unittest, blocks
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
    if type(key) is str:
        key = gr.pmt.string_to_symbol(key)
    return gr.tag_utils.python_to_tag((offset, key, value))

def untag(tag):
    if tag.value.is_null():
        return (tag.offset, str(tag.key))
    else:
        return (tag.offset, str(tag.key), gr.pmt.to_double(tag.value))

class qa_accumulate(gr_unittest.TestCase):

    def setUp(self):
        self.tb = gr.top_block()
        self.len_key = 'packet_len'
        self.reset_key = 'reset_sum'
        self.src_data = blocks.vector_source_f(
            data = [
                0,1,2, # packet 1
                3,4,5,6,7, # packet 2
                8,9,10, # packet 3
                11,12,13 # packet 4
            ],
            tags = [
                tag(0, self.len_key, 3),
                tag(2, self.reset_key, 0.5),
                tag(3, self.len_key, 5),
                tag(8, self.len_key, 3),
                tag(10, self.reset_key),
                tag(11, self.len_key, 3)
            ]
        )
        self.accf32 = accumulate_f32(True, self.len_key, self.reset_key)

    def tearDown(self):
        self.tb = None

    def test_instance(self):
        instance_c64 = accumulate_c64()
        instance_c32 = accumulate_c32()
        instance_f64 = accumulate_f64()
        instance_f32 = accumulate_f32()
        instance_s64 = accumulate_s64()
        instance_int = accumulate_int()

    def test_001_descriptive_test_name(self):
        sink_data = blocks.vector_sink_f()
        self.tb.connect(self.src_data, self.accf32, sink_data)

        self.tb.run()

        for t in sink_data.tags():
            t = gr.tag_utils.tag_to_python(t)
        self.assertSequenceEqual([
            (0, self.len_key, 3.0),
            (2, self.reset_key, 0.5),
            (3, self.len_key, 5.0),
            (8, self.len_key, 5.0),
            (10, self.reset_key),
            (13, self.len_key, 5.0)
        ], [untag(t) for t in sink_data.tags()])
        np.testing.assert_almost_equal([
            0,1,2,
            (3+0/2)/1.5,(4+1/2)/1.5,(5+2/2)/1.5,6,7,
            (8+3+0/2)/2.5,(9+4+1/2)/2.5,(10+5+2/2)/2.5,0,0,
            11,12,13,0,0,
        ], sink_data.data(), 7)

if __name__ == '__main__':
    gr_unittest.run(qa_accumulate)
