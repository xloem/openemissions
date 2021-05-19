#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2021 Free Software Foundation, Inc..
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

from gnuradio import gr, gr_unittest, blocks
try:
    from openemissions import pad_crop
except ImportError:
    import os
    import sys
    dirname, filename = os.path.split(os.path.abspath(__file__))
    sys.path.append(os.path.join(dirname, "bindings"))
    from openemissions import pad_crop

class qa_pad_crop (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()
        self.len_key = gr.pmt.string_to_symbol('packet_len')
        self.src = blocks.vector_source_i(
            data = [
                1,2,3, # packet 1
                4,5,6,7,8 # packet 2
            ],
            tags = [
                gr.tag_utils.python_to_tag((0, self.len_key, gr.pmt.from_long(3))),
                gr.tag_utils.python_to_tag((3, self.len_key, gr.pmt.from_long(5)))
            ]
        )

    def tearDown (self):
        self.tb = None
        self.src = None

    def test_instance(self):
        instance = pad_crop(gr.sizeof_gr_complex, 16)

    def test_pad_crop (self):
        block_pad_crop = pad_crop(gr.sizeof_int, 4, str(self.len_key), True, True)
        sink_pad_crop = blocks.vector_sink_i()

        self.tb.connect(self.src, block_pad_crop, sink_pad_crop)

        self.tb.run ()
        
        self.assertSequenceEqual (
            sink_pad_crop.data(), 
            (1,2,3,0,4,5,6,7)
        )

    def test_pad_only (self):
        block_pad_only = pad_crop(gr.sizeof_int, 4, str(self.len_key), True, False)
        sink_pad_only = blocks.vector_sink_i()

        self.tb.connect(self.src, block_pad_only, sink_pad_only)

        self.tb.run ()
        
        self.assertSequenceEqual (
            sink_pad_only.data(),
            (1,2,3,0,4,5,6,7,8)
        )

    def test_crop_only (self):
        block_crop_only = pad_crop(gr.sizeof_int, 4, str(self.len_key), False, True)
        sink_crop_only = blocks.vector_sink_i()

        self.tb.connect(self.src, block_crop_only, sink_crop_only)

        self.tb.run ()
        
        self.assertSequenceEqual (
            sink_crop_only.data(),
            (1,2,3,4,5,6,7)
        )


if __name__ == '__main__':
    gr_unittest.run(qa_pad_crop)
