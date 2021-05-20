#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2021 Free Software Foundation, Inc..
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

working code is in test.py

from gnuradio import gr, gr_unittest
from gnuradio import blocks
try:
    from openemissions import pigpio_sink_float, pigpio_sink_int, pigpio_sink_short, pigpio_sink_byte
except ImportError:
    import os
    import sys
    dirname, filename = os.path.split(os.path.abspath(__file__))
    sys.path.append(os.path.join(dirname, "bindings"))
    from openemissions import pigpio_sink

import select
import socket
import struct
import threading
import pigpio
class mock_pigpiod(socket.socket, threading.Thread):
    def __init__(self, ip = '127.0.0.1', port = 8888):
        socket.socket.__init__(self, socket.AF_INET, socket.SOCK_STREAM)
        threading.Thread.__init__(self)
        self.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.bind((ip, port))
        self.listen(2)
        self.start()
    def run(self):
        socks = [self]
        connected = False
        self.running = True
        while self.running:
            readable, writable, exceptional = select.select(socks, [], socks)
            for s in readable:
                if s is self:
                    connected = True
                    conn, addr = s.accept()
                    socks.append(conn)
                else:
                    cmd, p1, p2, ext = struct.unpack('<LLLL',s.recv(16))
                    print(pigpio.time.time(), 'received ', cmd, p1, p2, ext)
                    if cmd == pigpio._PI_CMD_NOIB:
                        # notify handle
                        print('sending reply')
                        s.sendall(struct.pack('12xI', 0))
                        print(pigpio.time.time(), 'sent')
                    else:
                        self.close()
                        raise Exception(cmd, p1, p2, ext)
            for s in exceptional:
                socks.remove(s)
                s.close()
                if connected and len(inputs) == 1:
                    break;
        self.close()


class qa_pigpio_sink (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()
        self.pigpiod = mock_pigpiod()

    def tearDown (self):
        self.pigpiod.running = False
        self.tb = None

    def test_instance (self):
        instance_float = pigpio_sink_float (1000)
        instance_int = pigpio_sink_int (1000)
        instance_short = pigpio_sink_short (1000)
        instance_byte = pigpio_sink_byte (1000)

    def test_001_descriptive_test_name (self):
        src = blocks.vector_source_f ([
            0, 0.2, 0.5, 0.4, 0.3
        ])
        sink = pigpio_sink_float (
            samp_rate = 1000,
            pin = 4,
            level = 0.3
        )
        
        self.tb.connect (src, sink)

        self.tb.run ()


if __name__ == '__main__':
    gr_unittest.run (qa_pigpio_sink)
