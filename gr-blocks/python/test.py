#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2021 Free Software Foundation, Inc..
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

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
import multiprocessing
import pigpio
from itertools import accumulate
class mock_pigpiod:
    def __init__(self, ip = '127.0.0.1', port = 8888):
        super().__init__()
        self.ip = ip
        self.port = port
        self.server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.server.bind((self.ip, self.port))
        self.server.listen(2)
    def run(self):
        self.errors = []
        self.modes = {}
        self.wavesendmodes = set()
        self.wave = []
        self.waves = {}
        self.wavequeue = []
        self.curwave = None
        self.pulsequeue = []
        self.pulses_sent = []
        socks = [self.server, self.server.accept()[0]]
        while len(socks) > 1:
            print('selecting')
            readable, writable, exceptional = select.select(socks, [], socks)
            print('got ', len(readable))

            # transmit queued waves
            print('pulsequeue', self.pulsequeue)
            print('wavequeue', self.wavequeue)
            while self.pulsequeue and self.us() >= self.pulsequeue[0][2]:
                nextpulse = self.pulsequeue[0]
                self.pulses_sent.append(nextpulse)
                self.start_time = nextpulse[2]
                self.pulsequeue = self.pulsequeue[1:]
            if not self.pulsequeue:
                if self.wavequeue:
                    print('setting curwave')
                    self.curwave = self.wavequeue[0]
                    self.pulsequeue = [
                        (gpioOn, gpioOff, time + self.start_time)
                        for gpioOn, gpioOff, time
                        in self.waves[self.curwave]
                    ]
                    self.wavequeue = self.wavequeue[1:]
                else:
                    print('unsetting curwave')
                    self.curwave = None

            for s in readable:
                if s is self.server:
                    conn, addr = s.accept()
                    socks.append(conn)
                else:
                    cmd = b''
                    while len(cmd) < 16:
                        nextdata = s.recv(16 - len(cmd))
                        if len(nextdata) == 0:
                            socks.remove(s)
                            break
                        print('recv', nextdata)
                        cmd += nextdata
                    if len(cmd) < 16:
                        continue
                    cmd, p1, p2, ext = struct.unpack('<4L',cmd)
                    #print(pigpio.time.time(), 'received ', cmd, p1, p2, ext)
                    if cmd in (pigpio._PI_CMD_NOIB, pigpio._PI_CMD_BR1, pigpio._PI_CMD_HC):
                        self.reply(s)
                    elif cmd == pigpio._PI_CMD_MODES:
                        self.modes[p1] = p2
                        self.reply(s)
                    elif cmd == pigpio._PI_CMD_WVNEW:
                        self.wave = []
                        self.reply(s)
                    elif cmd == pigpio._PI_CMD_WVAG:
                        pulses = s.recv(ext)
                        pulses = [pulses[i:i+12] for i in range(0, len(pulses), 12)]
                        pulses = [[*struct.unpack('<3L', pulse)] for pulse in pulses]
                        for index, pulse in enumerate(pulses):
                            if index > 0:
                                continue
                            pulse[2] += pulses[index - 1][2]
                        self.wave.extend(pulses)
                        self.wave.sort(key = lambda pulse: pulse[2])
                        self.reply(s)
                    elif cmd == pigpio._PI_CMD_WVCAP:
                        for wavenum in range(len(self.waves)+1):
                            if wavenum in self.waves:
                                continue
                            self.waves[wavenum] = self.wave
                            self.wave = []
                            #print('made wave', wavenum)
                            self.reply(s, wavenum)
                            break
                    elif cmd == pigpio._PI_CMD_WVTXM:
                        wave_id, mode = p1, p2
                        self.wavesendmodes.add(mode)
                        self.wavequeue.append(wave_id)
                        print('send wave', wave_id)
                        if not self.pulses_sent or not self.pulsequeue:
                            self.start_time = self.us()
                        self.reply(s)
                    elif cmd == pigpio._PI_CMD_WVTAT:
                        if self.curwave is None:
                            #print('no tx wave')
                            self.reply(s, pigpio.NO_TX_WAVE)
                        elif self.curwave not in self.waves:
                            #print('current wave not found')
                            self.reply(s, pigpio.WAVE_NOT_FOUND)
                        else:
                            #print('wave at', self.curwave)
                            self.reply(s, self.curwave)
                    elif cmd == pigpio._PI_CMD_WVDEL:
                        #print('delete wave', p1)
                        if p1 not in self.waves:
                            self.reply(s, pigpio.PI_BAD_WAVE_ID)
                        else:
                            self.waves[p1] = None
                            self.reply(s)
                    else:
                        self.server.close()
                        raise Exception(cmd, p1, p2, ext)
            for s in exceptional:
                socks.remove(s)
                s.close()
        self.server.close()
    def us(self):
        return int(pigpio.time.time() * 1000000)
    def reply(self, sock, code = 0):
        if code < 0 or code > 9000:
            self.errors.append(code)
        sock.sendall(struct.pack('<12xl', code))

class bg_flowgraph(multiprocessing.Process):
    def __init__(self, address = '127.0.0.1', port = 8888):
        super().__init__()
        self.tb = gr.top_block ()
        self.address = address
        self.port = port
        self.src = blocks.vector_source_f ([
            0, 0.2, 0.5, 0.4, 0.3
        ])
    def run(self):
        print('flowgraph starting')
        sink = pigpio_sink_float (
            samp_rate = 1000,
            pin = 4,
            level = 0.3,
            address = self.address + ":" + str(self.port)
        )
        print('made sink')
        
        self.tb.connect (self.src, sink)
        print('running')

        self.tb.run ()

class qa_pigpio_sink:

    def setUp (self):
        self.bgfg = bg_flowgraph()
        self.pigpiod = mock_pigpiod()

    def tearDown (self):
        self.pigpiod = None
        self.bgfg = None

    def test_instance (self):
        instance_float = pigpio_sink_float (1000)
        instance_int = pigpio_sink_int (1000)
        instance_short = pigpio_sink_short (1000)
        instance_byte = pigpio_sink_byte (1000)

    def test_001_descriptive_test_name (self):
        self.bgfg.start()
        self.pigpiod.run()

        print(len(self.pigpiod.errors), 0)
        print(self.pigpiod.modes.items(), (4, pigpio.OUTPUT))
        print(self.pigpiod.pulses_sent)



if __name__ == '__main__':
    a = qa_pigpio_sink()
    a.setUp()
    a.test_001_descriptive_test_name()
    a.tearDown()

