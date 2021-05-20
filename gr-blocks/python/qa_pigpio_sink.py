#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2021 Free Software Foundation, Inc..
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

import tracemalloc
tracemalloc.start()

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


# pigpiod is a daemon that runs on a raspberry pi to provide gpio access.
# to test this, a mock pigpiod is made below, using multiprocessing to
# run it in the background.

# this would be better replaced by an upstream patch to the pigpiod
# binary to have a mock mode itself, which would be simpler and provide
# much better testing.  such a patch is expected to be easy to implement.
# https://github.com/joan2937/pigpio/issues/469

# debugging multiprocess code in gdb:
# `set follow-fork-mode child` switches to the child process on fork
# `set detach-on-fork off` connects to all forks but might need tweaking
#                          to allocate process time properly

import select
import socket
import struct
import multiprocessing
import pigpio
from itertools import accumulate
class mock_pigpiod(multiprocessing.Process):
    def __init__(self, ip = '127.0.0.1', port = 8888):
        super().__init__()
        self.ip = ip
        self.port = port
        self.server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)

        while self.port < 9388:
            try:
                self.server.bind((self.ip, self.port))
                break
            except OSerror:
                self.port += 1

        self.address = self.ip + ':' + str(self.port)

        self.server.close()
    def start(self):
        super().start()
        # give time for server to listen
        pigpio.time.sleep(0.1)
    def run(self):
        self.errors = []
        self.modes = {}
        self.wavesendmodes = set()
        self.wave = []
        self.waves = {}
        self.wavequeue = []
        self.curwave = None
        self.start_time = None
        self.pulsequeue = []
        self.pulses_sent = []

        self.server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.server.bind((self.ip, self.port))
        #print('mock pigpiod listening on', (self.ip, self.port))
        self.server.listen()
        socks = [self.server, self.server.accept()[0]]
        while len(socks) > 1 or self.pulsequeue:
            readable, writable, exceptional = select.select(socks, [], socks, 1)

            # transmit queued waves
            #print('pulsequeue', self.pulsequeue)
            #print('wavequeue', self.wavequeue)
            #print('pulses_sent', self.pulses_sent)
            while self.pulsequeue and self.us() >= self.start_time:
                nextpulse = self.pulsequeue[0]
                self.pulsequeue = self.pulsequeue[1:]
                # pulses in pulses_sent are converted to have onset rather than finish times
                self.pulses_sent.append((
                    nextpulse[0], nextpulse[1], self.start_time
                ))
                self.start_time = nextpulse[2]
            while not self.pulsequeue and self.wavequeue:
                #print('setting curwave')
                self.curwave = self.wavequeue[0]
                self.pulsequeue = [
                    (gpioOn, gpioOff, time + self.start_time)
                    for gpioOn, gpioOff, time
                    in self.waves[self.curwave]
                ]
                self.wavequeue = self.wavequeue[1:]
            if not self.pulsequeue:
                #if self.curwave:
                    #print('unsetting curwave')
                self.curwave = None
                if self.start_time:
                    self.pulses_sent.append([
                        0, 0, self.start_time
                    ])
                    self.start_time = None

            for s in readable:
                if s is self.server:
                    conn, addr = s.accept()
                    socks.append(conn)
                else:
                    cmd = b''
                    while len(cmd) < 16:
                        nextdata = s.recv(16 - len(cmd))
                        if len(nextdata) == 0:
                            s.close()
                            socks.remove(s)
                            break
                        cmd += nextdata
                    if len(cmd) < 16:
                        continue
                    cmd, p1, p2, ext = struct.unpack('<4L',cmd)
                    #print('received ', cmd, p1, p2, ext)
                    if cmd == pigpio._PI_CMD_NOIB: # notify
                        self.reply(s)
                    elif cmd == pigpio._PI_CMD_NC: # notify close
                        self.reply(s)
                    elif cmd == pigpio._PI_CMD_BR1: # read_bank_1
                        self.reply(s)
                    elif cmd == pigpio._PI_CMD_HC: # hardware_clock
                        self.reply(s)
                    elif cmd == pigpio._PI_CMD_MODES: # set_mode
                        self.modes[p1] = p2
                        self.reply(s)
                    elif cmd == pigpio._PI_CMD_WVNEW: # wave_add_new
                        self.wave = []
                        self.reply(s)
                    elif cmd == pigpio._PI_CMD_WVAG: # wave_add_generic
                        pulses = s.recv(ext)
                        pulses = [pulses[i:i+12] for i in range(0, len(pulses), 12)]
                        pulses = [[*struct.unpack('<3L', pulse)] for pulse in pulses]
                        #print('PI_CMD_WVAG', pulses)
                        for index, pulse in enumerate(pulses):
                            if index == 0:
                                continue
                            pulse[2] += pulses[index - 1][2]
                        #print('accumulated pulse delays', pulses)
                        self.wave.extend(pulses)
                        self.wave.sort(key = lambda pulse: pulse[2])
                        self.reply(s)
                    elif cmd == pigpio._PI_CMD_WVCAP: # wave_create_and_pad
                        for wavenum in range(len(self.waves)+1):
                            if wavenum in self.waves:
                                continue
                            if self.wave:
                                self.waves[wavenum] = self.wave
                                self.wave = []
                                #print('PI_CMD_WVCAP', wavenum)
                                self.reply(s, wavenum)
                            else:
                                self.reply(s, pigpio.PI_EMPTY_WAVEFORM)
                            break
                    elif cmd == pigpio._PI_CMD_WVTXM: # wave_send_using_mode
                        wave_id, mode = p1, p2
                        self.wavesendmodes.add(mode)
                        self.wavequeue.append(wave_id)
                        #print('PI_CMD_WVTXM', wave_id)
                        if not self.start_time:
                            self.start_time = self.us()
                        self.reply(s)
                    elif cmd == pigpio._PI_CMD_WVTAT: # wave_tx_at
                        if self.curwave is None:
                            self.reply(s, pigpio.NO_TX_WAVE)
                        elif self.curwave not in self.waves:
                            self.reply(s, pigpio.WAVE_NOT_FOUND)
                        else:
                            #print('PI_CMD_WVTAT', self.curwave)
                            self.reply(s, self.curwave)
                    elif cmd == pigpio._PI_CMD_WVDEL: # wave_delete
                        #print('PI_CMD_WVDEL', p1)
                        if p1 not in self.waves:
                            self.reply(s, pigpio.PI_BAD_WAVE_ID)
                        else:
                            self.waves[p1] = None
                            self.reply(s)
                    else:
                        [sock.close() for sock in socks]
                        raise Exception('mock pigpiod received unimplemented command', cmd, p1, p2, ext)
            for s in exceptional:
                socks.remove(s)
                s.close()
        #print('mock pigpiod on', (self.ip, self.port), 'closing')
        socks[0].close()
        if self.start_time:
            self.pulses_sent.append([
                0, 0, self.start_time
            ])
            self.start_time = None
    def us(self):
        return int(pigpio.time.time() * 1000000)
    def reply(self, sock, code = 0):
        if code < 0 or code == pigpio.WAVE_NOT_FOUND:
            self.errors.append(code)
        sock.sendall(struct.pack('<12xl', code))

class bg_flowgraph(multiprocessing.Process):
    def __init__(self, samples_per_sec, data, max_nitems, pin, level, address = '127.0.0.1:8888'):
        super().__init__()
        self.tb = gr.top_block ()
        self.address = address
        self.sample_rate = samples_per_sec
        self.src = blocks.vector_source_f (data)
        self.max_nitems = max_nitems
        self.pin = pin
        self.level = level
        self.throttle = blocks.throttle(gr.sizeof_float, self.sample_rate)
    def run(self):
        sink = pigpio_sink_float (
            self.sample_rate,
            pin = self.pin,
            level = self.level,
            address = self.address
        )
        sink.set_max_noutput_items(self.max_nitems)
        
        self.tb.connect (self.src, self.throttle, sink)

        self.tb.run ()

class qa_pigpio_sink (gr_unittest.TestCase):

    def setUp (self):
        pass

    def tearDown (self):
        pass

    def test_instance (self):
        pigpiod = mock_pigpiod()
        pigpiod.start()
        instance_float = pigpio_sink_float (1000, address = pigpiod.address)
        instance_int = pigpio_sink_int (1000, address = pigpiod.address)
        instance_short = pigpio_sink_short (1000, address = pigpiod.address)
        instance_byte = pigpio_sink_byte (1000, address = pigpiod.address)

    def test_001_descriptive_test_name (self):
        pigpiod = mock_pigpiod()
        bgfg = bg_flowgraph(samples_per_sec = 10000,
                            data = [0, 0.2, 0.5, 0.4, 0.3],
                            max_nitems = 3,
                            pin = 4,
                            level = 0.3,
                            address = pigpiod.address)
        bgfg.start()
        pigpiod.run()

        self.assertSequenceEqual(pigpiod.errors, [])
        self.assertSequenceEqual([*pigpiod.modes.items()], [(4, pigpio.OUTPUT)])

        start_time = pigpiod.pulses_sent[0][2]
        self.assertSequenceEqual([
            (gpioOn, gpioOff, (time - start_time) * bgfg.sample_rate / 1000000)
            for gpioOn, gpioOff, time
            in pigpiod.pulses_sent
        ], [
            (0, 4, 0), # 4 dropped at sample 0
            (4, 0, 2), # 4 raised at sample 2
            (0, 4, 4), # 4 dropped at sample 4
            (0, 0, 5)  # stream ends
        ])

if __name__ == '__main__':
    gr_unittest.run (qa_pigpio_sink)
