# Copyright 2011 Free Software Foundation, Inc.
#
# This file was generated by gr_modtool, a tool from the GNU Radio framework
# This file is a part of gr-openemissions
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

install(FILES
    openemissions_pigpio_sink.block.yml
    openemissions_pad_crop.block.yml
    openemissions_histogram.block.yml
    openemissions_tagged_stream_histogram.block.yml
    openemissions_accumulate.block.yml
    openemissions_histogram_solve.block.yml DESTINATION share/gnuradio/grc/blocks
)
