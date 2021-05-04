#!/usr/bin/env python3

# Goal: Make a live plot of data from the radio as a function of a parameter.
# The user can explore the data, changing parameters and view via navigation.

# Parameter axes and planes:
# - Gymbal orientation (elevation, azimuth)
# - Radio tuning
# - Radio gain
# - Historical data [consolidate holes]
# - Frequency decomposition
# - Other parameters offered by interfaces, such as which antenna on the radio to select

# An axis can be cross-sectioned, accumulated and averaged, a spread and cropped for display, or spread and stretched for full display

# Data format proposal:
#  dense binary chunks, looked up by database indices

import sqlite3

class DB:
    def __init__(self, dbname):
        self.db = sqlite3.connect(dbname)
