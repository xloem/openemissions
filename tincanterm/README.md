A simple interchange protocol for transmitting data over just one or two wires.

96 unique values are used to represent just enough characters for a simple ASCII terminal.  All other values are discarded.

The protocol is such:
1. Rest state sends FALSE or LOW.
2. When a character is available, HIGH is sent.
3. The subsequent 7 bits represent the value.
4. Rest state is returned to.
