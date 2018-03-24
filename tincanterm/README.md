A protocol for transmitting data slowly and simply over very simple DIY connections consisting of just one- or two-bit communication channels.

96 unique values are used to represent just enough characters for a simple ASCII terminal.  All other values are discarded.

The intent of this project is to be as simple as possible, to make reviewing the code and system as easy as possible, and to allow
cheap communication from remote or highly isolated environments.

The protocol is such:
1. Rest state sends FALSE or LOW.
2. When a character is available, HIGH is sent.
3. The subsequent 7 bits represent the value.
4. Rest state is returned to.

For now, a single channel is unidirectional.
