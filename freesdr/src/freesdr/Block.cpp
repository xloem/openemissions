
using namespace freesdr;
using namespace std;

Block::Block(std::initializer_list<StreamType> inputStreams, std::initializer_list<StreamType> outputStreams)
: platform(*this, inputStreams, outputStreams)
{ }
