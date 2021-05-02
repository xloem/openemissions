#pragma once

// thoughts: I'm wrapping a linear algebra lib.
// I wrapped the pothos lib by having it provide an implementation of the BlockImplementation class.
// I'd actually rather allow for one program to have both pothos and gnuradio blocks, I think ...
// (although when would that ever be needed? ... maybe the point of it is to make it easy for me
//  to change my mind about something given my indecision)
// either way there's a little more flexibility if the implementation classes are different
// and more flexibility if they can be selected at runtime
//
// given the math takes advantage of compile-time optimizations, a lot of it is not going to be
// selected at runtime.
// but one could write implementations that use multiple backends, and try them out at runtime for
// e.g. profiling.
// to enable this we'd probably want some way of selecting the implementation in the client code.
// this means for one thing the implementations should have different class names.
// likely we'll want to take the implementation as a template parameter
//   alternatively we could provide the interface in each implementation class, but that would provide for
//   less code reuse.
//
// this doesn't work well with the block class, because it wraps the implementation class and already
// has template parameters for the channel types !
// I guess I'd want to add a template parameter for the implemented, and specify a default.

#include <string>
#include <tuple>

namespace freesdr
{

template <typename Block>
class BlockRegistration
{
public:
  BlockRegistration()
  {
    Block::StreamGraphImplementation::template registerBlock<Block>();
  }
};

/*
#include <limits>
class StreamType
{
public:
  template <typename T>
  StreamType(size_t index, size_t channels = 1)
  : _index(index),
    _chans(channels),
    _size(channels * sizeof(T)),
    _signed(std::numeric_limits<T>::is_signeda),
    _floating(std::numeric_limits<T>::is_iec559),
    _type(typeid(T))
  {}

  size_t index() const { return _index; }
  size_t chans() const { return _chans; }
  size_t size() const { return _size; }
  bool isSigned() const { return _signed; }
  bool isFloating() const { return _floating; }
  std::type_info & type() const { return _type; }

private:
  size_t _index;
  size_t _chans;
  size_t _size;
  bool _signed;
  bool _floating;
  std::type_info & _type;
};
*/

using StreamDirection = int;
constexpr StreamDirection SD_OUTPUT = true;
constexpr StreamDirection SD_INPUT = false;

// there's a special situation where when blocks are linked together end-to-end linearly,
// the buffers can be reused between them such that there is only 1 buffer for the whole link.
// my streams don't have a concept of index when they are first created: it's implicit in the
// order they are created.
// input and output streams use the same type class.
// to specify that an input stream will share buffer with an output stream, what is done?
// maybe i'lls pecify the index of the other stream as a reuse argument
// maybe only use it for input streams, and it'll send the buffer to that output stream

// do I need this reuse argument?
// I think I can just assume that buffers are forwarded when they are not reused.
// I'll comment it out

// bufShareStream is the index among all input/output stream types to share a buffer with
// should be the index of an output buffer, specified on the input type
template <typename _Type, StreamDirection _direction, size_t _chans = 1>//, int _bufShareStream = -1>
class StreamType
{
public:
  using Type = _Type;
  static constexpr StreamDirection direction = _direction;
  static constexpr size_t chans = _chans;
  static constexpr size_t size = _chans * sizeof(_Type);
//  static constexpr size_t bufShareStream = _bufShareStream;
};

template <typename Type, size_t chans = 1>
using InputStreamType = StreamType<Type, SD_INPUT, chans>;
template <typename Type, size_t chans = 1>
using OutputStreamType = StreamType<Type, SD_OUTPUT, chans>;

// man what form do buffers take
//
// end-user: I want a nice function signature that I can type easily and is clear.
//    'Buffer' had better be a name exposed in the main Block class
//    templated on something I have easy access to, such as the port number or type.
//
//    it would be nice to use the type.
//    it would also be nice to have the class definition in this header file, so that it is easy to look up how to use the buffer.

template <typename _Block, typename _StreamType>
class Buffer : private _Block::LinAlgImplementation::template Matrix<_StreamType::Type>
{
  friend class _Block::BlockImplementation;
public:
  using Block = _Block;
  using StreamType = _StreamType;
  using LinAlgImplementation = typename Block::LinAlgImplementation;
  using Matrix = typename LinAlgImplementation::template Matrix<StreamType::Type>;

  Matrix const & get(std::size_t count)
  {
    if (StreamType::direction == SD_OUTPUT && _got + count > this->sizeMajor())
    {
      // TODO: use buffer pool or somesuch
      if (_got == 0)
      {
        this->resizeUninit(count);
      }
      else
      {
        this->resizeCopy(_got + count);
      }
    }
    _getView = this->subview(_got, 0, _got+count, StreamType::chans);
    _got += count;
    return _getView;
  }

  void unget(Matrix & data, std::size_t offset = 0)
  {
    _got = (&(data[0]) - &((*this)[0])) / StreamType::chans + offset;
  }

  std::size_t gotten() const { return _got; }

  std::size_t available() const { return this->major() - _got; }

private:
  void operator=(Matrix && matrix)
  {
    _got = 0;
    Matrix::operator=(std::forward(matrix));
  }
  void reset()
  {
    _got = 0;
    Matrix::reset();
  }
  std::size_t _got;
  Matrix _getView;
};

template <typename _LinAlgImplementation, typename _StreamGraphImplementation, typename ..._StreamTypes>
class Block
{
  using StreamGraphImplementation = _StreamGraphImplementation;
  using LinAlgImplementation = _LinAlgImplementation;
  using BlockImplementation = typename StreamGraphImplementation::template Block<_StreamTypes...>;
  using StreamTypes = typename std::tuple<_StreamTypes...>

  friend BlockImplementation;

public:
  template <typename StreamType>
  using Buffer = freesdr::Buffer<Block, StreamType>;
  template <typename Element, std::size_t chans = 1>
  using InputBuffer = Buffer<StreamType<Element, SD_INPUT, chans>>;
  template <typename Element, std::size_t chans = 1>
  using OutputBuffer = Buffer<StreamType<Element, SD_OUTPUT, chans>>;


protected:
  Block()
  : platform(*this)
  {}

  virtual void work(Buffer<_StreamTypes> & ...) = 0;

private:
  BlockImplementation platform;
};
 
}

//
// okay, streaming
// pothos supports:
// - custom buffer forwarding
// - automatic buffer forwarding
// - memory copying
//
// I think gnuradio only supports memory copying with library-provided buffers for each stream connection.
// kind of a waste, produces tons of copying for complex graphs.
// the pothos approach is much closer to the ideal of allowing the graph to optimize down to a single calculation
//
// looks like gnuradio requires memory copying for each block.
// we'd like to default to zerocopy buffers.
//    so, input buffers come in.  they likely have a length.
//    be nice to allow the user to consume only some of them.
//
//    maybe stream buffer objects =/  I guess it's just one more object.
//
//    input:
//      - provide buffers matching specifications of class
//      - perhaps assume that all data is processed, but allow caller to 'unread' data
//
//    output:
//      - buffer reception ....
//        maybe we can provide member functions for writing to output ports
//        outputCopy<idx>(type *, len)
//        outputMove<idx>(type *, len)
//        outputGet<idx>(idx, len)
//        => these will become methods on output buffers
//
//        what are we likely to want to do?
//        if we're producing new data of a new type, we'll want some allocated buffer to store the data in
//        -> outputGet
//        if we're performing a required transformation on a stream, like filtering, we can just reuse the input buf
//        -> outputMove
//        if a stream has multiple transformations, a copy will be needed, but the user has no knowledge of that
//        -> outputCopy can move into outputMove
//        if the data requires multiple passes to generate the output, we'll need to store them side-by-side
//        -> outputGet
//        if we have raw device buffers from hardware, we'll first of all want to pass them to the system as-is, but second-of-all need to handle when they are finished, to reuse them in hardware 
//            pothos provides for this with SharedBuffers.  one can pass a shared pointer, which holds a destructor
//            to call when all references are closed.
//            could be easy to wrap this and accept a shared pointer
//
//        so maybe outputBuffer() and outputGetBuffer(), outputBufferFilled()
//
//          outputGetBuffer() will likely return a std::vector or somesuch
//          user will assign elements of the buffer
//          maybe they could pass the buffer back into outputBuffer() ?
//          outputBuffer would want to recognize that the buffer is the same, which it could do via pointer comparison, I suppose.
//          it's a little slow, but simple
//          we could allow for optimization later with two separate calls
//          outputUserBuffer()
//          outputSystemBuffer()
//          getSystemBuffer()
//
//          -> issue here that getSystemBuffer and outputUserBuffer can't be used together.
//          i guess that's okay, better document it
//
//          maybe nicer to just use output() and newBuffer() or somesuch
//          output takes:
//            - buffer allocated by system
//            - input port buffer, moved to output, copied automatically if there is reuse
//            - user-provided shared pointer buffer
//
//          what about gnuradio wrt reusing input buffers? gnuradio's input buffers are const =/
//          i guess maybe specify how the buffers will be used for the class
//          then gnuradio will know whether to duplicate the data into the output buffer for reuse or not
//      
//          okay, I've added to the stream specification another stream the buffer is shared with.
//          i think that'll work.
//          how does this affect output?
//          I guess when we output, we need not provide the buffer.  it will be an input buffer.
//          user can now handle the copy/move distinction by manually tracking which buffer to share
//          whoops!
//          this is a decision that shuold be made at runtime, not in the template =/
//          so the code will need to be able to decide to duplicate the buffer into the input
//          if the output it came from is wired to other inputs too
//
//          how would that work with threading?
//          we can't move the buffer on.  other inputs may use it.
//          i guess if there's mroe than 1 we'll need to copy it.
//
//          --
//          I've partly groked the above.
//          The new plan is possibly to pass a list of buffers to the work function; this reduces the
//          template juggling needed to write the work function signature.
//
//          buffer will need a way to acquire vector from other buffer, or from user
//          buffer may need to handle copying when needed, could be confusing
//
//          given that Vector() could possi
