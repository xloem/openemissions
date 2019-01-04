#pragma once

#include <stdexcept>

#include <Pothos/Framework.hpp>

namespace freesdr
{

class StreamGraphImplementationPothos
{
public:
  template <typename Block>
  static void registerBlock()
  {
    static Pothos::BlockRegistry registration(
      "/" + Block::blockCategory() + "/" + Block::blockName(), 
      []()->Block *
      {
        return new Block();
      }
    );
  }

  template <typename _FreeSDRBlock, typename... _StreamTypes>
  class Block : public Pothos::Block
  {
  public:
    using FreeSDRBlock = _FreeSDRBlock;
    using LinAlgImplementation = typename FreeSDRBlock::LinAlgImplementation;

    Block(FreeSDRBlock & freesdrBlock)
    : freesdrBlock(freesdrBlock)
    {
      // pre-c++17 paramater pack iteration
      // the parameter pack will expand inside an initializer list, so a dummy integer array is used
      (void)(int[]){ 
        setupStream(_StreamTypes::direction, typeid(_StreamTypes::Type), _StreamTypes::chans)...
      };
    }
  
    void work()
    {
      callWork(std::make_index_sequence<sizeof...(_StreamTypes)>());
    }
  
    FreeSDRBlock& freesdrBlock;
    using StreamTypes = std::tuple<_StreamTypes...>;
    using FreeSDRBuffers = std::tuple<typename FreeSDRBlock::template Buffer<_StreamTypes>...>;
    FreeSDRBuffers freesdrBuffers;
    std::array<Pothos::BufferChunk, sizeof...(_StreamTypes)> pothosBuffers;
    //  plan to convert pothosBuffers to shared_ptr's of buffers (perhaps?)
    //  question is where they will be allocated!
    // shared pointers are good for allowing pothos to manage buffer lifetime when may be sent across complex graph
    // seems like should be allocated at input ... ??
    //
    // possible searches:
    //  - does vector assignment provide for reseating of pointer?
    //    [X} yes if move semantics are used
    //  - does vector allow us to take ownership of vector-allocated pointer?
    //    -> ... well, you can use move semantics to another vector
    //    -> additionally you can have the vector wholly wrap preallocated memory
    //  - can we pass a shared ptr of a wrong object to sharedbuffer to get destructor called?
    //    [X] yes !
    //  - does vector allow us to provide an allocator?
    //    -> arma provides specialized aligned allocators.  could be good to use them and provide managed buffers.
    //  - does pothos provide general buffer allocation, not associated with output port?
    //    -> roughly, yes.  the buffer object constructors do this, there's also BufferPool and BufferManager
    //  - what happens when buffer is posted?
    //    -> just added to a queue; an existing buffer does not appear to be deallocated or anything
    //  - is it acceptable to allocate a buffer for all outputs, and ignore if reused?
    
    // concerns:
    //  - buffer from rtl-sdr will not be wrapped in vector
    //  - desire to pass raw vector objects to work handler
    //  - support eventually: no memory allocation during long run
    //  - support eventually: no memory copies in 1-to-1 linked graphs
    //        - unless the client tells us directly when forwarding happens, we might have to iterate
    //          all ports in order to identify it.
    //        - sometime, at runtime, zerocopy may be possible or impossible for a given input port
    //        - wasteful copy could happen if user does not need to mutate data in an input port (or even discards it)
    //
    // will need to wrap buffers in vectors
    // will likely need to deallocate or reuse memory smoewhere different from initial use
    //        -> shared ptr is reasonable solution to this
    //        -> looks like buffer manager would do this
    //
    //        we could use manager/pool to manage our own buffers
    //        may need _another_ manager/pool to manage buffers coming in from elsewhere
    //              note pothos has mechanisms to handle these buffers
    //        so pothos has its own wrappers for buffers that handle cleanup and reuse
    //
    //        sounds like two separate concepts: manage the vectors, and manage the buffers
    //
    //        when a buffer comes in, we don't really know whether it is already wrapped in a vector or not
    //        let's use the vectors in two diffrent ways:
    //            1. allocation/deallocation for memory alignment etc
    //            2. a structure to pass to the user
    //
    // forwarding by client:
    //        pothos manages buffer objects with shared pointers and pools and whatnot that have lifetimes
    //        but when the client moves data to the output port, they might just std::move() it.
    //        the output buffer now contains the same pointer as the input buffer
    //        but the object that controls its lifetime has not been moved over.
    //
    //        on the output, we'll have a set of pointers, really, and lengths.
    //        ideally we'd like to be in control of allocating all these pointers
    //        but they could come from multiple spots
    //
    //        an alternative would be to output via method call
    //
    // possible end-user experiences:
    //    - get a list of vectors, reuse the input vectors, and std::move them to output buffers
    //            -> what if I need to generate a new channel of output data?
    //            - the output buffers can be preallocated, and you use the preallocated data rather than std::moving
    //           X wastes copies if input data is discarded
    //    - get a set of vectors that cannot be written to, acquire new vectors from the backend, and write to those
    //            -> to support zerocopy, a function might be needed to mutate an input buffer into an output buffer
    //            -> this might be the way to go? due to concerns?
    //               incoming vector is unwritable
    //               to write, we call func to get mutable vector, passing the input vector
    //               requires pointer misuse for O(1) time, but could do that later ...
    //            -> could make this easier by passing input port # and output port #, or even references to them.
    //    - annotate streams as connecting to other streams.  inputs will be forwarded or copied to outputs.
    //            -> no aux function needed, buf user needs to keep track of port #'s in registration
    //            -> work function may be more intuitive
    //         if is connected, forward to output.  otherwise offer output as new buffer reference. user can resize if needed, I suppose, for now?

    // user paths
    // outputBuffer = std::move(inputBuffer) <-- want to use shared_ptr to translate ownership, or a buffer manager
    // outputBuffer.resize(39) <-- want to ensure sizes are retained across calls to reduce allocations
    // outputBuffer[17] = 3 <-- want to provide a buffer into allocated memory at start
    //
    // So, we'll certainly want output buffers to already be allocated from something.
    // We'll want to detect if they were replaced, and free/reuse the allocation if so.
    // We'll want to reallocate vectors when moving them along, or to store their sizes.
    //
    // - pool of bytes? need aligned data?
    // maybe just let vector allocate itself =/ skip the pool
    // problem of posting buffers.  pothos wants a shared pointer.
    // only way to steal memory from a vector is move operator
    // likely need a shared ptr to a vector, which is moved in.
    // guess I need some kind of vector pool
    //
    // i have a tuple of bufferchunks
    // each could wrap a sharedbuffer or a managedbuffer
    // each starts empty!
    // I get to convert them to vectors before passing them to the work handler
    // so, I can provide vectors that do not resize.
    //
    // I can use buffers from the pool, or I can make my own pool tht uses aligned memory
    // - [ ] add TODO for making aligned pool
    //
    // client can use std::move to forward input buffers 
    // input buffers, when copied, need to come from somewhere
    // since they can be moved along, perhaps use buffers
    // ... aligned alocation can't be managed by buffers unless inside shared ptr, which is reasonable
    //
    // -> buffer primitive: sharedbuffer wrapping shared_ptr to vector
    //
    // we want access to both the buffer and the vector, so it'll be split into bufferchunks and shared_ptrs
    // 
    
    // LEFT OFF HERE
    // - [X] fill buffers tuple
    //  -> I'll let the constructor fill them with default content
    // - [ ] provide buffer allocator that produces a bufferchunk wrapping a vector shared_ptr
    //        -> include TODO to use a pool
    // - [ ] generate input buffers
    // - [ ] call func
    // - [ ] hand off output buffers

    // primary goal is to reuse vector buffers
    // in a pool or somesuch
    // unfortunately pothos' pool only uses default allocation
  
    inline int setupStream(StreamDirection direction, std::type_info & type, size_t chans)
    {
      if (direction == SD_INPUT)
      {
        setupInput(inputPortInfo().size(), {type, chans});
      }
      else if (direction == SD_OUTPUT)
      {
        setupOutput(outputPortInfo().size(), {type, chans});
      }
      else
      {
        throw std::invalid_argument("unsupported stream direction");
      }
      return 0;
    }

    template <std::size_t... idxs>
    void callWork(std::index_sequence<idxs...>)
    {
      fillBuffersBeforeWork(*this, 0, 0);
      freesdrBlock.work(std::forward<std::tuple_element<idxs, FreeSDRBuffers>::type>((std::get<idxs>(freesdrBuffers)))...);
      sendBuffersAfterWork(*this, 0, 0);
    }

  };


private:
  template <typename FreeSDRBlock, typename... StreamTypes>
  friend class Block;

  template <typename Block, std::size_t... streamIdxs>
  static void fillBuffersBeforeWork(Block & block, std::size_t inputIdx, std::size_t outputIdx);

  template <typename Block, std::size_t... streamIdxs>
  static void sendBuffersAfterWork(Block & block, std::size_t inputIdx, std::size_t outputIdx);
};

template <typename Block, std::size_t idx, std::size_t... rest>
void StreamGraphImplementationPothos::fillBuffersBeforeWork<Block, idx, rest...>(Block & block, std::size_t inputIdx, std::size_t outputIdx)
{
  using StreamType = typename std::tuple_element<idx, typename Block::StreamTypes>;
  using Element = typename StreamType::Type;
  //using Buffer = typename FreeSDRBlock::template Buffer<StreamType>;
  using Matrix = LinAlgImplementation::template Matrix<Element>;
  auto & freesdrBuf = std::get<idx>(block.freesdrBuffers);
  //auto & pothosBuf = pothosBuffers[idx];

  //pothosBuffers[idx].clear();
  if (decltype(vec)::direction == SD_INPUT)
  {
    auto * pothosPort = block.input(inputIdx);
    ++ inputIdx;

    if (StreamType::bufShareStream != -1)
    {
      auto & freesdrShared = std::get<StreamType::bufShareStream>(block.freesdrBuffers);
      auto & pothosShared = pothosBuffers[StreamType::bufShareStream];

      // input can share buffer with output
      if (pothosPort->buffer().unique())
      {
        // here the buffer is unique, so we can take it
        // writing to it will work fine
        pothosShared = std::move(pothosPort->takeBuffer());
        freesdrShared = std::move(Matrix::view(pothosShared, pothosShared.elements(), StreamType::chans);
        freesdrBuf = std::move(freesdrShared.reference());
      }
    }
    else
    {
      // input buffer gives a view of input port
      freesdrBuf = std::move(Buffer::view(pothosPort->buffer(), pothosPort->elements(), StreamType::chans));
    }
  }
  else
  {
    ++ outputIdx;

    // rewind output buffer

    freesdrBuf.unget(freesdrBuf);
  }
  fillBuffersBeforeWork<rest...>(inputIdx, outputIdx);
}

template <typename Block>
void StreamGraphImplementationPothos::fillBuffersBeforeWork<>(Block & block, std::size_t inputIdx, std::size_t outputIdx) {}

template <typename Block, std::size_t idx, std::size_t... rest>
void StreamGraphImplementationPothos::sendBuffersAfterWork<idx, rest...>(Block & block, std::size_t inputIdx, std::size_t outputIdx)
{
  using FreeSDRBuf = typename std::tuple_element<idx, FreeSDRBuffers>::type;
  auto & freesdrBuf = std::get<idx>(freesdrBuffers);
  auto & pothosBuf = pothosBuffers[idx];
  if (decltype(vec)::direction == SD_OUTPUT)
  {
    auto * pothosPort = this->output(outputIdx);
    ++ outputIdx;
    // pothos will want a shared pointer to the vector, to manage lifecycle/time
    if (&freesdrBuf[0] == pothosBuffers[idx])
    {
      // buffer was forwarded from input port
      pothosBuf.setElements(freesdrBuf.gotten());
    }
    else
    {
      // buffer allocated by linalg implementation
      // TODO: use a buffer pool / manager
      pothosBuf = std::move(Pothos::BufferChunk(std::move(Pothos::SharedBuffer(&freesdrBuf[0], freesdrBuf.gotten(), std::shared_ptr<FreeSDRBuf>(new FreeSDRBuf(std::move(freesdrBuf)))))));
    }
    pothosPort->postBuffer(std::move(pothosBuf));
  }
  else
  {
    ++ inputIdx;
  }
  freesdrBuf.reset();
  sendBuffersAfterWork<rest...>)(inputIdx, outputIdx);
}

template <typename FreeSDRBlock, typename... StreamTypes>
void StreamGraphImplementationPothos::Block<FreeSDRBlock, StreamTypes...>::sendBuffersAfterWork<>(std::size_t inputIdx, std::size_t outputIdx) {}

/*
class StreamTypeImplementation
{
  friend class BlockImplementation;
public:
  template <typename StreamType>
  StreamTypeImplementation()
  : _direction(StreamType::direction),
    _type(typeid(StreamType::Type)),
    _chans(StreamType::chans)
  { }

private:
  StreamDirection _direction;
  Pothos::DType _type;
  size_t _chans;
};
*/

}
