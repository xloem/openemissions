#pragma once

#include <functional>

#include "types.hpp"
#include "recbufmeta.hpp"

template <typename Scalar, class DataSource, class DataProcessor>
class DataFeeder
{
public:
  DataFeeder(DataSource & source,
             DataProcessor & processor,
             std::function<bool()> checkDataInterruptedHook = [](){return false;})
  : _source(source),
    _processor(processor),
    _checkInterruptedHook(checkDataInterruptedHook)
  {
    _lastMeta.sampleTime = 0;
    _lastMeta.freq = -1;
    _lastMeta.rate = -1;
  }

  template <typename Derived>
  void add(Eigen::PlainObjectBase<Derived> const & chunk, RecBufMeta const & meta)
  {
    using Eigen::numext::mini;

    size_t minNeeded = _processor.bufferSizeMin();
    size_t multipleNeeded = _processor.bufferSizeMultiple();
    size_t maxNeeded = _processor.bufferSizeMax();

    size_t chunkOffset = 0;
    size_t nextSize;
    
    // exhaust buffer
    if (_lastMeta.sampleTime + _buffer.size() != meta.sampleTime ||
        _lastMeta.freq != meta.freq ||
        _lastMeta.rate != meta.rate)
    {
      // buffer doesn't contain enough samples, and next chunk is not contiguous with it
      // drop it
      _buffer.conservativeResize(0);
      _lastMeta = meta;
    }
    if (_buffer.size())
    {
      nextSize = mini((size_t)_buffer.size() + (size_t)chunk.size(), maxNeeded);
      if (multipleNeeded > 1)
      {
        nextSize -= nextSize % multipleNeeded;
      }
      if (nextSize < minNeeded)
      {
        _buffer.conservativeResize(_buffer.size() + chunk.size());
        _buffer.tail(chunk.size()) = chunk;
        return;
      }
      size_t chunkOffset = nextSize - _buffer.size();
      _buffer.conservativeResize(nextSize);
      _buffer.tail(chunkOffset) = chunk.head(chunkOffset);
      _processor.process(*this, _buffer, _lastMeta);
      _lastMeta.sampleTime += _buffer.size();
    }

    // pass segments of chunk
    while (!checkInterrupted())
    {
      minNeeded = _processor.bufferSizeMin();
      multipleNeeded = _processor.bufferSizeMultiple();
      maxNeeded = _processor.bufferSizeMax();
      nextSize = mini(chunk.size() - chunkOffset, maxNeeded);
      if (multipleNeeded > 1)
      {
        nextSize -= nextSize % multipleNeeded;
      }
      if (nextSize < minNeeded)
      {
        break;
      }

      _processor.process(*this, chunk.segment(chunkOffset, nextSize), _lastMeta);
      _lastMeta.sampleTime += nextSize;
      chunkOffset += nextSize;
    }

    nextSize = chunk.size() - chunkOffset;
    assert(_lastMeta.sampleTime == meta.sampleTime + chunk.size() - nextSize);
    _buffer.conservativeResize(nextSize);
    _buffer = chunk.tail(nextSize);
  }

  bool checkInterrupted() { return _checkInterruptedHook(); }

  DataSource & source() { return _source; }

  DataProcessor & processor() { return _processor; }

private:
  DataSource & _source;
  DataProcessor & _processor;
  HeapVector<Scalar> _buffer;
  std::function<bool()> _checkInterruptedHook;
  RecBufMeta _lastMeta;
};
