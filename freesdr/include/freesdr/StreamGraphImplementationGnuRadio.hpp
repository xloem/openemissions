#pragma once

#include <stdexcept>

#include <gr_block.h>

namespace freesdr
{

class StreamGraphImplementationGnuRadio
{

  template <typename _FreeSDRBlock, typename... _StreamTypes>
  class Block : public gr_block
  {
  public:
    using FreeSDRBlock = _FreeSDRBlock;
    using LinAlgImplementation = typename FreeSDRBlock::LinAlgImplementation;

    Block(FreeSDRBlock & freesdrBlock)
    : freesdrBlock(freesdrBlock)
    {
      // setupstreams, pack expansion in pothos implementation
    }

    int general_work(int noutput_items,
                     gr_vector_int &ninput_items,
                     gr_vector_const_void_start &input_items,
                     gr_vector_void_start &output_items)
    {
      for (int streamtype
    }

    inline int setupStream(StreamDirection
  };

};

}
