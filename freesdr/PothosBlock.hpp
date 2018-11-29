#pragma once

#include <Pothos/Framework.hpp>

using namespace opensdr;

template <typename Block>
BlockRegistration<Block>::BlockRegistration()
{
  static Pothos::BlockRegister register(
      Block::blockCategory() + "/" + Block::blockName(), 
      []()->Block *
      {
        return new Block();
      }
      );
}
