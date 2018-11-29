#pragma once

#include <string>

namespace opensdr
{

class Block;

using blockFactor_t = Block * 

void registerBlock(std::string name, blockMake_t

template <typename Block>
class BlockRegistration
{
public:
  BlockRegistration();
};

class Block
{
public:
protected:
  Block()
private:

};
 
}
