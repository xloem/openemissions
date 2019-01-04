#pragma once



namespace freesdr
{

template <typename LinAlgImplementation, typename Element>
class Vector : public LinAlgImplementation::template Vector<Element>
{
  using VectorImplementation = class LinAlgImplementation::template Vector<Element>;
public:
  using VectorImplementation::VectorImplementation;
  using VectorImplementation::operator[];
};
 
}
