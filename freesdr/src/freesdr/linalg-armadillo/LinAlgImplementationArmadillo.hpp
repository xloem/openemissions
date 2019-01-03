#include <armadillo>

namespace freesdr
{

class LinAlgImplementationArmadillo
{
public:
  template <typename Element>
  using Vector = arma::template Col<Element>;
  
  template <typename Element>
  Vector<Element> vectorView(Element * ptr, std::size_t count)
  {
    return {ptr, count, false, true};
  }

  template <typename Element>
  Vector<Element> vectorView(Element const * ptr, std::size_t count)
  {
    return {ptr, count};
  }
};

}
