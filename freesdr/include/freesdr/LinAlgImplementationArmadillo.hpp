#include <armadillo>

namespace freesdr
{

class LinAlgImplementationArmadillo
{
public:
  // can likely move this class out of this file if desired using a 'using' line.
  template <typename Element>
  class Vector : private arma::Col<Element>
  {
    using Col = arma::Col<Element>;
  public:
    Vector(std::size_t size = 0) : Col(size) {}

    std::size_t size() const { return this->n_cols; }

    void resizeUninit(std::size_t size) { this->set_size(size); }
    void resizeCopy(std::size_t size) { this->resize(size); }

    void reset() { this->operator=(std::move(Vector())); }

    Element & operator[](std::size_t idx) { return this->at(idx); }
    Element const & operator[](std::size_t idx) const { return this->at(idx); }

    Vector subview(std::size_t start, std::size_t end) { return this->subvec(start, end); }
    Vector subview(std::size_t start, std::size_t end) const { return this->subvec(start, end); }
    Vector reference() { return subview(0, size()); }
    Vector reference() const { return subview(0, size()); }
    using Col::head;
    using Col::tail;

    static Vector view(Element * ptr, std::size_t count)
    {
      return {ptr, count, false, true};
    }

    static Vector copy(Element const * ptr, std::size_t count)
    {
      return {const_cast<Element *>(ptr), count, true, false};
    }

  private:
    Vector(Element * ptr, std::size_t count, bool copy, bool no_realloc_on_resize)
    : Col(ptr, count, copy, no_realloc_on_resize)
    { }
  };

  // in armadillo, matrices are stored column-first
  // and row # is the first index
  template <typename Element>
  class Matrix : private arma::Mat<Element>
  {
    using Mat = arma::Mat<Element>;
  public:
    Matrix(std::size_t major = 0, std::size_t minor = 0) : Mat(minor, major) {}

    std::size_t sizeMajor() const { return this->n_cols; }
    std::size_t sizeMinor() const { return this->n_rows; }

    void resizeUninit(std::size_t major, std::size_t minor) { this->set_size(minor, major); }
    void resizeCopy(std::size_t major, std::size_t minor) { this->resize(minor, major); }

    void reset() { this->operator=(std::move(Matrix())); }

    Vector<Element> operator[](std::size_t major) { return this->unsafe_col(major); }
    Vector<Element> operator[](std::size_t major) const { return this->unsafe_col(major); }

    Element & operator()(std::size_t major, std::size_t minor) { return this->at(minor, major); }
    Element const & operator()(std::size_t major, std::size_t minor) const { return this->at(minor, major); }

    Matrix subview(std::size_t startMajor, std::size_t startMinor, std::size_t endMajor, std::size_t endMinor)
    {
      return this->submat(startMinor, startMajor, endMinor, endMajor);
    }
    Matrix subview(std::size_t startMajor, std::size_t startMinor, std::size_t endMajor, std::size_t endMinor) const
    {
      return this->submat(startMinor, startMajor, endMinor, endMajor);
    }

    Matrix reference() { return submat(0, 0, sizeMajor(), sizeMinor()); }

    static Matrix view(Element * ptr, std::size_t major, std::size_t minor)
    {
      return {ptr, major, minor, false, true};
    }

    static Matrix copy(Element const * ptr, std::size_t major, std::size_t minor)
    {
      return {const_cast<Element *>(ptr), major, minor, true, false};
    }

  private:
    Matrix(Element * ptr, std::size_t major, std::size_t minor, bool copy, bool no_realloc_on_resize)
    : Mat(ptr, minor, major, copy, no_realloc_on_resize)
    { }
  };
};

}
