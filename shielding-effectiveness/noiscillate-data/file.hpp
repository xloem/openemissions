#pragma once

// TODO: finish converting to use unique ptr
#include <iostream>
class File
{
public:
  File(std::unique_ptr<std::iostream> && stream)
  : _stream(path, ios_base::in|ios_base::out|ios_base::binary)
  {
    stream.ignore(std::numeric_limits<std::streamsize>::max());
    length = stream.gcount();
    stream.clear();
    seek(0);
  }

  std::streamsize size() const { return length; }

  void seek(std::streamsize offset)
  {
    stream.seekg(offset);
    stream.seekp(offset);
  }

  std::streamsize tell()
  {
    auto ret1 = stream.tellg();
    auto ret2 = stream.tellp();
    return ret1 > ret2 ? ret1 : ret2;
  }

  template <typename T>
  void readReal(T & out)
  {
    checkAssign(out, _read<double>());
  }

  template <typename T>
  void writeReal(T const & in)
  {
    double in2;
    checkAssign(in2, in);
    _write(in2);
  }

  template <typename T>
  void readUInt(T & out)
  {
    uint8_t nums[8];
    _read(nums);
    uint64_t result =
      (nums[0] << 56ULL) |
      (nums[1] << 48ULL) |
      (nums[2] << 40ULL) |
      (nums[3] << 32ULL) |
      (nums[4] << 24ULL) |
      (nums[5] << 16ULL) |
      (nums[6] << 8ULL) |
      (nums[7]);
    checkAssign(out, result);
  }


  template <typename T>
  void writeUInt(T const & in)
  {
    uint64_t v;
    checkAssign(v, in);
    uint8_t nums[8] = {
      v >> 56ULL,
      v >> 48ULL,
      v >> 40ULL,
      v >> 32ULL,
      v >> 24ULL,
      v >> 16ULL,
      v >> 8ULL
    };
    _write(nums);
  }

  void readString(std::string & out)
  {
    decltype(std::string::size()) len;
    readUInt(len);
    out.resize(len);
    _read(out.data(), len);
  }

  void writeString(std::string const & in)
  {
    writeUInt(in.size());
    _write(in.data(), in.size());
  }

private:
  std::fstream stream;
  std::streamsize length;

  template <typename T> T _read() { T ret; _read(ret); return ret; }
  template <typename T> void _read(T &out) { _read(&out, sizeof(out)); }
  void _read(void *ptr, size_t ct) { stream.read(ptr, ct); }

  template <typename T> void _write(T const & in) { _write(&in, sizeof(in)); }
  void _write(void *ptr, size_t ct)
  {
    stream.write(ptr, ct);
    if (stream.tellp() > length)
    {
      length = stream.tellp();
    }
  }

  template <typename A, typename B>
  static void checkAssign(A & a, B const & b)
  {
    if (b < std::numeric_limits<A>::min)
    {
      throw std::runtime_error("value too negative for system conversion");
    }
    if (b > std::numeric_limits<A>::max)
    {
      throw std::runtime_error("value too large for system conversion");
    }
    a = b;
  }
};
