#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

/**
 * Exceptions for every occasion. Error "haiku" from
 * http://www.gnu.org/fun/jokes/error-haiku.html. 
 */

#include <exception>

class DimensionException : public std::exception {
public:
  const char *what() const throw() { return "Result must have at least one dimension."; }
};

class FileException : public std::exception {
public:
  const char *what() const throw() { return "Error using file."; }
};

class MemoryException : public std::exception {
public:
  const char *what() const throw() { return "Out of memory. \
We wish to hold the whole sky, \
But we never will."; }
};

class DataWriterException : public std::exception {
private:
  const char *buf_; 
public:
  DataWriterException(const char *buf) : buf_(buf) {}
  virtual const char *what() const throw() { return buf_; }
};

class NoNetCDFException : public std::exception {
public:
  const char *what() const throw() { return "NetCDF is not availble."; }
};

class NoHDF4Exception : public std::exception {
public:
  const char *what() const throw() { return "HDF4 is not availble."; }
};

class NoHDF5Exception : public std::exception {
public:
  const char *what() const throw() { return "HDF5 is not availble."; }
};

class StabilityException : public std::exception {
public:
  const char *what() const throw() { return "Delta x, y, z, and time settings do not statisfy Courant stability condition."; }
};

#endif // EXCEPTIONS_H
