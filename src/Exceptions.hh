#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

/**
 * Exceptions for every occasion. Error "haiku" from
 * http://www.gnu.org/fun/jokes/error-haiku.html. 
 */

#include <exception>

class DimensionException : public std::exception {
public:
  const char *what() { return "Result must have at least one dimension."; }
};

class FileException : public std::exception {
public:
  const char *what() { return "Error using file."; }
};

class MemoryException : public std::exception {
public:
  const char *what() { return "Out of memory. \
We wish to hold the whole sky, \
But we never will."; }
};

class DataWriterException : public std::exception {
private:
  const char *buf_; 
public:
  DataWriterException(const char *buf) : buf_(buf) {}
  const char *what() { return buf_; }
};

#endif // EXCEPTIONS_H
