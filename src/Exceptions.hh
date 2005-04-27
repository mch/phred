/* 
   phred - Phred is a parallel finite difference time domain
   electromagnetics simulator.

   Copyright (C) 2004-2005 Matt Hughes <mhughe@uvic.ca>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

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

class ResultException : public std::exception {
private:
  const char *buf_; 
public:
  ResultException(const char *buf) : buf_(buf) {}
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

class NoVtkException : public std::exception {
public:
  const char *what() const throw() { return "VTK is not availble."; }
};

class NoPythonException : public std::exception {
public:
  const char *what() const throw() { return "Python is not availble."; }
};

class StabilityException : public std::exception {
public:
  const char *what() const throw() { return "Delta x, y, z, and time settings do not statisfy Courant stability condition."; }
};

class PyInterpException : public std::exception {
private:
  const char *buf_;
public:
  PyInterpException(const char *buf) : buf_(buf) {}
  PyInterpException() : buf_(0) {}
  virtual const char *what() const throw() 
  { 
    if (buf_)
      return buf_; 
    else
      return "Python interpreter exception.";
  }
};

class FDTDException : public std::exception {
private:
  const char *buf_;
public:
  FDTDException(const char *buf) : buf_(buf) {}
  const char *what() const throw() { return buf_; }
};

class ParserException : public std::exception {
private:
  const char *buf_;
public:
  ParserException(const char *buf) : buf_(buf) {}
  const char *what() const throw() { return buf_; }
};

class BoundaryConditionException : public std::exception {
private:
  const char *buf_;
public:
  BoundaryConditionException(const char *buf) : buf_(buf) {}
  const char *what() const throw() { return buf_; }
};

class RegionException : public std::exception {
private:
  const char *buf_;
public:
  RegionException(const char *buf) : buf_(buf) {}
  const char *what() const throw() { return buf_; }
};

class UnknownMaterialException : public std::exception {
private:
  const char *buf_;
public:
  UnknownMaterialException(const char *buf) : buf_(buf) {}
  const char *what() const throw() { return buf_; }
};

class MaterialPropertyException : public std::exception {
private:
  const char *buf_;
public:
  MaterialPropertyException(const char *buf) : buf_(buf) {}
  const char *what() const throw() { return buf_; }
};

class CSGException : public std::exception {
private:
  const char *buf_;
public:
  CSGException(const char *buf) : buf_(buf) {}
  const char *what() const throw() { return buf_; }
};

#endif // EXCEPTIONS_H
