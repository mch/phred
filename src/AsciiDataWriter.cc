#include "AsciiDataWriter.hh"

AsciiDataWriter::AsciiDataWriter()
  : DataWriter(rank, size), open_(false)
{}

AsciiDataWriter::~AsciiDataWriter()
{
  if (open_)
    file_.close();
}

AsciiDataWriter::init()
{

}

AsciiDataWriter::deinit()
{

}

void AsciiDataWriter::add_variable(const Result &result)
{

}

void AsciiDataWriter::handle_data(Data &data)
{

}
