bin_PROGRAMS=phred

phred_cond = 
phred_cond_bs = 

# Conditional on having Python and boost.python
if PY_BINDINGS
phred_cond += python-bindings/PyExcitations.cc
#	python-bindings/PyGrid.cc \
#	python-bindings/PyMaterial.cc

phred_LDADD = ../thirdparty/boost.python/src/libboost_python.la
phred_LDFLAGS = @PYTHON_LDFLAGS@
endif

BUILT_SOURCES = jan-grammer.hh

AM_YFLAGS = -d

phred_SOURCES=phred.cc \
	jan-grammer.yy \
	jan-grammer-flex.ll \
	system.h \
	Grid.cc \
	Grid.hh \
	GridInfo.hh \
	GridInfo.cc \
	GridPlane.hh \
	FreqGrid.hh \
	FreqGrid.cc \
	Material.cc \
	MaterialLib.cc \
	Material.hh \
	MaterialLib.hh \
	Constants.hh \
	Gaussm.cc \
	Gaussm.hh \
	Excitation.hh \
	Excitation.cc \
	TimeExcitation.hh \
	TimeExcitation.cc \
	SpaceExcitation.hh \
	BoundaryCondition.hh \
	BoundaryCondition.cc \
	Ewall.hh \
	Ewall.cc \
	Hwall.hh \
	Hwall.cc \
	Pml.hh \
	Pml.cc \
	PmlCommon.hh \
	PmlCommon.cc \
	SubdomainBc.hh \
	SubdomainBc.cc \
	Types.hh \
	Exceptions.hh \
	Data.hh \
	SubdomainAlg.hh \
	SimpleSDAlg.cc \
	SimpleSDAlg.hh \
	Result.hh \
	PointResult.hh \
	PointResult.cc \
	PlaneResult.hh \
	PlaneResult.cc \
	DataWriter.hh \
	DataWriter.cc \
	AsciiDataWriter.hh \
	AsciiDataWriter.cc \
	NetCDFDataWriter.hh \
	NetCDFDataWriter.cc \
	HdfDataWriter.hh \
	Hdf5DataWriter.hh \
	PointDFTResult.hh \
	PointDFTResult.cc \
	SourceTimeResult.hh \
	SourceTimeResult.cc \
	SourceDFTResult.hh \
	SourceDFTResult.cc $(phred_cond)

#	SpaceExcitation.cc 
