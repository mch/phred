bin_PROGRAMS=phred

# Conditional on having Python and boost.python
if PY_BINDINGS
phred_cond = python-bindings/PyExcitations.cc
#	python-bindings/PyGrid.cc \
#	python-bindings/PyMaterial.cc

phred_LDADD = ../thirdparty/boost.python/src/libboost_python.la
phred_LDFLAGS = @PYTHON_LDFLAGS@
endif

phred_SOURCES=phred.cc \
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
	SourceTimeResult.hh \
	SourceTimeResult.cc \
	SourceDFTResult.hh \
	SourceDFTResult.cc $(phred_cond)
#	jan-grammer.yy \
#	jan-grammer-flex.ff \
#	lex.yy.cc

#	SpaceExcitation.cc 

#BUILT_SOURCES = jan-grammer.h

AM_YFLAGS = -d

