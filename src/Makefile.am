bin_PROGRAMS=phred

phred_cond = 
phred_LDFLAGS=
phred_LDADD=

# Conditional on having Python and boost.python
if PY_BINDINGS
phred_cond += python-bindings/PyExcitations.hh \
	python-bindings/PyResults.hh \
	python-bindings/PyGrid.hh \
	python-bindings/PyFDTD.hh \
	python-bindings/PyTypes.hh \
	python-bindings/PyInterpreter.hh \
	python-bindings/PyInterpreter.cc
#	python-bindings/PyMaterial.cc

phred_LDADD += ../thirdparty/boost.python/src/libboost_python.la
phred_LDFLAGS += @PYTHON_LDFLAGS@
endif

BUILT_SOURCES = jan-grammer.h

AM_YFLAGS = -d

phred_LDADD += @LEXLIB@

phred_SOURCES=\
	jan-grammer.yy \
	jan-grammer-flex.ll \
	system.h \
	LifeCycle.hh \
	LifeCycle.cc \
	FDTD.hh \
	FDTD.cc \
	Geometry.hh \
	Box.hh \
	Box.cc \
	SphereGeom.hh \
	SphereGeom.cc \
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
	SourceFunction.hh \
	AltiVec.hh \
	Excitation.hh \
	WindowedExcitation.hh \
	BartlettExcitation.hh \
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
	BlockResult.hh \
	BlockResult.cc \
	PointResult.hh \
	PointResult.cc \
	PlaneResult.hh \
	PlaneResult.cc \
	FarfieldResult.cc \
	FarfieldResult.hh \
	fdtd_complex.hh \
	fdtd_complex.cc \
	DataWriter.hh \
	DataWriter.cc \
	AsciiDataWriter.hh \
	AsciiDataWriter.cc \
	MatlabDataWriter.hh \
	MatlabDataWriter.cc \
	NetCDFDataWriter.hh \
	NetCDFDataWriter.cc \
	HdfDataWriter.hh \
	Hdf5DataWriter.hh \
	VtkDataWriter.hh \
	VtkDataWriter.cc \
	PyDataWriter.hh \
	PyDataWriter.cc \
	PointDFTResult.hh \
	PointDFTResult.cc \
	SourceTimeResult.hh \
	SourceTimeResult.cc \
	SourceDFTResult.hh \
	SourceDFTResult.cc \
	phred.cc $(phred_cond)

#	SpaceExcitation.cc 
