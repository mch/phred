SUBDIRS=python-bindings Results Excitations DataWriters Boundaries

bin_PROGRAMS=phred

phred_cond = 
phred_LDFLAGS=
phred_LDADD=

# Conditional on having Python and boost.python
if PY_BINDINGS
phred_LDADD += python-bindings/libpyphread.a
phred_LDFLAGS += @PYTHON_LDFLAGS@
endif

BUILT_SOURCES = jan-grammer.h

AM_YFLAGS = -d

phred_LDADD += @LEXLIB@ \
	Results/libresults.a \
	Excitations/libexcitations.a \
	DataWriters/libdatawriters.a \
	Boundaries/libboundaries.a

phred_SOURCES=\
	jan-grammer.yy \
	jan-grammer-flex.ll \
	system.h \
	LifeCycle.hh \
	LifeCycle.cc \
	FDTD.hh \
	FDTD.cc \
	JanFDTD.hh \
	JanFDTD.cc \
	Geometry.hh \
	Geometry.cc \
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
	ExpSine.hh \
	ExpSine.cc \
	SourceFunction.hh \
	AltiVec.hh \
	Types.hh \
	Exceptions.hh \
	Data.hh \
	SubdomainAlg.hh \
	SimpleSDAlg.cc \
	SimpleSDAlg.hh \
	fdtd_complex.hh \
	fdtd_complex.cc \
	counted_ptr.hh \
	phred.cc $(phred_cond)

#	SpaceExcitation.cc 
