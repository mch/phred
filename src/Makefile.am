SUBDIRS=PythonBindings Results Excitations DataWriters Boundaries \
	Sources CSG

bin_PROGRAMS=phred

AM_CPPFLAGS=-I./Results -I./PythonBindings -I./Excitations\
	 -I./DataWriters -I./Boundaries -I./Sources -I./CSG

phred_cond = 
phred_LDFLAGS=
phred_LDADD=

# Conditional on having Python and boost.python
if PY_BINDINGS
phred_LDADD += PythonBindings/libpyphread.a
phred_LDFLAGS += @PYTHON_LDFLAGS@
endif

BUILT_SOURCES = jan-grammer.h jan-grammer.cc jan-grammer-flex.cc

AM_YFLAGS = -d

# phred_LDADD += @LEXLIB@ 
phred_LDADD += \
	Results/libresults.a \
	Excitations/libexcitations.a \
	DataWriters/libdatawriters.a \
	Boundaries/libboundaries.a \
	Sources/libsources.a \
	CSG/libcsg.a

#	jan-grammer.yy \
#	jan-grammer-flex.ll \
#	JanFDTD.hh \
#	JanFDTD.cc 

phred_SOURCES=\
	system.h \
	ProblemGeometry.hh \
	ProblemGeometry.cc \
	Region.hh \
	Region.cc \
	LifeCycle.hh \
	LifeCycle.cc \
	FDTD.hh \
	FDTD.cc \
	Geometry.hh \
	Geometry.cc \
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
#	Tests.hh \
#	Tests.cc 
