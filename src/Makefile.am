SUBDIRS=PythonBindings Results Excitations DataWriters \
	Signals CSG Boundaries

bin_PROGRAMS=phred

AM_CPPFLAGS=-I./Results -I./PythonBindings -I./Excitations\
	 -I./DataWriters -I./Boundaries -I./Signals -I./CSG

phred_cond = 
phred_LDFLAGS=
phred_LDADD=

# Conditional on having Python and boost.python
if PY_BINDINGS
phred_LDADD += PythonBindings/libpyphread.a
phred_LDFLAGS += @PYTHON_LDFLAGS@
endif

AM_YFLAGS = -d

phred_LDADD += \
	Results/libresults.a \
	Excitations/libexcitations.a \
	DataWriters/libdatawriters.a \
	Boundaries/libboundaries.a \
	Signals/libsignals.a \
	CSG/libcsg.a

phred_SOURCES=\
	system.h \
	ProblemGeometry.hh \
	ProblemGeometry.cc \
	LifeCycle.hh \
	LifeCycle.cc \
	FDTD.hh \
	FDTD.cc \
	Block.hh \
	Block.cc \
	CellSet.hh \
	Grid.cc \
	Grid.hh \
	GridInfo.hh \
	GridInfo.cc \
	GridPlane.hh \
	PlaneTiling.hh \
	FreqGrid.hh \
	FreqGrid.cc \
	Material.cc \
	MaterialLib.cc \
	Material.hh \
	MaterialLib.hh \
	Constants.hh \
	Globals.hh \
	AltiVec.hh \
	Types.hh \
	Types.cc \
	Exceptions.hh \
	Data.hh \
	SubdomainAlg.hh \
	SimpleSDAlg.cc \
	SimpleSDAlg.hh \
	Contiguous.hh \
	Interval.hh \
	Tests.hh \
	Tests.cc \
	phred.cc $(phred_cond)

