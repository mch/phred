bin_PROGRAMS=phred
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
	PML.hh \
	PML.cc \
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
	DataWriter.hh \
	AsciiDataWriter.hh \
	AsciiDataWriter.cc \
	jan-grammer.yy \
	jan-grammer-flex.ff \
	lex.yy.cc

#	SpaceExcitation.cc 

phred_LDADD=-lfl

BUILT_SOURCES = jan-grammer.h

AM_YFLAGS = -d

EXTRA_DIST=phred.spec.in
