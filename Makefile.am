SUBDIRS=src doc tests

EXTRA_DIST=phred.spec.in 

# Run doxygen if it is available
doc:
	if test ! -z "@DOXYGEN@"; then cd doc; @DOXYGEN@; else echo "Doxygen is not installed. Cannot generate documentation."; fi
