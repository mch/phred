EXTRA_DIST=phred.lsm.in phred.texinfo html latex

#  build and install the .info pages
info_TEXINFOS = phred.texinfo
phred_TEXINFOS = gpl.texinfo

#  install the man pages
man_MANS=phred.1

dist-hook:
	if test ! -z "@DOXYGEN@"; then @DOXYGEN@ Doxyfile; fi
