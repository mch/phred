EXTRA_DIST=phred.lsm.in phred.texinfo

#if WANT_DOCS
#EXTRA_DIST += html latex
#endif

#EXTRA_DIST += @DX_OUT_FILES@ @DX_OUT_DIRS@

#  build and install the .info pages
info_TEXINFOS = phred.texinfo
phred_TEXINFOS = gpl.texinfo

#  install the man pages
man_MANS=phred.1

#dist-hook:
#	if test ! -z "@DOXYGEN@"; then @DOXYGEN@ Doxyfile; fi
