#!/bin/bash

export GS_OPTIONS="-sPAPERSIZE=a4"
latex slides.tex
dvips slides.dvi -o slides.ps
ps2pdf slides.ps slides.pdf
