#!/bin/bash

echo "Running an xterm on node `hostname`"
xterm -hold -sl 1000 -e $*
exit 0
