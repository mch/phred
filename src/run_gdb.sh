#!/bin/bash

echo "Running GDB on node `hostname`"
xterm -sl 256 -e gdb $*
exit 0
