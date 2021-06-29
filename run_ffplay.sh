#!/bin/bash 

./demo -f 91.8e6 -s 2400000 -r 48000 - | ffplay -f s16le -ar 48000 -showmode 1 -i -


