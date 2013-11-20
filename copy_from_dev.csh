#!/bin/tcsh

rsync -a --copy-links --delete ../*.c ../*.h ../*.tex ../*.pdf ../Makefile* ../README* .
