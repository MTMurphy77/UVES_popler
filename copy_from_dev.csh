#!/bin/tcsh

rsync -a --copy-links --delete --progress ../*.c ../*.h ../*.tex ../*.pdf ../Makefile* ../README* .
