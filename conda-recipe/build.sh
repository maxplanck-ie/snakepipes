#!/bin/bash
# Installation essentially includes copying shared/ and workflows/ to $PREFIX/shared and making symlinks
mkdir -p $PREFIX/share/snakepipes
mkdir -p $PREFIX/bin
cp -r shared $PREFIX/share/snakepipes/
cp -r workflows $PREFIX/share/snakepipes/
for d in ls -d workflows/* ; do
    ln -s $PREFIX/share/snakepipes/$d/$d $PREFIX/bin/$d
done
