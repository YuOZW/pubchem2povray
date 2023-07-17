#!/bin/sh

cp temp.pov scratch/
cd scratch/
povray temp.pov +W1280 +H720 +A +FN +Q11 +J > povray.log 2>&1
cp temp.png ../
