#!/bin/sh

cp temp.sdf scratch/
cd scratch/
xtb temp.sdf --gfn 2 --opt > xtb.log 2>&1
cp xtbopt.sdf ../
