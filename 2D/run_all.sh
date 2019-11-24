#!/bin/bash

cd ./isotl/
mpactrun.sh x30z30.inp ~/mpact_odd_xs
cd ../

cd ./oddtl/
mpactrun.sh x30z30.inp ~/mpact_odd_xs
cd ../

cd ./momtl/
mpactrun.sh x30z30.inp ~/mpact_odd_xs
cd ../

cd ./sntl/
mpactrun.sh x30z30.inp ~/mpact_odd_xs
cd ../
