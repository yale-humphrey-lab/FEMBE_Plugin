#!/bin/bash

g++ -fPIC -shared FEMbeCmm.cpp dllmain.cpp -o FEMbeCmm.o -std=c++11 -I/home/eschwarz/FEBio/ -L/home/eschwarz/FEBio-build/lib -lfebiomech -lfecore
