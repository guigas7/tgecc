#!/bin/bash
#indent -bap -bli0 -i4 -l79 -ncs -npcs -npsl -lc79 -fc1 -ts4 main.c
clear
gcc -o forkids kmain.c -lgmp
./forkids 
