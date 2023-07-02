#!/bin/bash
MAPFILE=$1.yaml 
ARCHFILE=arch/$4.yaml 
PROBFILE=prob/attention.yaml
MACROFILE=$2
python parser.py map-raw/$MAPFILE map/$MAPFILE 
tileflow map/$MAPFILE $ARCHFILE $PROBFILE $MACROFILE > $3 2>&1
test $? -eq 0 || echo tileflow map/$MAPFILE $ARCHFILE $PROBFILE $MACROFILE > $3 2>&1 >> error.sh