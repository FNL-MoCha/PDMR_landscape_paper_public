#!/usr/bin/env bash


python3.9 /data/MoCha/patidarr/oncokb-annotator-1.1.0/MafAnnotator.py -i $1 -o $2 -t $3 -u http://ec2-34-234-234-196.compute-1.amazonaws.com:8080/oncokb/

chmod 775 $2
