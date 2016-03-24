#!/usr/bin/env bash

pigz -p $1 -dc $3 | head -n $(echo $2 | awk '{print $1*4}') | pigz -p 4 -9 > $4
#zcat $2 | head -n $(echo $1 | awk '{print $1*4}') | gzip > $3
