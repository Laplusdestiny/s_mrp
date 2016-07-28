#!/bin/sh
make clean
make
make job
sleep 1s
tail -f nohup.out
