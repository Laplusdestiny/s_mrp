#!/bin/sh
make clean
make
make job25
sleep 1s
tail -f nohup.out
