#!/bin/sh
make clean
make
make job ; tail -f nohup.out
// sleep 10s
// tail -f nohup.out
