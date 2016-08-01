#! /bin/sh
make clean
make
make job
sleep 10s
tail -f nohup.out
