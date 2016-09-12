#! /bin/sh
make clean
make
sjob run.sh &
sleep 10s
tail -f nohup.out
