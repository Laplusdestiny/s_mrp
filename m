#! /bin/sh
make clean
make
# make job
# sjob -v run.sh &
sjob -h c25 run.sh &
sleep 10s
# `sjob run.sh &` && `tail -f nohup.out`
tail -f nohup.out
