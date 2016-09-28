#! /bin/sh
make clean
make
if [ `hostname` == 'itohws10' ] ; then
	sjob -v run.sh &
	sleep 10s
	tail -f nohup.out
else
	sh run.sh &
fi

