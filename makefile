#Makefile

CC = gcc
ifeq ($(CC), icc)
  CFLAG = -O3 -ip -axsse4.1 -msse3 -g -unroll
else
  CFLAG = -O3 -funroll-loops -Wall -g -fopenmp
endif

ENCODER = ENCMRP
ENCOBJ = encmrp.o common.o rc.o log.o
DECODER = DECMRP
DECOBJ = decmrp.o common.o rc.o log.o


all : $(ENCODER) $(DECODER)

$(ENCODER) : $(ENCOBJ)
	$(CC) $(CFLAG) -o $@ $(ENCOBJ) -lm

$(DECODER) : $(DECOBJ)
	$(CC) $(CFLAG) -o $@ $(DECOBJ) -lm

.c.o :
	$(CC) $(CFLAG) -c $<

encmrp.o : encmrp.c mrp.h
decmrp.o : decmrp.c mrp.h
common.o : common.c mrp.h
rc.o : rc.c mrp.h
log.o : log.c common.c mrp.h

clean :
	rm -f $(ENCODER) $(DECODER) $(ENCOBJ) $(DECOBJ) core.* *~

job :
	sjob -h compute-0-20 run.sh &

cleanlog :
	rm -rf Log/*
