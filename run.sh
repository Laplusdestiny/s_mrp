#! /bin/sh

for image in  1 #1 2 3 4 5 6 7 8 9 10 11 12 13
do

	DIR="/rda2/DATABASE/TMW"

	FLAG="-o"$@""

	## IMAGE ##
	if   [ "$image" = 1 ]; then IMG=$DIR/airplane.pgm;		NAME="airplane"
	elif [ "$image" = 2 ]; then IMG=$DIR/baboon.pgm;		NAME="baboon"
	elif [ "$image" = 3 ]; then IMG=$DIR/balloon.pgm;		NAME="balloon"
	elif [ "$image" = 4 ]; then IMG=$DIR/barb.pgm;			NAME="barb"
	elif [ "$image" = 5 ]; then IMG=$DIR/barb2.pgm;		NAME="barb2"
	elif [ "$image" = 6 ]; then IMG=$DIR/camera.pgm;		NAME="camera"
	elif [ "$image" = 7 ]; then IMG=$DIR/couple.pgm;		NAME="couple"
	elif [ "$image" = 8 ]; then IMG=$DIR/goldhill.pgm;		NAME="goldhill"
	elif [ "$image" = 9 ]; then IMG=$DIR/lena.pgm;			NAME="lena"
	elif [ "$image" = 10 ]; then IMG=$DIR/lennagrey.pgm;		NAME="lennagrey"
	elif [ "$image" = 11 ]; then IMG=$DIR/noisesquare.pgm;		NAME="noisesquare"
	elif [ "$image" = 12 ]; then IMG=$DIR/peppers.pgm;		NAME="peppers"
	elif [ "$image" = 13 ]; then IMG=$DIR/shapes.pgm;		NAME="shapes"
	fi

	LOG="Log/$NAME-`date +%y%m%d_%H%M`.txt"
	MRP="Encoded_File/$NAME.mrp"
	PGM="Decoded_File/$NAME.pgm"

	echo $FLAG > $LOG
	echo `hostname` >> $LOG

	date | tee -a $LOG
	./ENCMRP $FLAG $IMG $MRP | tee -a $LOG
	date | tee -a $LOG
	./DECMRP $MRP $PGM | tee -a $LOG
	if cmp -s $IMG $PGM;
	then
		echo "OK!" | tee -a $LOG
		rm -f $PGM
	else
		echo "ERROR!" | tee -a $LOG
		rm -f $PGM
		exit
	fi
done
