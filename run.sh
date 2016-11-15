#! /bin/sh

for image in  6 #1 2 3 4 5 6 7 8 9 10 11 12 13
		#21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
do

	DIR="/rda2/DATABASE/TMW"
	DIR2="/rda1/users/sumito/image/gray8bit"

	FLAG="-o $@"

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

# http://www.imagecompression.info
	elif [ "$image" = 21 ]; then IMG=$DIR2/artificial.pgm;		NAME="artificial"
	elif [ "$image" = 22 ]; then IMG=$DIR2/big_building.pgm;		NAME="big_building"
	elif [ "$image" = 23 ]; then IMG=$DIR2/big_tree.pgm;		NAME="big_tree"
	elif [ "$image" = 24 ]; then IMG=$DIR2/bridge.pgm;		NAME="bridge"
	elif [ "$image" = 25 ]; then IMG=$DIR2/cathedral.pgm;		NAME="cathedral"
	elif [ "$image" = 26 ]; then IMG=$DIR2/deer.pgm;		NAME="deer"
	elif [ "$image" = 27 ]; then IMG=$DIR2/fireworks.pgm;		NAME="fireworks"
	elif [ "$image" = 28 ]; then IMG=$DIR2/flower_foveon.pgm;		NAME="flower_foveon"
	elif [ "$image" = 29 ]; then IMG=$DIR2/hdr.pgm;		NAME="hdr"
	elif [ "$image" = 30 ]; then IMG=$DIR2/leaves_iso_200.pgm;		NAME="leaves_iso_200"
	elif [ "$image" = 31 ]; then IMG=$DIR2/leaves_iso_1600.pgm;		NAME="leaves_iso_1600"
	elif [ "$image" = 32 ]; then IMG=$DIR2/nightshot_iso_100.pgm;		NAME="nightshot_iso_100"
	elif [ "$image" = 33 ]; then IMG=$DIR2/nightshot_iso_1600.pgm;		NAME="nightshot_iso_1600"
	elif [ "$image" = 34 ]; then IMG=$DIR2/spider_web.pgm;		NAME="spider_web"
	elif [ "$image" = 35 ]; then IMG=$DIR2/zone_plate.pgm;		NAME="zone_plate"
	fi

	LOG="Log/$NAME-`date +%y%m%d_%H%M`.log"
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
