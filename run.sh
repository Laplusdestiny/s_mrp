#! /bin/sh

	TMW="/rda2/DATABASE/TMW"
	DIR2="/rda1/users/sumito/image/gray8bit"
	KODAK="/rda2/DATABASE/kodak/gray"
	ITE="/rda2/DATABASE/ITE_ARIB_NEW/B_Series_HD"

	FLAG="-o$@"

for image in 13 #1 2 3 4 5 6 7 8 9 10 11 12 13
		#21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
		#40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62
do
	## IMAGE ##
	if   [ "$image" = 1 ]; then IMG=$TMW/airplane.pgm;		NAME="airplane"
	elif [ "$image" = 2 ]; then IMG=$TMW/baboon.pgm;		NAME="baboon"
	elif [ "$image" = 3 ]; then IMG=$TMW/balloon.pgm;		NAME="balloon"
	elif [ "$image" = 4 ]; then IMG=$TMW/barb.pgm;			NAME="barb"
	elif [ "$image" = 5 ]; then IMG=$TMW/barb2.pgm;			NAME="barb2"
	elif [ "$image" = 6 ]; then IMG=$TMW/camera.pgm;		NAME="camera"
	elif [ "$image" = 7 ]; then IMG=$TMW/couple.pgm;		NAME="couple"
	elif [ "$image" = 8 ]; then IMG=$TMW/goldhill.pgm;		NAME="goldhill"
	elif [ "$image" = 9 ]; then IMG=$TMW/lena.pgm;			NAME="lena"
	elif [ "$image" = 10 ]; then IMG=$TMW/lennagrey.pgm;	NAME="lennagrey"
	elif [ "$image" = 11 ]; then IMG=$TMW/noisesquare.pgm;	NAME="noisesquare"
	elif [ "$image" = 12 ]; then IMG=$TMW/peppers.pgm;		NAME="peppers"
	elif [ "$image" = 13 ]; then IMG=$TMW/shapes.pgm;		NAME="shapes"

# http://www.imagecompression.info
	elif [ "$image" = 21 ]; then IMG=$DIR2/artificial.pgm;			NAME="artificial"
	elif [ "$image" = 22 ]; then IMG=$DIR2/big_building.pgm;		NAME="big_building"
	elif [ "$image" = 23 ]; then IMG=$DIR2/big_tree.pgm;			NAME="big_tree"
	elif [ "$image" = 24 ]; then IMG=$DIR2/bridge.pgm;				NAME="bridge"
	elif [ "$image" = 25 ]; then IMG=$DIR2/cathedral.pgm;			NAME="cathedral"
	elif [ "$image" = 26 ]; then IMG=$DIR2/deer.pgm;				NAME="deer"
	elif [ "$image" = 27 ]; then IMG=$DIR2/fireworks.pgm;			NAME="fireworks"
	elif [ "$image" = 28 ]; then IMG=$DIR2/flower_foveon.pgm;		NAME="flower_foveon"
	elif [ "$image" = 29 ]; then IMG=$DIR2/hdr.pgm;					NAME="hdr"
	elif [ "$image" = 30 ]; then IMG=$DIR2/leaves_iso_200.pgm;		NAME="leaves_iso_200"
	elif [ "$image" = 31 ]; then IMG=$DIR2/leaves_iso_1600.pgm;		NAME="leaves_iso_1600"
	elif [ "$image" = 32 ]; then IMG=$DIR2/nightshot_iso_100.pgm;	NAME="nightshot_iso_100"
	elif [ "$image" = 33 ]; then IMG=$DIR2/nightshot_iso_1600.pgm;	NAME="nightshot_iso_1600"
	elif [ "$image" = 34 ]; then IMG=$DIR2/spider_web.pgm;			NAME="spider_web"
	elif [ "$image" = 35 ]; then IMG=$DIR2/zone_plate.pgm;			NAME="zone_plate"

# kodak
	elif [ "$image" = 40 ]; then IMG=$KODAK/kodak01.pgm;	NAME="kodak01"
	elif [ "$image" = 41 ]; then IMG=$KODAK/kodak02.pgm;	NAME="kodak02"
	elif [ "$image" = 42 ]; then IMG=$KODAK/kodak03.pgm;	NAME="kodak03"
	elif [ "$image" = 43 ]; then IMG=$KODAK/kodak04.pgm;	NAME="kodak04"
	elif [ "$image" = 44 ]; then IMG=$KODAK/kodak05.pgm;	NAME="kodak05"
	elif [ "$image" = 45 ]; then IMG=$KODAK/kodak06.pgm;	NAME="kodak06"
	elif [ "$image" = 46 ]; then IMG=$KODAK/kodak07.pgm;	NAME="kodak07"
	elif [ "$image" = 47 ]; then IMG=$KODAK/kodak08.pgm;	NAME="kodak08"
	elif [ "$image" = 48 ]; then IMG=$KODAK/kodak09.pgm;	NAME="kodak09"
	elif [ "$image" = 49 ]; then IMG=$KODAK/kodak10.pgm;	NAME="kodak10"
	elif [ "$image" = 50 ]; then IMG=$KODAK/kodak11.pgm;	NAME="kodak11"
	elif [ "$image" = 51 ]; then IMG=$KODAK/kodak12.pgm;	NAME="kodak12"
	elif [ "$image" = 52 ]; then IMG=$KODAK/kodak13.pgm;	NAME="kodak13"
	elif [ "$image" = 53 ]; then IMG=$KODAK/kodak14.pgm;	NAME="kodak14"
	elif [ "$image" = 54 ]; then IMG=$KODAK/kodak15.pgm;	NAME="kodak15"
	elif [ "$image" = 55 ]; then IMG=$KODAK/kodak16.pgm;	NAME="kodak16"
	elif [ "$image" = 56 ]; then IMG=$KODAK/kodak17.pgm;	NAME="kodak17"
	elif [ "$image" = 57 ]; then IMG=$KODAK/kodak18.pgm;	NAME="kodak18"
	elif [ "$image" = 58 ]; then IMG=$KODAK/kodak19.pgm;	NAME="kodak19"
	elif [ "$image" = 59 ]; then IMG=$KODAK/kodak20.pgm;	NAME="kodak20"
	elif [ "$image" = 60 ]; then IMG=$KODAK/kodak21.pgm;	NAME="kodak21"
	elif [ "$image" = 61 ]; then IMG=$KODAK/kodak22.pgm;	NAME="kodak22"
	elif [ "$image" = 62 ]; then IMG=$KODAK/kodak23.pgm;	NAME="kodak23"

# ITE ARIB NEW HD
	elif [ "$image" = 71 ]; then IMG=$ITE/s201_000060.pgm;	NAME="ite201"
	elif [ "$image" = 72 ]; then IMG=$ITE/s202_000060.pgm;	NAME="ite202"
	elif [ "$image" = 73 ]; then IMG=$ITE/s203_000060.pgm;	NAME="ite203"
	elif [ "$image" = 74 ]; then IMG=$ITE/s204_000060.pgm;	NAME="ite204"
	elif [ "$image" = 75 ]; then IMG=$ITE/s205_000060.pgm;	NAME="ite205"
	elif [ "$image" = 76 ]; then IMG=$ITE/s206_000060.pgm;	NAME="ite206"
	elif [ "$image" = 77 ]; then IMG=$ITE/s207_000060.pgm;	NAME="ite207"
	elif [ "$image" = 78 ]; then IMG=$ITE/s208_000060.pgm;	NAME="ite208"
	elif [ "$image" = 79 ]; then IMG=$ITE/s209_000060.pgm;	NAME="ite209"
	elif [ "$image" = 80 ]; then IMG=$ITE/s210_000060.pgm;	NAME="ite210"
	elif [ "$image" = 81 ]; then IMG=$ITE/s211_000060.pgm;	NAME="ite211"
	elif [ "$image" = 82 ]; then IMG=$ITE/s212_000060.pgm;	NAME="ite212"
	elif [ "$image" = 83 ]; then IMG=$ITE/s213_000060.pgm;	NAME="ite213"
	elif [ "$image" = 84 ]; then IMG=$ITE/s214_000060.pgm;	NAME="ite214"
	elif [ "$image" = 85 ]; then IMG=$ITE/s215_000060.pgm;	NAME="ite215"
	elif [ "$image" = 86 ]; then IMG=$ITE/s216_000060.pgm;	NAME="ite216"
	elif [ "$image" = 87 ]; then IMG=$ITE/s217_000060.pgm;	NAME="ite217"
	elif [ "$image" = 88 ]; then IMG=$ITE/s218_000060.pgm;	NAME="ite218"
	elif [ "$image" = 89 ]; then IMG=$ITE/s251_000060.pgm;	NAME="ite251"
	elif [ "$image" = 90 ]; then IMG=$ITE/s252_000060.pgm;	NAME="ite252"
	elif [ "$image" = 91 ]; then IMG=$ITE/s255_000060.pgm;	NAME="ite255"
	elif [ "$image" = 92 ]; then IMG=$ITE/s256_000060.pgm;	NAME="ite256"
	elif [ "$image" = 93 ]; then IMG=$ITE/s259_000060.pgm;	NAME="ite259"
	elif [ "$image" = 94 ]; then IMG=$ITE/s260_000060.pgm;	NAME="ite260"
	elif [ "$image" = 95 ]; then IMG=$ITE/s261_000060.pgm;	NAME="ite261"
	elif [ "$image" = 96 ]; then IMG=$ITE/s262_000060.pgm;	NAME="ite262"
	elif [ "$image" = 97 ]; then IMG=$ITE/s263_000060.pgm;	NAME="ite263"
	elif [ "$image" = 98 ]; then IMG=$ITE/s264_000060.pgm;	NAME="ite264"
	elif [ "$image" = 99 ]; then IMG=$ITE/s265_000060.pgm;	NAME="ite265"
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

	for file in `find Info/ -name "$NAME*.p*m" -type f`; do
		convert $file $file.png
		# echo $file converted $file.png
	done
	echo "Log file converted" | tee -a $LOG
done

