/*
Arithmetic coder based on "carryless rangecoder" by D. Subbotin.
* http://www.softcomplete.com/algo/pack/rus-range.txt
(Original C++ code by Subbotin)
* http://hem.spray.se/mikael.lundqvist/
(C impelementation by M. Lundqvist)
* H. Okumura: C magazine, Vol.14, No.7, pp.13-35, Jul. 2002
(Tutorial of data compression with sample codes, in japanese).
*/
#if defined(_MSC_VER)
#	pragma warning (disable : 4146) //disregard a warning C4146.
#endif

#include <stdio.h>
#include <stdlib.h>
#include "mrp.h"

RANGECODER *rc_init(void)
{
	RANGECODER *rc;

	rc = alloc_mem(sizeof(RANGECODER));
	rc->range = (range_t) -1;
	rc->low = rc->code = 0;
	return(rc);
}

void rc_encode(FILE *fp, RANGECODER *rc, uint cumfreq, uint freq, uint totfreq)
{
	rc->range /= totfreq;
	// if(rc->y == check_y && rc->x == check_x) printf("range\t%lld | totfreq\t%lld\n", rc->range, totfreq);
	rc->low += cumfreq * rc->range;
	// if(rc->y == check_y && rc->x == check_x) printf("low\t%lld | cumfreq\t%lld\n", rc->low, cumfreq);
	rc->range *= freq;
	// if(rc->y == check_y && rc->x == check_x) printf("range\t%lld | freq\t%lld\n", rc->range, freq);
	// if(rc->y == check_y && rc->x == check_x) printf("%lld\n%lld\n", (rc->low ^ (rc->low + rc->range)), RANGE_TOP);
	while((rc->low ^ (rc->low + rc->range)) < RANGE_TOP){
		putc(rc->low >> (RANGE_SIZE - 8), fp);
		// if(rc->y == check_y && rc->x == check_x)	printf("putc\t%lld\n", rc->low >> (RANGE_SIZE - 8));
		rc->code += 8;
		if (rc->code > 1E8) {
			fprintf(stderr, "Too large!\n");
			exit(1);
		}
		rc->range <<= 8;
		rc->low <<= 8;
		// if(rc->y == check_y && rc->x == check_x) printf("range\t%lld | low\t%lld | code\t%lld\n", rc->range, rc->range, rc->code);
	}
	while(rc->range < RANGE_BOT){
		putc(rc->low >> (RANGE_SIZE - 8), fp);
		// if(rc->y == check_y && rc->x == check_x)	printf("putc2\t%lld\n", rc->low >> (RANGE_SIZE - 8));
		rc->code += 8;
		if (rc->code > 1E8) {
			fprintf(stderr, "Too large!\n");
			exit(1);
		}
		rc->range = ((-rc->low) & (RANGE_BOT - 1)) << 8;
		rc->low <<= 8;
		// if(rc->y == check_y && rc->x == check_x) printf("range\t%lld | low\t%lld | code\t%lld\n", rc->range, rc->range, rc->code);
	}
	return;
}

void rc_finishenc(FILE *fp, RANGECODER *rc)
{
	int i;

	for (i = 0; i < RANGE_SIZE; i += 8) {
		putc(rc->low >> (RANGE_SIZE - 8), fp);
		rc->code += 8;
		rc->low <<= 8;
	}
	return;
}

int rc_decode(FILE *fp, RANGECODER *rc, PMODEL *pm, int min, int max)
{
	int i, j, k;
	range_t offset, totfreq, rfreq;

	offset = pm->cumfreq[min];
	totfreq = pm->cumfreq[max] - offset;
	rc->range /= (range_t)totfreq;
	rfreq = (rc->code - rc->low) / rc->range;
	if (rfreq >= totfreq) {
		fprintf(stderr, "Data is corrupted!\n");
		exit(1);
	}
	i = min;
	j = max - 1;
	while (i < j) {
		k = (i + j) / 2;
		// if(rc->y == check_y && rc->x == check_x)	printf("%d %d %d\n", k, i, j);
		if ((pm->cumfreq[k+1] - offset) <= rfreq) {
			i = k+1;
		} else {
			j = k;
		}
	}
	// if(rc->y == check_y && rc->x == check_x)	printf("(%d,%d)\n", i, j);
	rc->low += (pm->cumfreq[i] - offset) * rc->range;
	rc->range *= pm->freq[i];
	// if(rc->y == check_y && rc->x == check_x)	printf("low\t%lld | range\t%lld\n", rc->low, rc->range);
	while((rc->low ^ (rc->low + rc->range)) < RANGE_TOP){
		rc->code = (rc->code << 8) | getc(fp);
		rc->range <<= 8;
		rc->low <<= 8;
	}
	while(rc->range < RANGE_BOT){
		rc->code = (rc->code << 8) | getc(fp);
		rc->range = ((-rc->low) & (RANGE_BOT - 1)) << 8;
		rc->low <<= 8;
	}

	return (i);
}

void rc_startdec(FILE *fp, RANGECODER *rc)
{
	int i;

	for (i = 0; i < RANGE_SIZE; i += 8) {
		rc->code = (rc->code << 8) | getc(fp);
	}
	return;
}
