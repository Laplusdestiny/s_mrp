/***** Decoder *****/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mrp.h"

extern CPOINT dyx[];
extern double sigma_a[];
extern double qtree_prob[];
extern double zerocoef_prob[];
extern int mask_y[], mask_x[];
extern int win_sample[], win_dis[];

MASK *mask;
int ***tempm_array;

uint getbits(FILE *fp, int n)
{
	static int bitpos = 0;
	static uint bitbuf = 0;
	int x = 0;

	if (n <= 0) return (0);
	while (n > bitpos) {
		n -= bitpos;
		x = (x << bitpos) | bitbuf;
		bitbuf = getc(fp) & 0xff;
		bitpos = 8;
	}
	bitpos -= n;
	x = (x << n) | (bitbuf >> bitpos);
	bitbuf &= ((1 << bitpos) - 1);
	return (x);
}

DECODER *init_decoder(FILE *fp)
{
	DECODER *dec;
	int i;

	dec = (DECODER *)alloc_mem(sizeof(DECODER));
	if (getbits(fp, 16) != MAGIC_NUMBER) {
		fprintf(stderr, "Not a compressed file!\n");
		exit(1);
	}
	dec->version = getbits(fp, 8);
	dec->width = getbits(fp, 16);
	dec->height = getbits(fp, 16);
	dec->maxval = getbits(fp, 16);
	dec->num_class = getbits(fp, 6);
	dec->num_group = getbits(fp, 6);
	dec->max_prd_order = getbits(fp, 7);
	dec->num_pmodel = getbits(fp, 6) + 1;
	dec->coef_precision = getbits(fp, 4) + 1;
	dec->max_coef = (2 << dec->coef_precision);
	dec->pm_accuracy = getbits(fp, 3);
	dec->quadtree_depth = (getbits(fp, 1))? QUADTREE_DEPTH : -1;
	dec->maxprd = dec->maxval << dec->coef_precision;
	dec->predictor = (int **)alloc_2d_array(dec->num_class, dec->max_prd_order,
		sizeof(int));
	dec->num_nzcoef = (int *)alloc_mem(dec->num_class * sizeof(int));
	dec->nzconv = (int **)alloc_2d_array(dec->num_class, dec->max_prd_order, sizeof(int));
	dec->th = (int **)alloc_2d_array(dec->num_class, dec->num_group - 1,
		sizeof(int));
	dec->err = (int **)alloc_2d_array(dec->height, dec->width, sizeof(int));
	dec->ctx_weight = init_ctx_weight();
	if (dec->quadtree_depth > 0) {
		int x, y, xx, yy;
		yy = (dec->height + MAX_BSIZE - 1) / MAX_BSIZE;
		xx = (dec->width + MAX_BSIZE - 1) / MAX_BSIZE;
		for (i = dec->quadtree_depth - 1; i >= 0; i--) {
			dec->qtmap[i] = (char **)alloc_2d_array(yy, xx, sizeof(char));
			for (y = 0; y < yy; y++) {
				for (x = 0; x < xx; x++) {
					dec->qtmap[i][y][x] = 0;
				}
			}
			yy <<= 1;
			xx <<= 1;
		}
	}
	dec->class = (char **)alloc_2d_array(dec->height, dec->width,
		sizeof(char));
	  dec->mask = (char **)alloc_2d_array(dec->height, dec->width,
				       sizeof(char));
	if (dec->num_pmodel > 1) {
		dec->pm_idx = (int *)alloc_mem(dec->num_group * sizeof(int));
	} else {
		dec->pm_idx = NULL;
	}
	dec->spm.freq = alloc_mem((MAX_SYMBOL * 2 + 1) * sizeof(uint));
	dec->spm.cumfreq = &(dec->spm.freq[MAX_SYMBOL]);

	dec->mult_pm.freq = alloc_mem(((dec->maxval + 1) * 2 + 1) * sizeof(uint));
	dec->mult_pm.cumfreq = &(dec->mult_pm.freq[dec->maxval + 1]);

	dec->sigma = sigma_a;

	dec->mtfbuf = (int *)alloc_mem(dec->num_class * sizeof(int));
#if AUTO_PRD_ORDER
	dec->ord2mhd = (int *)alloc_mem(dec->max_prd_order * sizeof(int));
	for (i = 0; i < dec->max_prd_order; i++) {
		dec->ord2mhd[i] = (int)((sqrt(1 + 4 * i) - 1) / 2);
	}
	dec->prd_mhd = dec->ord2mhd[dec->max_prd_order - 1] + 1;
	dec->zero_fr = (int *)alloc_mem(NUM_ZMODEL * sizeof(int));
	for (i = 0; i < NUM_ZMODEL; i++) {
		dec->zero_fr[i] = (int)(zerocoef_prob[i] * (double)TOT_ZEROFR);
	}
#endif
	return (dec);
}

#if AUTO_PRD_ORDER

void decode_predictor(FILE *fp, DECODER *dec)
{
	int k, cl, coef, sgn, d, zero_m, coef_m;
	PMODEL *pm;
	int b;

	pm = &dec->spm;
	pm->size = dec->max_coef + NUM_ZMODEL + 5;
	pm->cumfreq[dec->max_coef + 5] = 0;
	for(k = dec->max_coef + 5; k < pm->size; k++) {
		pm->freq[k] = 1;
		pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
	}
	b = dec->max_coef + 2;
	for (d = 0; d < dec->prd_mhd; d++) {
		zero_m = rc_decode(fp, dec->rc, pm, dec->max_coef + 5, dec->max_coef + NUM_ZMODEL + 5)
			- (dec->max_coef + 5);
		coef_m = rc_decode(fp, dec->rc, pm, dec->max_coef + 5, dec->max_coef + 13)
			- (dec->max_coef + 5);
		pm->cumfreq[b] = 0;
		pm->freq[b] = TOT_ZEROFR - dec->zero_fr[zero_m];
		pm->freq[b + 1] = dec->zero_fr[zero_m];
		pm->cumfreq[b + 1] = pm->freq[b];
		pm->cumfreq[b + 2] = TOT_ZEROFR;
		set_spmodel(pm, dec->max_coef + 1, coef_m);
		for (k = d * (d + 1); k < (d + 1) * (d + 2); k++) {
			for (cl = 0; cl < dec->num_class; cl++) {
				coef = rc_decode(fp, dec->rc, pm, dec->max_coef + 2, dec->max_coef + 4)
					- (dec->max_coef + 2);
				if (coef == 1) {
					coef = rc_decode(fp, dec->rc, pm, 1, dec->max_coef + 1);
					sgn = rc_decode(fp, dec->rc, pm, dec->max_coef + 5, dec->max_coef + 7)
						- (dec->max_coef + 5);
					if (sgn) {
						coef = -coef;
					}
					dec->predictor[cl][k] = coef;
				} else {
					dec->predictor[cl][k] = 0;
				}
			}
		}
	}
	for (cl = 0; cl < dec->num_class; cl++) {
		d = 0;
		for (k = 0; k < dec->max_prd_order; k++) {
			if (dec->predictor[cl][k] != 0) {
				dec->nzconv[cl][d++] = k;
			}
		}
		dec->num_nzcoef[cl] = d;
	}
	return;
}

#else

void decode_predictor(FILE *fp, DECODER *dec)	//when AUTO_PRD_ORDER 0
{
	int k, m, cl, coef, sgn, d;
	PMODEL *pm;

	pm = &dec->spm;
	pm->size = dec->max_coef + 18;
	pm->cumfreq[dec->max_coef + 2] = 0;
	for(k = dec->max_coef + 2; k < pm->size; k++) {
		pm->freq[k] = 1;
		pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
	}
	for (k = 0; k < dec->max_prd_order; k++) {
		m = rc_decode(fp, dec->rc, pm, dec->max_coef + 2, dec->max_coef + 18) - (dec->max_coef + 2);
		set_spmodel(pm, dec->max_coef + 1, m);
		for (cl = 0; cl < dec->num_class; cl++) {
			coef = rc_decode(fp, dec->rc, pm, 0, dec->max_coef + 1);
			if (coef > 0) {
				sgn = rc_decode(fp, dec->rc, pm, dec->max_coef+2, dec->max_coef+4)
					- (dec->max_coef + 2);
				if (sgn) {
					coef = -coef;
				}
			}
			dec->predictor[cl][k] = coef;
		}
	}
	for (cl = 0; cl < dec->num_class; cl++) {
		d = 0;
		for (k = 0; k < dec->max_prd_order; k++) {
			if (dec->predictor[cl][k] != 0) {
				dec->nzconv[cl][d++] = k;
			}
		}
		dec->num_nzcoef[cl] = d;
	}
	return;
}

#endif

void decode_threshold(FILE *fp, DECODER *dec)
{
	int cl, gr, m, k;
	PMODEL *pm;

	pm = &dec->spm;
	pm->size = 16;
	pm->cumfreq[0] = 0;
	for (k = 0; k < pm->size; k++) {
		pm->freq[k] = 1;
		pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
	}
	m = rc_decode(fp, dec->rc, pm, 0, pm->size);
	set_spmodel(pm, MAX_UPARA + 2, m);
	for (cl = 0; cl < dec->num_class; cl++) {
		k = 0;
		for (gr = 1; gr < dec->num_group; gr++) {
			if (k <= MAX_UPARA) {
				k += rc_decode(fp, dec->rc, pm, 0, pm->size - k);
			}
			dec->th[cl][gr - 1] = k;
		}
	}

	if (dec->num_pmodel > 1) {
		pm->size = dec->num_pmodel;
		pm->freq[0] = 0;
		for (k = 0; k < pm->size; k++) {
			pm->freq[k] = 1;
			pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
		}
		for (gr = 0; gr < dec->num_group; gr++) {
			dec->pm_idx[gr] = rc_decode(fp, dec->rc, pm, 0, pm->size);
		}
	}
	return;
}

void decode_qtindex(FILE *fp, DECODER *dec, PMODEL *cpm,
					int tly, int tlx, int blksize, int width, int level)
{
	int i, cl, y, x, bry, brx, ctx;
	char **qtmap;
	PMODEL *pm;

	brx = (tlx + blksize < dec->width) ? (tlx + blksize) : dec->width;
	bry = (tly + blksize < dec->height) ? (tly + blksize) : dec->height;
	if (tlx >= brx || tly >= bry) return;
	if (level > 0) {
		ctx = 0;
		qtmap = dec->qtmap[level - 1];
		y = (tly / MIN_BSIZE) >> level;
		x = (tlx / MIN_BSIZE) >> level;
		if (y > 0) {
			if (qtmap[y - 1][x] == 1) ctx++;
			if (brx < width && qtmap[y - 1][x + 1] == 1) ctx++;
		}
		if (x > 0 && qtmap[y][x - 1] == 1) ctx++;
		ctx = ((level - 1) * 4 + ctx) << 1;

		pm = &dec->spm;
		i = rc_decode(fp, dec->rc, pm, ctx, ctx + 2) - ctx;

		if (i == 1) {
			qtmap[y][x] = 1;
			blksize >>= 1;
			decode_qtindex(fp, dec, cpm, tly, tlx,
				blksize, width, level - 1);
			decode_qtindex(fp, dec, cpm, tly, tlx + blksize,
				blksize, width, level - 1);
			decode_qtindex(fp, dec, cpm, tly + blksize, tlx,
				blksize, width, level - 1);
			decode_qtindex(fp, dec, cpm, tly + blksize, tlx + blksize,
				blksize, brx, level - 1);
			return;
		}
	}
	i = rc_decode(fp, dec->rc, cpm, 0, cpm->size);

	mtf_classlabel(dec->class, dec->mtfbuf, tly, tlx,
		blksize, width, dec->num_class);
	for (cl = 0; cl < dec->num_class; cl++) {
		if (dec->mtfbuf[cl] == i) break;
	}
	for (y = tly; y < bry; y++) {
		for (x = tlx; x < brx; x++) {
			dec->class[y][x] = cl;
		}
	}
	return;
}

void decode_class(FILE *fp, DECODER *dec)
{
	int i, j, x, y, blksize, level;
	PMODEL *pm, cpm[1];
	double p;
	int ctx;
	int qtree_code[QUADTREE_DEPTH << 2], mtf_code[MAX_CLASS];

	if (dec->quadtree_depth >= 0) {
		level = dec->quadtree_depth;
		blksize = MAX_BSIZE;
	} else {
		level = 0;
		blksize = BASE_BSIZE;
	}

	pm = &dec->spm;
	if (dec->quadtree_depth > 0) {
		set_spmodel(pm, 7, -1);
		for (ctx = 0; ctx < QUADTREE_DEPTH << 2; ctx++) {
			qtree_code[ctx] = rc_decode(fp, dec->rc, pm, 0, pm->size);
		}
	}
	set_spmodel(pm, PMCLASS_LEVEL, -1);
	for (i = 0; i < dec->num_class; i++) {
		mtf_code[i] = rc_decode(fp, dec->rc, pm, 0, pm->size);
		if (pm->cumfreq[pm->size] < MAX_TOTFREQ) {
			for (j = 0; j < pm->size; j++) {
				if (j < mtf_code[i]) {
					pm->freq[j] /= 2;
				} else {
					pm->freq[j] *= 2;
				}
				if (pm->freq[j] <= 0) pm->freq[j] = 1;
				pm->cumfreq[j + 1] = pm->cumfreq[j] + pm->freq[j];
			}
		}
	}
	/* set prob. models */
	if (level > 0) {
		pm->size = QUADTREE_DEPTH << 3;
		for (ctx = 0; ctx < QUADTREE_DEPTH << 2; ctx++) {
			i = qtree_code[ctx];
			p = qtree_prob[i];
			pm->freq[(ctx << 1) + 1] = (uint)(p * (1 << 10));
			p = 1.0 - p;
			pm->freq[(ctx << 1)] = (uint)(p * (1 << 10));
		}
		for (i = 0; i < pm->size; i++) {
			pm->cumfreq[i + 1] = pm->cumfreq[i] + pm->freq[i];
		}
	}
	cpm->size = dec->num_class;
	cpm->freq = (uint *)alloc_mem((cpm->size * 2 + 1) * sizeof(uint));
	cpm->cumfreq = &cpm->freq[cpm->size];
	cpm->cumfreq[0] = 0;
	for (i = 0; i < dec->num_class; i++) {
		p = exp(-log(2.0) * ((double)mtf_code[i] + 0.5)
			* PMCLASS_MAX/PMCLASS_LEVEL);
		cpm->freq[i] = (uint)(p * (1 << 10));
		if (cpm->freq[i] <= 0) cpm->freq[i] = 1;
		cpm->cumfreq[i + 1] = cpm->cumfreq[i] + cpm->freq[i];
	}

	for (i = 0; i < dec->num_class; i++) {
		dec->mtfbuf[i] = i;
	}
	for (y = 0; y < dec->height; y += blksize) {
		for (x = 0; x < dec->width; x += blksize) {
			decode_qtindex(fp, dec, cpm, y, x,
				blksize, dec->width, level);
		}
	}
	return;
}


int calc_udec(DECODER *dec, int y, int x)
{
	int rx, ry, k, u;
	int **err, *wt_p;

	u = 0;
	err = dec->err;
	wt_p = dec->ctx_weight;
	if (y > UPEL_DIST && x > UPEL_DIST && x <= dec->width - UPEL_DIST) {
		for (k = 0; k < NUM_UPELS; k++) {
			ry = y + dyx[k].y;
			rx = x + dyx[k].x;
			u += err[ry][rx] * (*wt_p++);
		}
	} else if (y == 0) {
		if (x == 0) {
			for (k = 0; k < NUM_UPELS; k++) {
				u += ((dec->maxval + 1) >> 2) * (*wt_p++);
			}
		} else {
			ry = 0;
			for (k =0; k < NUM_UPELS; k++) {
				rx = x + dyx[k].x;
				if (rx < 0) rx = 0;
				else if (rx >= x) rx = x - 1;
				u += err[ry][rx] * (*wt_p++);
			}
		}
	} else {
		if (x == 0) {
			for (k = 0; k < NUM_UPELS; k++) {
				ry = y + dyx[k].y;
				if (ry < 0) ry = 0;
				else if (ry >= y) ry = y - 1;
				rx = x + dyx[k].x;
				if (rx < 0) rx = 0;
				u += err[ry][rx] * (*wt_p++);
			}
		} else {
			for (k = 0; k < NUM_UPELS; k++) {
				ry = y + dyx[k].y;
				if (ry < 0) ry = 0;
				rx = x + dyx[k].x;
				if (rx < 0) rx = 0;
				else if (rx >= dec->width) rx = dec->width - 1;
				u += err[ry][rx] * (*wt_p++);
			}
		}
	}
	u >>= 6;
	if (u > MAX_UPARA) u = MAX_UPARA;
	return (u);
}

int calc_prd(IMAGE *img, DECODER *dec, int cl, int y, int x)
{
	int k, prd, prd_order, rx, ry, *coef_p, *nzc_p, i;

	prd_order = dec->num_nzcoef[cl];
	prd = 0;
	coef_p = dec->predictor[cl];
	nzc_p = dec->nzconv[cl];
	if (y == 0) {
		if (x == 0) {
			for (k = 0; k < prd_order; k++) {
				prd += coef_p[nzc_p[k]];
			}
			prd *= ((img->maxval + 1) >> 1);
		} else {
			ry = 0;
			for (k = 0; k < prd_order; k++) {
				i = nzc_p[k];
				rx = x + dyx[i].x;
				if (rx < 0) rx = 0;
				else if (rx >= x) rx = x - 1;
				prd += coef_p[i] * img->val[ry][rx];
			}
		}
	} else {
		if (x == 0) {
			for (k = 0; k < prd_order; k++) {
				i = nzc_p[k];
				ry = y + dyx[i].y;
				if (ry < 0) ry = 0;
				else if (ry >= y) ry = y - 1;
				rx = x + dyx[i].x;
				if (rx < 0) rx = 0;
				prd += coef_p[i] * img->val[ry][rx];
			}
		} else {
			for (k = 0; k < prd_order; k++) {
				i = nzc_p[k];
				ry = y + dyx[i].y;
				if (ry < 0) ry = 0;
				rx = x + dyx[i].x;
				if (rx < 0) rx = 0;
				else if (rx >= img->width) rx = img->width - 1;
				prd += coef_p[i] * img->val[ry][rx];
			}
		}
	}
        prd = CLIP(0, dec->maxprd, prd);
	return (prd);
}

#if MULT_PEAK_MODE
#if 0
void decode_mask(FILE *fp, DECODER *dec)
{
	int cl, tly, tlx, bry, brx, yy, xx, k, y, x, i, j, ty, tx;
	PMODEL *pm;

	pm = &dec->spm;
	pm->size = NUM_MASK;
	pm->freq[0] = 0;

	for (k = 0; k < pm->size; k++) {
		pm->freq[k] = 1;
		pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
	}

	for(y = 0; y < dec->height; y += WIN_BSIZE){
		for(x = 0; x < dec->width; x += WIN_BSIZE){
			k = 0;
			cl = dec->class[y][x];
			tly = y - (NUM_MASK-1);
			tlx = x - (NUM_MASK-1);
			bry = y + WIN_BSIZE + (NUM_MASK-1);
			brx = x + WIN_BSIZE + (NUM_MASK-1);
			if(tly < 0) tly = 0;
			if(tlx < 0) tlx = 0;
			if(bry > dec->height) bry = dec->height;
			if(brx > dec->width) brx = dec->width;
			for(yy = tly; yy < bry; yy++){
				for(xx = tlx; xx < brx; xx++){
					if (dec->class[yy][xx] != cl) {
						k = rc_decode(fp, dec->rc, pm, 0, pm->size);
						goto dec;
					}
				}
			}
dec:
			if (y + WIN_BSIZE > dec->height)
				ty = dec->height % WIN_BSIZE;
			else
				ty = WIN_BSIZE;
			if(x + WIN_BSIZE > dec->width)
				tx  = dec->width % WIN_BSIZE;
			else
				tx = WIN_BSIZE;

			for(i = 0; i < ty; i++){
				for(j = 0; j < tx; j++){
					dec->mask[y + i][x + j] = k;
				}
			}
		}
	}
	return;
}
#endif
void decode_mask(FILE *fp, DECODER *dec)
{
	int i, j, x, y, k, tlx, tly, brx, bry, blksize;
	PMODEL *pm, cpm[1];
	double p;
	int mtf_code[NUM_MASK];

	blksize = WIN_BSIZE;

	pm = &dec->spm;

	set_spmodel(pm, PMMASK_LEVEL, -1);
	for (i = 0; i < NUM_MASK; i++) {
		mtf_code[i] = rc_decode(fp, dec->rc, pm, 0, pm->size);
		if (pm->cumfreq[pm->size] < MAX_TOTFREQ) {
			for (j = 0; j < pm->size; j++) {
				if (j < mtf_code[i]) {
					pm->freq[j] /= 2;
				} else {
					pm->freq[j] *= 2;
				}
				if (pm->freq[j] <= 0) pm->freq[j] = 1;
				pm->cumfreq[j + 1] = pm->cumfreq[j] + pm->freq[j];
			}
		}
	}
	/* set prob. models */
	cpm->size = NUM_MASK;
	cpm->freq = (uint *)alloc_mem((cpm->size * 2 + 1) * sizeof(uint));
	cpm->cumfreq = &cpm->freq[cpm->size];
	cpm->cumfreq[0] = 0;
	for (i = 0; i < NUM_MASK; i++) {
		p = exp(-log(2.0) * ((double)mtf_code[i] + 0.5)
			* PMMASK_MAX/PMMASK_LEVEL);
		cpm->freq[i] = (uint)(p * (1 << 10));
		if (cpm->freq[i] <= 0) cpm->freq[i] = 1;
		cpm->cumfreq[i + 1] = cpm->cumfreq[i] + cpm->freq[i];
	}

	for (i = 0; i < NUM_MASK; i++) {
		dec->mtfbuf[i] = i;
	}
	for (tly = 0; tly < dec->height; tly += blksize) {
		for (tlx = 0; tlx < dec->width; tlx += blksize) {
			brx = (tlx + blksize < dec->width) ? (tlx + blksize) : dec->width;
			bry = (tly + blksize < dec->height) ? (tly + blksize) : dec->height;
			i = rc_decode(fp, dec->rc, cpm, 0, cpm->size);
			mtf_classlabel(dec->mask, dec->mtfbuf, tly, tlx,
				blksize, dec->width, NUM_MASK);
			for (k = 0; k < NUM_MASK; k++) {
				if (dec->mtfbuf[k] == i) break;
			}
	  		for (y = tly; y < bry; y++) {
				for (x = tlx; x < brx; x++) {
					dec->mask[y][x] = k;
				}
			}
		}
	}
	return;
}

void init_mask()
{
	mask = (MASK *)alloc_mem(sizeof(MASK));
	mask->weight = (int *)alloc_mem(MAX_PEAK_NUM * sizeof(int));
	mask->class = (char *)alloc_mem(MAX_PEAK_NUM * sizeof(char));
	mask->base = (int *)alloc_mem(MAX_PEAK_NUM * sizeof(int));
	mask->pm =(PMODEL **)alloc_mem(MAX_PEAK_NUM * sizeof(PMODEL *));
}

int set_mask_parameter(IMAGE *img, DECODER *dec,int y, int x, int u, int bmask, int shift)
{
	int ty, tx, cl, i, peak, sample,m_gr,m_prd,m_base,r_cl,r_prd;
	int count_cl[dec->num_class];
	int *th_p;

	r_prd = -1;
	sample = win_sample[(int)dec->mask[y][x]];
	for (cl = 0; cl < dec->num_class; cl++) {
		count_cl[cl]  = 0;
	}

	for (i = 0; i < sample; i++){
		ty = y + mask_y[i];
		tx = x + mask_x[i];
		if(ty < 0) ty = 0;
		else if(ty >= dec->height) ty = dec->height - 1;
		if(tx < 0) tx = 0;
		else if(tx >= dec->width) tx = dec->width - 1;
		cl = dec->class[ty][tx];
		count_cl[cl]++;		//マスク内のクラスの数のカウント
	}

	r_cl = dec->class[y][x]; //real_class
	for(cl = peak = 0; cl < dec->num_class; cl++){
		if (count_cl[cl]!=0){
			// mask->class[peak] = cl;
			mask->weight[peak] = ( (count_cl[cl] << W_SHIFT) / sample);
			th_p = dec->th[cl];
			for (m_gr = 0; m_gr < dec->num_group - 1; m_gr++) {
				if (u < *th_p++) break;
			}

			m_prd = calc_prd(img, dec, cl, y, x);
			if (cl == r_cl) r_prd = m_prd;
			m_base = (dec->maxprd - m_prd + (1 << shift) / 2) >> shift;
			mask->pm[peak] = dec->pmodels[m_gr][0] + (m_base & bmask);
			m_base >>= dec->pm_accuracy;
			mask->base[peak] = m_base;
			peak++;
		}
	}
	mask->num_peak = peak;	//ピークの数
	return(r_prd);
}


IMAGE *decode_image(FILE *fp, DECODER *dec)		//when MULT_PEEK_MODE 1
{
	int x, y, cl, gr, prd, u, e, p, bitmask, shift, base;
	int *th_p;
	IMAGE *img;
	PMODEL *pm;

	img = alloc_image(dec->width, dec->height, dec->maxval);

	bitmask = (1 << dec->pm_accuracy) - 1;
	shift = dec->coef_precision - dec->pm_accuracy;
	for (y = 0; y < dec->height; y++) {
		for (x = 0; x < dec->width; x++) {
			u = calc_udec(dec, y, x);

			if (dec->mask[y][x] == 0){
				cl = dec->class[y][x];
				th_p = dec->th[cl];
				for (gr = 0; gr < dec->num_group - 1; gr++) {
					if (u < *th_p++) break;
				}
				prd = calc_prd(img, dec, cl, y, x);
				base = (dec->maxprd - prd + (1 << shift) / 2) >> shift;
				pm = dec->pmodels[gr][0] + (base & bitmask);
				base >>= dec->pm_accuracy;
				p = rc_decode(fp, dec->rc, pm, base, base+dec->maxval+1)
					- base;
			}else{ //mult_peak

				prd = set_mask_parameter(img, dec, y, x, u, bitmask, shift);
				if (mask->num_peak == 1){ // single_peak
					base = mask->base[0];
					pm = mask->pm[0];
					p = rc_decode(fp, dec->rc, pm, base, base+dec->maxval+1)
						- base;
				}else{
					pm = &dec->mult_pm;
					set_pmodel_mult(pm,mask,dec->maxval+1);
					p = rc_decode(fp, dec->rc, pm, 0, dec->maxval+1);
				}
			}
			img->val[y][x] = p;
			prd >>= (dec->coef_precision - 1);
			e = (p << 1) - prd - 1;
			if (e < 0) e = -(e + 1);
			dec->err[y][x] = e;
		}
	}
	return (img);
}
#else

IMAGE *decode_image(FILE *fp, DECODER *dec)		//when MULT_PEEK_MODE 0
{
	int x, y, cl, gr, prd, u, e, p, mask, shift, base;
	int *th_p;
	IMAGE *img;
	PMODEL *pm;

	img = alloc_image(dec->width, dec->height, dec->maxval);

	mask = (1 << dec->pm_accuracy) - 1;
	shift = dec->coef_precision - dec->pm_accuracy;
	for (y = 0; y < dec->height; y++) {
		for (x = 0; x < dec->width; x++) {
			cl = dec->class[y][x];
			u = calc_udec(dec, y, x);
			th_p = dec->th[cl];
			for (gr = 0; gr < dec->num_group - 1; gr++) {
				if (u < *th_p++) break;
			}
			prd = calc_prd(img, dec, cl, y, x);
			base = (dec->maxprd - prd + (1 << shift) / 2) >> shift;
			pm = dec->pmodels[gr][0] + (base & mask);
			base >>= dec->pm_accuracy;
			p = rc_decode(fp, dec->rc, pm, base, base+dec->maxval+1)
				- base;
			img->val[y][x] = p;
			prd >>= (dec->coef_precision - 1);
			e = (p << 1) - prd - 1;
			if (e < 0) e = -(e + 1);
			dec->err[y][x] = e;
		}

	}
	return (img);
}

#endif
void write_pgm(IMAGE *img, char *filename)
{
	int i, j;
	FILE *fp;
	fp = fileopen(filename, "wb");
	fprintf(fp, "P5\n%d %d\n%d\n", img->width, img->height, img->maxval);
	for (i = 0; i < img->height; i++) {
		for (j = 0; j < img->width; j++) {
			putc(img->val[i][j], fp);
		}
	}
	fclose(fp);
	return;
}

int main(int argc, char **argv)
{
	int i;
	IMAGE *img;
	DECODER *dec;
	char *infile, *outfile;
	FILE *fp;

	cpu_time();
	setbuf(stdout, 0);
	infile = outfile = NULL;
	for (i = 1; i < argc; i++) {
		if (infile == NULL) {
			infile = argv[i];
		} else {
			outfile = argv[i];
		}
	}
	if (infile == NULL || outfile == NULL) {
		printf(BANNER"\n", 0.01 * VERSION);
		printf("usage: decmrp infile outfile\n");
		printf("infile:     Input file\n");
		printf("outfile:    Output file\n");
		exit(0);
	}
	fp = fileopen(infile, "rb");
	dec = init_decoder(fp);

	dec->rc = rc_init();
	rc_startdec(fp, dec->rc);

	decode_class(fp, dec);
	decode_predictor(fp, dec);
	decode_threshold(fp, dec);
	dec->pmodels = init_pmodels(dec->num_group, dec->num_pmodel,
		dec->pm_accuracy, dec->pm_idx, dec->sigma,
		dec->maxval + 1);
#if  MULT_PEAK_MODE
	decode_mask(fp, dec);
	init_mask();
#endif
	img = decode_image(fp, dec);
	fclose(fp);

	write_pgm(img, outfile);
	printf("cpu time: %.2f sec.\n", cpu_time());


#if LOG_PUT_OUT_DEC
#if defined(_WIN32)
	if (set_directory()) {
		fprintf(stderr, "check \"DIR\"!\n");
		exit(1);
	}
#endif
	print_predictor(dec->predictor, dec->max_prd_order, dec->num_class, dec->max_coef, outfile);
	print_threshold(dec->th, dec->num_group, dec->num_class, NULL, dec->pm_idx, outfile);
//	print_class(dec->class, dec->num_class, dec->height, dec->width, outfile);
	print_class_color(dec->class, dec->num_class, dec->height, dec->width, outfile);
	print_class_and_block(dec->class, dec->num_class, dec->qtmap, dec->quadtree_depth, dec->height, dec->width, outfile);
//  print_mask(dec->mask, dec->height, dec->width, outfile);
	print_amp_chara(dec->predictor, dec->max_prd_order, dec->num_class, dec->height, dec->width, outfile);
#endif
	return(0);
}
