
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "mrp.h"
#include <omp.h>
#include <time.h>

extern CPOINT dyx[];
extern double sigma_a[];
extern double qtree_prob[];
extern double zerocoef_prob[];
extern int mask_x[], mask_y[];
extern int win_sample[], win_dis[];

MASK *mask;
int ***tempm_array;
int ***exam_array;

float ****calc_entropy_of_conditional_probability(PMODEL ***pmodels, int num_group,
	int num_pmodel, int pm_accuracy, int maxval)
{
	int i, j, k, l, gr, total;
	uint *freq_p, *cumfreq_p;
	PMODEL *pm_p, *pm;
	double *logfreq, logtotal, log2 = log(2.0);
	float **ptr1, *ptr2, ****c_prob, entropy;

	/* alloc 4D memory */
	c_prob = (float ****)alloc_2d_array(num_group, num_pmodel, sizeof(float **));
	ptr1 = (float **)alloc_mem(num_group * num_pmodel * (1 << pm_accuracy)
		* sizeof(float *));
	for (gr = 0; gr < num_group; gr++) {
		for (i = 0; i < num_pmodel; i++) {
			c_prob[gr][i] = ptr1;
			ptr1 += (1 << pm_accuracy);
		}
	}
	ptr2 = (float *)alloc_mem(num_group * num_pmodel * (1 << pm_accuracy)
		* (maxval + 1) * sizeof(float));
	for (gr = 0; gr < num_group; gr++) {
		for (i = 0; i < num_pmodel; i++) {
			for (j = 0; j < (1 << pm_accuracy); j++) {
				c_prob[gr][i][j] = ptr2;
				ptr2 += maxval + 1;
			}
		}
	}
	/* calc entropy of conditional probability */
	logfreq = (double *)alloc_mem(((maxval << 1) + 1) * sizeof(double));
	for (gr = 0; gr < num_group; gr++) {
		for (i = 0; i < num_pmodel; i++) {
			pm_p = pmodels[gr][i];
			for (j = 0; j < (1 << pm_accuracy); j++) {
				pm = pm_p + j;
				freq_p = pm->freq;
				cumfreq_p = pm->cumfreq;
				for (k = 0; k < (maxval << 1) + 1; k++) {
					logfreq[k] = log(freq_p[k]);
				}
				for (k = 0; k < maxval + 1; k++) {
					total = cumfreq_p[k + maxval + 1] - cumfreq_p[k];
					logtotal = log(total);
					entropy = 0;
					for (l = k; l < k + maxval + 1; l++) {
						entropy += (float)(freq_p[l] * (logtotal - logfreq[l]));
					}
					entropy /= (float)(total * log2);
					c_prob[gr][i][j][k] = entropy;
				}
			}
		}
	}
	return c_prob;
}

void calc_ratio_of_model_to_rate(ENCODER *enc)
{
	int y, x, prd, gr, e, base, frac, totpel, *gr_pel;
	double *entropy, *cost, *ratio, totent, totcost, totratio;
	PMODEL *pm;
	static int calc_entropy_flag = 0;
	static float ****entropy_of_conditional_probability;


	if (calc_entropy_flag == 0) {
		entropy_of_conditional_probability =
			calc_entropy_of_conditional_probability(enc->pmodels, enc->num_group,
			enc->num_pmodel,enc->pm_accuracy,
			enc->maxval);
		calc_entropy_flag = 1;
	}
	entropy = (double *)alloc_mem(enc->num_group * sizeof(double));
	cost = (double *)alloc_mem(enc->num_group * sizeof(double));
	ratio = (double *)alloc_mem(enc->num_group * sizeof(double));
	gr_pel = (int *)alloc_mem(enc->num_group * sizeof(int));
	for (gr = 0; gr < enc->num_group; gr++) {
		entropy[gr] = 0;
		cost[gr] = 0;
		gr_pel[gr] = 0;
	}
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			// cl = enc->class[y][x];
			gr = enc->group[y][x];
			e = enc->encval[y][x];
			prd = enc->prd[y][x];
			if (prd < 0) prd = 0;
			else if (prd > enc->maxprd) prd = enc->maxprd;
			base = enc->bconv[prd];
			frac = enc->fconv[prd];
			pm = enc->pmlist[gr] + frac;
			cost[gr] += pm->cost[base + e] + pm->subcost[base];
			entropy[gr] += entropy_of_conditional_probability[gr][pm->id][frac][base];
			gr_pel[gr]++;
		}
	}
	/* calc ratio */
	totpel = 0;
	totcost = totent = 0;
	for (gr = 0; gr < enc->num_group; gr++){
		if (entropy[gr] != 0.0)
			ratio[gr] = (1.0 - fabs(entropy[gr] - cost[gr]) / entropy[gr]);
		else
			ratio[gr] = 1.0;
		totent += entropy[gr];
		totcost += cost[gr];
		totpel += gr_pel[gr];
	}
	totratio = (1.0 - fabs(totent - totcost) / totent);
	/* print */
	printf("******* differences between entropy and rate *******\n");
	printf("(gr)  [shape]\tentropy\t\t|rate\t\t|fitness\t|pel\n");
	printf("------------------------------------------------------------------------------\n");
	for (gr = 0; gr < enc->num_group; gr++) {
		printf("(%2d)  [%.1f]\t%10d\t|%10d\t|   %.3f\t|%10d\n", gr,
			0.2 * (enc->pmlist[gr]->id + 1), (int)entropy[gr],
			(int)cost[gr], ratio[gr], gr_pel[gr]);
	}
	printf("------------------------------------------------------------------------------\n");
	printf("all         \t%10d\t|%10d\t|   %.3f\t|%10d\n", (int)totent, (int)totcost, totratio, totpel);
	free(entropy);
	free(cost);
	free(ratio);
	free(gr_pel);
}

IMAGE *read_pgm(char *filename)
{
	int i, j, width, height, maxval;
	char tmp[256];
	IMAGE *img;
	FILE *fp;

	fp = fileopen(filename, "rb");
	fgets(tmp, 256, fp);
	if (tmp[0] != 'P' || tmp[1] != '5') {
		fprintf(stderr, "Not a PGM file!\n");
		exit(1);
	}
	while (*(fgets(tmp, 256, fp)) == '#');
	sscanf(tmp, "%d %d", &width, &height);
	while (*(fgets(tmp, 256, fp)) == '#');
	sscanf(tmp, "%d", &maxval);

	if (maxval > 255) {
		fprintf(stderr, "Sorry, this version only supports 8bpp images!\n");
		exit(1);
	}
	img = alloc_image(width, height, maxval);
	for (i = 0; i < img->height; i++) {
		for (j = 0; j < img->width; j++) {
			img->val[i][j] = (img_t)fgetc(fp);
		}
	}
	fclose(fp);
	return (img);
}

int ***init_ref_offset(int height, int width, int prd_order, CPOINT *yx_array)
{
	int ***roff, *ptr;
	int x, y, dx, dy, k;
	int order, min_dx, max_dx, min_dy;

	min_dx = max_dx = min_dy = 0;
	order = (prd_order > NUM_UPELS)? prd_order : NUM_UPELS;
	for (k = 0; k < order; k++) {
		dy = yx_array[k].y;
		dx = yx_array[k].x;
		if (dy < min_dy) min_dy = dy;
		if (dx < min_dx) min_dx = dx;
		if (dx > max_dx) max_dx = dx;
	}
	roff = (int ***)alloc_2d_array(height, width, sizeof(int *));
	if (min_dy != 0) {
		ptr = (int *)alloc_mem((1 - min_dy) * (1 + max_dx - min_dx) * order * sizeof(int));
	}else {
		ptr = (int *)alloc_mem((2 + max_dx - min_dx) * order * sizeof(int));
	}
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			if (y == 0) {
				if (x == 0) {
					roff[y][x] = ptr;
					dx = 0;
					dy = height;
					for (k = 0; k < order; k++) {
						*ptr++ = dy * width + dx;
					}
				} else if (x + min_dx <= 0 || x + max_dx >= width) {
					roff[y][x] = ptr;
					dy = 0;
					for (k = 0; k < order; k++) {
						dx = yx_array[k].x;
						if (x + dx < 0) dx = -x;
						else if (dx >= 0) dx = -1;
						*ptr++ = dy * width + dx;
					}
				} else {
					roff[y][x] = roff[y][x - 1];
				}
				// for K = 1 and NUM_UPELS = 1
			} else if (min_dy == 0 && y == 1 && x == 0) {
				roff[y][x] = ptr;
				dy = -1;
				dx = 0;
				*ptr++ = dy * width + dx;
			} else if (y + min_dy <= 0) {
				if (x == 0) {
					roff[y][x] = ptr;
					for (k = 0; k < order; k++) {
						dy = yx_array[k].y;
						if (y + dy < 0) dy = -y;
						else if (dy >= 0) dy = -1;
						dx = yx_array[k].x;
						if (x + dx < 0) dx = -x;
						*ptr++ = dy * width + dx;
					}
				} else if (x + min_dx <= 0 || x + max_dx >= width) {
					roff[y][x] = ptr;
					for (k = 0; k < order; k++) {
						dy = yx_array[k].y;
						if (y + dy < 0) dy = -y;
						dx = yx_array[k].x;
						if (x + dx < 0) dx = -x;
						else if (x + dx >= width) {
							dx = width - x - 1;
						}
						*ptr++ = dy * width + dx;
					}
				} else {
					roff[y][x] = roff[y][x - 1];
				}
			} else {
				roff[y][x] = roff[y - 1][x];
			}
		}
	}
	return (roff);
}

ENCODER *init_encoder(IMAGE *img, int num_class, int num_group,
					  int prd_order, int coef_precision, int quadtree_depth,
					  int num_pmodel, int pm_accuracy)
{
	ENCODER *enc;
	int x, y, i, j;
	double c;
#if AUTO_PRD_ORDER
	cost_t *ptr;
	int k;
#endif

	enc = (ENCODER *)alloc_mem(sizeof(ENCODER));
	enc->height = img->height;
	enc->width = img->width;
	enc->maxval = img->maxval;
	enc->num_class = num_class;
	enc->num_group = num_group;
#if AUTO_PRD_ORDER
	enc->prd_order = BASE_PRD_ORDER;
	enc->max_prd_order = MAX_PRD_ORDER;
#else
	enc->prd_order = prd_order;
	enc->max_prd_order = prd_order;
#endif
	enc->coef_precision = coef_precision;
	enc->max_coef = (2 << coef_precision);
	enc->quadtree_depth = quadtree_depth;
	enc->num_pmodel = num_pmodel;
	enc->pm_accuracy = pm_accuracy;
	enc->maxprd = enc->maxval << enc->coef_precision;
	enc->predictor = (int **)alloc_2d_array(enc->num_class, enc->max_prd_order,
		sizeof(int));
	enc->num_nzcoef = (int *)alloc_mem(enc->num_class * sizeof(int));
	enc->optimize_loop = 0;
#if AUTO_PRD_ORDER
	for (i = 0; i < enc->num_class; i++) {
		enc->num_nzcoef[i] = BASE_PRD_ORDER;
	}
#else
	for (i = 0; i < enc->num_class; i++) {
		enc->num_nzcoef[i] = prd_order;
	}
#endif
	enc->nzconv = (int **)alloc_2d_array(enc->num_class, enc->max_prd_order, sizeof(int));
	for (i = 0; i < enc->num_class; i++) {
		for (j = 0; j < enc->max_prd_order; j++) {
			enc->nzconv[i][j] = j;
		}
	}
#if AUTO_PRD_ORDER
	enc->num_search = (int *)alloc_mem(enc->num_class * sizeof(int));
	for (i = 0; i < enc->num_class; i++) {
		enc->num_search[i] = 30;
	}
	enc->ord2mhd = (int *)alloc_mem(enc->max_prd_order * sizeof(int));
	for (i = 0; i < enc->max_prd_order; i++) {
		enc->ord2mhd[i] = (int)((sqrt(1 + 4 * i) - 1) / 2);
	}
	enc->prd_mhd = enc->ord2mhd[enc->max_prd_order - 1] + 1;
#endif
	enc->th = (int **)alloc_2d_array(enc->num_class, enc->num_group, sizeof(int));
	for (i = 0; i < enc->num_class; i++) {
		for (j = 0; j < enc->max_prd_order; j++) {
			enc->predictor[i][j] = 0;
		}
		for (j = 0; j < enc->num_group - 1; j++) {
			enc->th[i][j] = 0;
		}
		enc->th[i][enc->num_group - 1] = MAX_UPARA + 1;
	}

	enc->upara = (int **)alloc_2d_array(enc->height, enc->width, sizeof(int));
	enc->prd = (int **)alloc_2d_array(enc->height, enc->width, sizeof(int));
//***********
	enc->prd_class = (int ***)alloc_2d_array(enc->height, enc->width,
			sizeof(int *));
	enc->weight = (int ***)alloc_2d_array(enc->height, enc->width,
			 sizeof(int *));

	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			enc->prd_class[y][x] = (int *)alloc_mem(enc->num_class * sizeof(int));
			enc->weight[y][x] = (int *)alloc_mem(enc->num_class * sizeof(int));
		}
	}

//*************
	enc->roff = init_ref_offset(enc->height, enc->width, enc->max_prd_order, dyx);
	enc->org = (int **)alloc_2d_array(enc->height+1, enc->width, sizeof(int));
	enc->err = (int **)alloc_2d_array(enc->height+1, enc->width, sizeof(int));
	enc->cost = (cost_t **)alloc_2d_array(enc->height+1, enc->width, sizeof(cost_t));
	for(y=0; y<enc->height; y++){
		for(x=0; x<enc->width; x++){
			enc->cost[y][x] = 0;
		}
	}
	enc->cost[enc->height][0] = 0;
	if (enc->quadtree_depth > 0) {
		y = (enc->height + MAX_BSIZE - 1) / MAX_BSIZE;
		x = (enc->width + MAX_BSIZE - 1) / MAX_BSIZE;
		for (i = enc->quadtree_depth - 1; i >= 0; i--) {
			enc->qtmap[i] = (char **)alloc_2d_array(y, x, sizeof(char));
			y <<= 1;
			x <<= 1;
		}
	}
	enc->ctx_weight = init_ctx_weight();
	enc->class = (char **)alloc_2d_array(enc->height, enc->width, sizeof(char));
	enc->group = (char **)alloc_2d_array(enc->height, enc->width, sizeof(char));

	enc->mask = (char **)alloc_2d_array(enc->height, enc->width, sizeof(char));

	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			enc->group[y][x] = 0;
			enc->org[y][x] = img->val[y][x];
			enc->mask[y][x] = INIT_MASK;
		}
	}
	enc->org[enc->height][0] = (enc->maxval + 1) >> 1;
	enc->err[enc->height][0] = (enc->maxval + 1) >> 2;
	enc->uquant = (char **)alloc_2d_array(enc->num_class, MAX_UPARA + 1,
		sizeof(char));
	for (i = 0; i < enc->num_class; i++) {
		for (j = 0; j <= MAX_UPARA; j++) {
			enc->uquant[i][j] = enc->num_group - 1;
		}
	}
	enc->econv = (int **)alloc_2d_array(enc->maxval+1, (enc->maxval<<1)+1,
		sizeof(int));
	enc->bconv = (img_t *)alloc_mem((enc->maxprd + 1) * sizeof(img_t));
	enc->fconv = (img_t *)alloc_mem((enc->maxprd + 1) * sizeof(img_t));
	enc->pmlist = (PMODEL **)alloc_mem(enc->num_group * sizeof(PMODEL *));
	enc->spm.freq = alloc_mem((MAX_SYMBOL * 2 + 1) * sizeof(uint));
	enc->spm.cumfreq = &(enc->spm.freq[MAX_SYMBOL]);
	enc->mult_pm.freq = alloc_mem(((enc->maxval + 1) * 2 + 1) * sizeof(uint));
	enc->mult_pm.cumfreq = &(enc->mult_pm.freq[enc->maxval + 1]);

	enc->sigma = sigma_a;

	enc->mtfbuf = (int *)alloc_mem(enc->num_class * sizeof(int));

	enc->coef_m = (int *)alloc_mem( enc->prd_order * sizeof(int));
	for (i = 0; i < enc->prd_order; i++) {
		enc->coef_m[i] = 0;
	}
#if AUTO_PRD_ORDER
	enc->zero_m = (int *)alloc_mem(enc->prd_mhd * sizeof(int));
	for (j = 0; j < enc->prd_mhd; j++) {
		enc->zero_m[j] = NUM_ZMODEL >> 1;
	}
	enc->zero_fr = (int *)alloc_mem(NUM_ZMODEL * sizeof(int));
	for (i = 0; i < NUM_ZMODEL; i++) {
		enc->zero_fr[i] = (int)(zerocoef_prob[i] * (double)TOT_ZEROFR);
	}
	enc->coef_cost = (cost_t ***)alloc_2d_array(NUM_ZMODEL, 16, sizeof(cost_t *));
	ptr = (cost_t *)alloc_mem(NUM_ZMODEL * 16 * (enc->max_coef + 1) * sizeof(cost_t));
	for (i = 0; i < NUM_ZMODEL; i++) {
		for (j = 0; j < 16; j++) {
			enc->coef_cost[i][j] = ptr;
			ptr += (enc->max_coef + 1);
		}
	}
	for (k = 0; k < NUM_ZMODEL; k++) {
		for (i = 0; i < 16; i++) {
#if OPT_SIDEINFO
			double p, zero, nonz;
			uint cumb = 0;
			nonz = log((double)TOT_ZEROFR / (double)enc->zero_fr[k]) / log(2.0);
			zero = log((double)TOT_ZEROFR / ((double)TOT_ZEROFR - (double)enc->zero_fr[k])) / log(2.0);
			set_spmodel(&enc->spm, enc->max_coef + 1, i);
			cumb = enc->spm.freq[0];
			p = log(enc->spm.cumfreq[enc->max_coef + 1] - cumb);
			for (j = 1; j <= enc->max_coef; j++) {
				enc->coef_cost[k][i][j] = (p - log(enc->spm.freq[j])) / log(2.0);
				enc->coef_cost[k][i][j] += 1.0 + nonz;
			}
			enc->coef_cost[k][i][0] = zero;
#else
			for (j = 0; j <= enc->max_coef; j++) {
				enc->coef_cost[k][i][j] = 0;
			}
#endif
		}
	}
#else
	enc->coef_cost = (cost_t **)alloc_2d_array(16, enc->max_coef + 1,
		sizeof(cost_t));
	for (i = 0; i < 16; i++) {
#if OPT_SIDEINFO
		double p;
		set_spmodel(&enc->spm, enc->max_coef + 1, i);
		p = log(enc->spm.cumfreq[enc->max_coef + 1]);
		for (j = 0; j <= enc->max_coef; j++) {
			enc->coef_cost[i][j] = (p - log(enc->spm.freq[j])) / log(2.0);
			if (j > 0) enc->coef_cost[i][j] += 1.0;
		}
#else
		for (j = 0; j <= enc->max_coef; j++) {
			enc->coef_cost[i][j] = 0;
		}
#endif
	}
#endif
	enc->th_cost = (cost_t *)alloc_mem((MAX_UPARA + 2) * sizeof(cost_t));
	for (i = 0; i < MAX_UPARA + 2; i++) {
		enc->th_cost[i] = 0;
	}
	enc->class_cost = (cost_t *)alloc_mem(enc->num_class * sizeof(cost_t));
	c = log((double)enc->num_class) / log(2.0);
	for (i = 0; i < enc->num_class; i++) {
		enc->class_cost[i] = c;
	}
	for (i = 0; i < (QUADTREE_DEPTH << 3); i++) {
		enc->qtflag_cost[i] = 1.0;
	}
#if AUTO_DEL_CL
	i =  (enc->width + ( MAX_BSIZE >> QUADTREE_DEPTH) - 1) / ( MAX_BSIZE >> QUADTREE_DEPTH)
		* (enc->height + ((MAX_BSIZE >> QUADTREE_DEPTH) -1 ) / (MAX_BSIZE >> QUADTREE_DEPTH));
	enc->err_cost = (cost_t **)alloc_2d_array(i, enc->num_class, sizeof(cost_t));
#endif
	enc->cl_hist = (int *)alloc_mem(enc->num_class * sizeof(int));

#if TEMPLATE_MATCHING_ON
	enc->temp_num = (int **)alloc_2d_array(enc->height, enc->width, sizeof(int));
	enc->array = (int ***)alloc_3d_array(enc->height, enc->width, MAX_DATA_SAVE_DOUBLE, sizeof(int));
	init_3d_array(enc->array, enc->height, enc->width, MAX_DATA_SAVE_DOUBLE, 0);
	enc->w_gr = (int *)alloc_mem(enc->num_group * sizeof(int));
#endif
	return (enc);
}

void init_class(ENCODER *enc)
{
	int k, x, y, ly, lx, i, j, v, cl, sum, num_block;
	int *var, *tmp, **ptr;

	num_block = ((enc->height + BASE_BSIZE - 1) / BASE_BSIZE)
		* ((enc->width + BASE_BSIZE - 1) / BASE_BSIZE);
	var = (int *)alloc_mem(num_block * sizeof(int));
	ptr = (int **)alloc_mem(num_block * sizeof(int *));
	for (k = 0; k < num_block; k++) {
		y = (k / ((enc->width + BASE_BSIZE - 1) / BASE_BSIZE)) * BASE_BSIZE;
		x = (k % ((enc->width + BASE_BSIZE - 1) / BASE_BSIZE)) * BASE_BSIZE;
		ly = y + BASE_BSIZE;
		if (ly > enc->height)
			ly = enc->height;
		lx = x + BASE_BSIZE;
		if (lx > enc->width)
			lx = enc->width;
		var[k] = sum = 0;
		for (i = y; i < ly; i++) {
			for (j = x; j < lx; j++) {
				v = enc->org[i][j];
				sum += v;
				var[k] += v * v;
			}
		}
		var[k] -= sum * sum / ((ly -y) * (lx - x));
		ptr[k] = &(var[k]);
	}
	/* sort */
	for (i = num_block - 1; i > 0; i--) {
		for (j = 0; j < i; j++) {
			if (*ptr[j] > *ptr[j + 1]) {
				tmp = ptr[j];
				ptr[j] = ptr[j + 1];
				ptr[j + 1] = tmp;
			}
		}
	}
	for (k = 0; k < num_block; k++) {
		cl = (k * enc->num_class) / num_block;
		v = (int)(ptr[k] - var);
		y = (v / ((enc->width + BASE_BSIZE - 1) / BASE_BSIZE)) * BASE_BSIZE;
		x = (v % ((enc->width + BASE_BSIZE - 1) / BASE_BSIZE)) * BASE_BSIZE;
		ly = y + BASE_BSIZE;
		if (ly > enc->height) ly = enc->height;
		lx = x + BASE_BSIZE;
		if (lx > enc->width) lx = enc->width;
		for (i = y; i < ly; i++) {
			for (j = x; j < lx; j++) {
				enc->class[i][j] = cl;
			}
		}
	}
	free(ptr);
	free(var);
}

void init_mask()
{
	int peak_num = MAX_PEAK_NUM*2 + TEMPLATE_CLASS_NUM;
	mask = (MASK *)alloc_mem(sizeof(MASK));
	mask->weight = (int *)alloc_mem(peak_num * sizeof(int));
	mask->class = (char *)alloc_mem(peak_num * sizeof(char));
	mask->base = (int *)alloc_mem(peak_num * sizeof(int));
	mask->pm =(PMODEL **)alloc_mem(peak_num * sizeof(PMODEL *));
}

void set_cost_model(ENCODER *enc, int f_mmse)
{
	int gr, i, j, k;
	double a, b, c, var;
	int gauss_index;
	PMODEL *pm;

	/* temporary repairs */
	if( enc->num_pmodel == 1) {
		gauss_index = 0;
	}else if( enc->num_pmodel == 8 || enc->num_pmodel == 16 || enc->num_pmodel == 32 ) {
		gauss_index = ( 5 * enc->num_pmodel) / 8 - 1;
	}else {
		fprintf(stderr, "Cannot build Gaussian Distribution!\n");
		exit(1);
	}

	for (i = 0; i <= enc->maxval; i++) {
		for (j = 0; j <= (enc->maxval << 1); j++) {
			k = (i << 1) - j - 1;
			if (k < 0) k = -(k + 1);
			enc->econv[i][j] = k;
			// printf("[%d][%d]%d\n",i,j,k);
		}
	}
	enc->encval = enc->err;
	for (gr = 0; gr < enc->num_group; gr++) {
		var = enc->sigma[gr] * enc->sigma[gr];
		if (f_mmse) {
			a = 0;
			b = 1.0;
		} else {
			a = 0.5 * log(2 * M_PI * var) / log(2.0);
			b = 1.0 / (2.0 * log(2.0) * var);
		}
		enc->pmlist[gr] = pm = enc->pmodels[gr][gauss_index];
		for (k = 0; k < pm->size; k++) {
			c = (double)k * 0.5 + 0.25;
			pm->cost[k] = (float)(a + b * c * c);
		}
		pm->subcost[0] = 0.0;
	}
	for (k = 0; k <= enc->maxprd; k++) {
		enc->bconv[k] = 0;
		enc->fconv[k] = 0;
	}
	return;
}

void set_cost_rate(ENCODER *enc)	//分散ごとの確率モデルにおける符号量の配列を更新
{
	int gr, k, i, j, mask, shift, num_spm;
	double a, c;
	PMODEL *pm;

	enc->encval = enc->org;
	mask = (1 << enc->pm_accuracy) - 1;
	shift = enc->coef_precision - enc->pm_accuracy;
	for (k = 0; k <= enc->maxprd; k++) {
		i = (enc->maxprd - k + (1 << shift) / 2) >> shift;
		enc->fconv[k] = (i & mask);
		enc->bconv[k] = (i >> enc->pm_accuracy);
		 // printf("[%5d]i:%d | fconv:%d | bconv:%d\n", k, i, enc->fconv[k], enc->bconv[k]);
	}
	num_spm = 1 << enc->pm_accuracy;

	a = 1.0 / log(2.0);
	for (gr = 0; gr < enc->num_group; gr++) {
		for (i = 0; i < enc->num_pmodel; i++) {
			pm = enc->pmodels[gr][i];
			for (j = 0; j < num_spm; j++) {
				for (k = 0; k < pm->size; k++) {
					pm->cost[k] = (float)(-a * log(pm->freq[k]));
				}
				for (k = 0; k <= enc->maxval; k++) {
					c = pm->cumfreq[k + enc->maxval + 1] - pm->cumfreq[k];
					pm->subcost[k] = (float)(a * log(c));
				}
				pm++;
			}
		}
	}
}

void set_pmodel_mult_cost(MASK *mask,int size, int e)
{
	int p,base;
	PMODEL *pm;

	mask->freq = MIN_FREQ;
	mask->cumfreq = 0;
	for (p = 0; p < mask->num_peak; p++){
		base = mask->base[p];
		pm = mask->pm[p];
		mask->freq += (mask->weight[p] * (pm->freq[base + e] - MIN_FREQ)) >> W_SHIFT;
		mask->cumfreq += (mask->weight[p] * (pm->cumfreq[base + size] - pm->cumfreq[base])) >> W_SHIFT;
	}
	return;
}


void predict_region(ENCODER *enc, int tly, int tlx, int bry, int brx)	//予測値，予測誤差の再計算
{
	int x, y, k, l, cl, prd, org;
	int *coef_p, *nzc_p;
	int *prd_p;
	int *roff_p, **roff_pp;
	int *err_p, *org_p;
	char *class_p;

	for (y = tly; y < bry; y++) {
		class_p = &enc->class[y][tlx];
		org_p = &enc->org[y][tlx];
		roff_pp = &enc->roff[y][tlx];
		err_p = &enc->err[y][tlx];
		prd_p = &enc->prd[y][tlx];
		for (x = tlx; x < brx; x++) {
			cl = *class_p++;
			roff_p = *roff_pp++;
			coef_p = enc->predictor[cl];
			nzc_p = enc->nzconv[cl];
			prd = 0;
			if(enc->optimize_loop == 1){
				if(nzc_p[0] == -1){
					prd = exam_array[y][x][0];
				} else {
					for (k = 0; k < enc->num_nzcoef[cl]; k++) {
						l = nzc_p[k];
						prd += org_p[roff_p[l]] * (coef_p[l]);
					}
				}
			} else if(enc->optimize_loop==2){
				if(enc->num_nzcoef[cl] == -1){
					prd = exam_array[y][x][0];
				} else {
					for(k = 0; k < enc->num_nzcoef[cl]; k++){
						l = nzc_p[k];
						prd += org_p[roff_p[l]] * coef_p[l];
						// if(y == check_y && x== check_x && enc->function_number == 100)	printf("prd;%d|org;%d|roff;%d|coef;%d\n", prd, org_p[roff_p[l]], roff_p[l], coef_p[l]);
					}
				}
			}
			org = *org_p++;
			*prd_p++ = prd;
			prd = CLIP(0, enc->maxprd, prd);
			prd >>= (enc->coef_precision - 1);
			*err_p++ = enc->econv[org][prd];
			if(enc->function_number == F_NUM)	printf("%d|%d|err:|%d|org:|%d|prd:|%d\n", y, x, enc->err[y][x], org, prd);
		}
	}
}

void save_prediction_value(ENCODER *enc)
{
	int x, y, k, l, cl, prd;
	int *psave_p, *org_p;
	int **roff_pp, *roff_p, *coef_p, *nzc_p;

	for (y = 0; y < enc->height; y++) {
		roff_pp = enc->roff[y];
		for (x = 0; x < enc->width; x++) {
			roff_p = *roff_pp++;
			psave_p = enc->prd_class[y][x];
			org_p = &enc->org[y][x];
			for (cl = 0; cl < enc->num_class; cl++) {
				coef_p = enc->predictor[cl];
				nzc_p = enc->nzconv[cl];
				prd = 0;
				if(enc->num_nzcoef[cl] == -1){
					prd = exam_array[y][x][0];
				} else {
					for (k = 0; k < enc->num_nzcoef[cl]; k++) {
						l = nzc_p[k];
						prd += org_p[roff_p[l]] * (coef_p[l]);
					}
				}
				*psave_p++ = prd;	//クラス情報毎の予測値を一元保存
			}
		}
	}
}

int calc_uenc(ENCODER *enc, int y, int x)		//特徴量算出
{
	int u=0, k, *roff_p, *wt_p;
	#if CONTEXT_COST_MOUNT
		cost_t cost=0, *cost_p;
		cost_p = &enc->cost[y][x];
	#elif CONTEXT_ERROR
		int *err_p;
		err_p = &enc->err[y][x];
	#endif
	roff_p = enc->roff[y][x];
	wt_p = enc->ctx_weight;

	for (k =0; k < NUM_UPELS; k++) {
		#if CONTEXT_COST_MOUNT
			cost += cost_p[roff_p[k]] * wt_p[k];
		#elif CONTEXT_ERROR
			// u += err_p[*roff_p++] * (*wt_p++);
			u += err_p[roff_p[k]] * wt_p[k];
			#if CHECK_DEBUG
				if(y==check_y && x==check_x && enc->function_number == F_NUM)	printf("u: %d | err: %d(%3d) | wt_p: %d\n", u, err_p[roff_p[k]], roff_p[k], wt_p[k]);
			#endif
		#endif
	}
	#if CONTEXT_COST_MOUNT
		u = (int)(cost / NUM_UPELS);
	#elif CONTEXT_ERROR
		u >>= 6;
		// u >>= enc->coef_precision;
	#endif
	if (u > MAX_UPARA) u = MAX_UPARA;
	#if CHECK_DEBUG
		if(y==check_y && x==check_x && enc->function_number == F_NUM)	printf("u: %d\n", u);
	#endif
	return (u);
}

#if TEMPLATE_MATCHING_ON
void*** TemplateM (ENCODER *enc, char *outfile) {
	clock_t start=0, end=0;
	start = clock();
#if TEMPLATEM_LOG_OUTPUT
	int x , y , bx , by , i , j=0 , k , count , area1[AREA] , area_o[AREA] , *tm_array ,
		*roff_p , *org_p , x_size = X_SIZE , sum1 , sum_o, temp_x, temp_y, break_flag=0,
		**encval, temp_peak_num=0;
	double ave1=0 , ave_o , nas ;
	TM_Member tm[Y_SIZE * (X_SIZE * 2 + 1) + X_SIZE ], *tm_save;
#if ZNCC
	double dist1=0, dist_o=0, *area1_d=0, *area_o_d=0;
	area1_d = (double * )alloc_mem(AREA * sizeof(double));
	area_o_d = (double * )alloc_mem(AREA * sizeof(double));
	printf("ZNCC\tON\n");
#endif

#if MANHATTAN_SORT
	int *mcost_num, max_nas=0, before_nas_num=0;
	printf("MANHATTAN_SORT\tON\n");
#endif

/////////////////////////
///////メモリ確保////////
/////////////////////////

	tm_array = (int *)alloc_mem((Y_SIZE * (X_SIZE * 2 + 1) + X_SIZE) * 4 * sizeof(int)) ;
	encval = (int **)alloc_2d_array(enc->height+1, enc->width, sizeof(int));

	/*for(y=0; y<enc->height; y++){
		for(x=0; x<enc->width; x++){
			encval[y][x] = enc->org[y][x] << enc->coef_precision;
		}
	}*/
	init_2d_array(encval, enc->height, enc->width, 0);
	encval[enc->height][0] = (int)(enc->maxval + 1) << enc->coef_precision;

///////////////////////////
////////画像の走査/////////
///////////////////////////
printf("Calculating Template Matching\r");
for(y = 0 ; y < enc->height ; y++){
	for (x = 0; x < enc->width; x++){
		if(y==0 && x==0) {
			enc->temp_num[y][x] = 0;
			// encval[y][x] = enc->org[y][x] << enc->coef_precision;
			continue;
		}

		init_2d_array(encval, enc->height, enc->width, 0);
		for(i=0; i<enc->height; i++){
			if( i != y){
				for(j=0; j<enc->width; j++){
					encval[i][j] = enc->org[i][j] << enc->coef_precision;
				}
			} else if(i == y){
				for(j=0; j<x; j++){
					encval[i][j] = enc->org[i][j] << enc->coef_precision;
				}
			}
			if(i == y)break;
		}

		// bzero(&tm, sizeof(tm));
		memset(&tm, 0, sizeof(tm));

		roff_p = enc->roff[y][x];//init_ref_offsetが入っている．予測器の範囲指定と番号付け
		org_p = &encval[y][x];

		sum1 = 0;
		for(i=0;i < AREA; i++){//市街地距離AREA個分
			area1[i] = 0;
			area1[i] = org_p[roff_p[i]];
			sum1 += area1[i];
				#if CHECK_TM_DETAIL
					if(y==check_y && x==check_x)	printf("sum1: %d | area1[%d]:%d\n", sum1, i, area1[i]);
				#endif
		}
		ave1 = (double)sum1 / AREA;

	#if ZNCC
		dist1=0;
		for(i=0; i<AREA; i++){
			dist1 += ((double)area1[i] - ave1) * ((double)area1[i] - ave1);
			#if CHECK_TM_DETAIL
				if(y==check_y && x==check_x)	printf("dist1: %f | area1: %d | ave1: %f\n", dist1, area1[i], ave1);
			#endif
		}
		dist1 = sqrt(dist1);
		if(dist1 == 0) dist1 = 1;

		for(i=0; i<AREA; i++){
			area1_d[i] = ((double)area1[i] -ave1) / dist1;
			#if CHECK_TM_DETAIL
				if(y==check_y && x==check_x)	printf("area1_d: %f | area1: %d | ave1: %f | dist1: %f\n", area1_d[i], area1[i], ave1, dist1);
			#endif
		}
	#endif

////////////////////////////
//テンプレートマッチング///
///////////////////////////

		j = 0;
		break_flag=0;
		#if MANHATTAN_SORT
			max_nas = 0;
		#endif

		if(y == 0 || y == 1 || y == 2){
			x_size = 50;
		}else{
			x_size = X_SIZE;
		}

		for (by = y - Y_SIZE ; by <= y ; by++) {
			if((by < 0) || (by > enc->height))continue;
			for (bx = x - x_size ; bx <= x + x_size  ; bx++) {
				if((bx < 0) || (bx > enc->width))continue;
				if(by==y && bx >= x) break_flag=1;
				if ( break_flag )	break;

				roff_p = enc->roff[by][bx];
				org_p = &encval[by][bx];

				sum_o = 0;
				for(i=0;i < AREA; i++){//市街地距離AREA個分
					area_o[i] = 0;
					area_o[i] = org_p[roff_p[i]];
					sum_o += area_o[i];
					#if CHECK_TM_DETAIL
						if(y==check_y && x==check_x)	printf("sum_o: %d | area_o[%d]: %d\n",sum_o, i, area_o[i]);
					#endif
				}

				ave_o = (double)sum_o / AREA;

			#if ZNCC
				dist_o = 0;
				for(i=0; i<AREA; i++){
					dist_o += ((double)area_o[i] - ave_o) * ((double)area_o[i] - ave_o);
					#if CHECK_TM_DETAIL
						if(y==check_y && x==check_x)	printf("dist_o: %f | area_o: %d | ave_o: %f\n", dist_o, area_o[i], ave_o);
					#endif
				}
				dist_o = sqrt(dist_o);
				if(dist_o == 0)	dist_o = 1;

				for(i=0; i<AREA; i++){
					area_o_d[i] = ((double)area_o[i] - ave_o) / dist_o;
					#if CHECK_TM_DETAIL
						if(y==check_y && x==check_x)	printf("area_o_d: %f | area_o: %d | ave_o: %f | dist_o: %f\n", area_o_d[i], area_o[i], ave_o, dist_o);
					#endif
				}
			#endif

				nas = 0;

				for(i = 0; i < AREA ; i++){//マッチングコストの計算
					#if ZNCC
						// nas += fabs(area1_d[i] - area_o_d[i]);
						nas += (area1_d[i] - area_o_d[i]) * (area1_d[i] - area_o_d[i]);
						#if CHECK_TM_DETAIL
							if(y==check_y && x==check_x)	printf("nas: %f | area1: %f | area_o: %f\n", nas, area1_d[i], area_o_d[i]);
						#endif
					#elif ZSAD
						nas += fabs( ((double)area1[i] - ave1) - ((double)area_o[i] - ave_o) );
						#if CHECK_TM_DETAIL
							if(y==check_y && x==check_x)	printf("nas: %f | area1: %d | area_o: %d | ave1: %f | ave_o: %f\n", nas, area1[i], area_o[i], ave1, ave_o);
						#endif
					#else
						nas += (area1[i] - area_o[i] ) * (area1[i] - area_o[i]);
					#endif

				}

				tm[j].id = j;
				tm[j].by = by;
				tm[j].bx = bx;
				tm[j].ave_o = (int)ave_o;
				tm[j].sum = (int)(nas * NAS_ACCURACY);
				if(tm[j].sum < 0)	tm[j].sum = 0;

				#if ZNCC
					tm[j].s_devian = dist_o;
				#endif

				#if MANHATTAN_SORT
					tm[j].mhd = abs(x-bx) + abs(y-by);
					if(tm[j].sum > max_nas)	max_nas = tm[j].sum;
					// printf("(%3d,%3d)max_nas: %d | sum:%d\n", y, x, max_nas, tm[j].sum);
				#endif

				#if CHECK_TM
					if(y == check_y && x == check_x){
						printf("B[%3d](%3d,%3d) sum: %d | ave: %d", tm[j].id, tm[j].by, tm[j].bx, tm[j].sum, tm[j].ave_o);
						#if ZNCC
							printf(" | devian: %f", tm[j].s_devian);
						#endif
						printf("\n");
					}
				#endif

				j++;
			}//bx fin
			if( break_flag )break;
		}//by fin

/////////////////////////
///////ソートの実行//////
/////////////////////////

		enc->temp_num[y][x] = j;
		tm_save = (TM_Member *)alloc_mem(j * sizeof(TM_Member));
		for(i=0; i<j; i++){
			tm_save[i] = tm[i];
		}
		qsort(tm_save, j , sizeof(TM_Member), cmp);
		/*TM_Member temp;
		for (g = 0; g < j - 1; g++) {
			for (h = j - 1; h > g; h--) {
				if(tm[h - 1].sum > tm[h].sum) {	//前の要素の方が大きかったら
					temp = tm[h];		//交換する
					tm[h] = tm[h - 1];
					tm[h - 1] = temp;
				}
			}
		}*/

	#if MANHATTAN_SORT
		mcost_num = (int *)alloc_mem( (max_nas + 1) *sizeof(int));
		for(i=0; i <= max_nas; i++){
			mcost_num[i] = 0;
		}
		// mcost_num = memset(mcost_num, 0, sizeof(memset));
		before_nas_num=0;
		for( g = 0 ; g < j ; g++){
			mcost_num[tm[g].sum]++;
		}

		for( g = 0 ; g <= max_nas ; g++){
			if(mcost_num[g] == 0)continue;
			for(h=0; h<mcost_num[g]-1; h++){
				for(i=mcost_num[g] -1; i>h; i--){
					if(tm[i+before_nas_num-1].mhd > tm[i+before_nas_num].mhd){
						temp= tm[i+before_nas_num];
						tm[i+before_nas_num] = tm[i+before_nas_num-1];
						tm[i+before_nas_num-1] = temp;
					}
				}

			}
			before_nas_num += mcost_num[g];
		}

		free(mcost_num);
	#endif
		for(i=0; i<j; i++){
			tm[i] = tm_save[i];
		}
		free(tm_save);

		for(k = 0 ; k < j  ; k++){
			count = 0;
			tm_array[k * 4 + count] = 0;
			tm_array[k * 4 + count] = tm[k].id;
			count++;
			tm_array[k * 4 + count] = 0;
			tm_array[k * 4 + count] = tm[k].by;
			count++;
			tm_array[k * 4 + count] = 0;
			tm_array[k * 4 + count] = tm[k].bx;
			count++;
			tm_array[k * 4 + count] = 0;
			tm_array[k * 4 + count] = tm[k].sum;
			#if CHECK_TM
				if(y == check_y && x == check_x){
					printf("A[%3d](%3d,%3d) sum: %d | ave: %d", tm[k].id, tm[k].by, tm[k].bx, tm[k].sum, tm[k].ave_o);
					#if ZNCC
						printf(" | devian: %f\n", tm[k].s_devian);
					#endif
					printf("\n");
				}
			#endif
		}

		for(k = 0 ; k < MAX_DATA_SAVE_DOUBLE ; k++){
			tempm_array[y][x][k] = tm_array[k];
		}
		for(k = 0 ; k < MAX_DATA_SAVE ; k++){
			enc->array[y][x][k] = tm[k].ave_o;
		}

//マッチングコストが小さいものをTEMPLATE_CLASS_NUMの数だけ用意
		if(enc->temp_num[y][x] < TEMPLATE_CLASS_NUM){
			temp_peak_num = enc->temp_num[y][x];
		} else {
			temp_peak_num = TEMPLATE_CLASS_NUM;
		}
		for(i=0; i<temp_peak_num; i++){
			temp_y = tempm_array[y][x][i*4 + 1];
			temp_x = tempm_array[y][x][i*4 + 2];
			ave_o = enc->array[y][x][i];
			#if ZNCC
				dist_o = tm[i].s_devian;
			#endif

			if(y == 0 && x< 3){
				exam_array[y][x][i] = (enc->maxprd >> 1);
			} else {
				#if ZNCC
					exam_array[y][x][i] = (int)( ((double)encval[temp_y][temp_x] - ave_o) * dist1 / dist_o + ave1);
				#elif ZSAD
					exam_array[y][x][i] = (int)((double)encval[temp_y][temp_x] - ave_o + ave1);
				#else
					exam_array[y][x][i] = encval[temp_y][temp_x];
				#endif
				if(exam_array[y][x][i] < 0 || exam_array[y][x][i] > enc->maxprd)	exam_array[y][x][i] = (int)ave1;
			}
			#if CHECK_TM
				if(y==check_y && x==check_x)	printf("exam_array[%d]: %d[%3d] | (%3d,%3d) ave1: %f | ave_o: %f\n", i, exam_array[y][x][i], encval[temp_y][temp_x], temp_y, temp_x, ave1, ave_o);
			#endif
		}
		#if CHECK_TM
			if(y==check_y && x==check_x)	printf("(%3d,%3d)org: %d\n", y, x, encval[y][x]);
		#endif
		// encval[y][x] = enc->org[y][x] << enc->coef_precision;
	}//x fin
}//y fin
TemplateM_Log_Output(enc, outfile, tempm_array, exam_array);
free(tm_array);
free(encval);
#else
printf("Restoring Template Matching\r");
TemplateM_Log_Input(enc, outfile, tempm_array, exam_array);
#endif
end = clock();
printf("Calculating Template Matching Fin[%f sec]\n", (float)(end-start)/CLOCKS_PER_SEC);

	return(0);
}

double continuous_GGF(ENCODER *enc, double e,int w_gr){
	int lngamma(double), cn=WEIGHT_CN, num_pmodel=enc->num_pmodel;
	double sigma,delta_c,shape,eta,p;
	double accuracy = 1 / (double)NAS_ACCURACY;

	sigma = enc->sigma[w_gr];
	delta_c = 3.2 / (double)num_pmodel;//cnは0.2ずつ動くので. num_pmodel is 16 . delta_c = 0.2
	shape = delta_c * (double)(cn);//cnを実際の数値に変換
	eta = exp(0.5*(lgamma(3.0/shape)-lgamma(1.0/shape))) / sigma;//一般化ガウス関数.ηのみ

	if(e <= accuracy){
		p = 1.0;
	}else{
		p = exp(-pow(eta * (e), shape));
	}

	return(p);
}

int temp_mask_parameter(ENCODER *enc, int y, int x, int u, int peak, int cl, int weight_all, int w_gr)
{
	int i, m_gr, m_prd, m_frac, template_peak = enc->temp_peak_num;
	double weight[template_peak], sum_weight = 0, weight_coef=0;
	if(template_peak > enc->temp_num[y][x])	template_peak = enc->temp_num[y][x];
	if(template_peak == 0){
		mask->class[peak] = cl;
		mask->weight[peak] = weight_all;
		m_gr = enc->uquant[cl][u];
		m_prd = enc->maxprd > 1;
		m_prd = CLIP(0, enc->maxprd, m_prd);
		mask->base[peak] = enc->bconv[m_prd];
		m_frac = enc->fconv[m_prd];
		mask->pm[peak] = enc->pmlist[m_gr] + m_frac;
		peak++;
		return(peak);
	}
	// if(y==check_y && x==check_x && enc->function_number == F_NUM)	printf("w_gr: %d\n", w_gr);
	for(i=0; i<template_peak; i++){
		weight[i] = continuous_GGF(enc, (double)(tempm_array[y][x][i*4+3] >> enc->coef_precision) / NAS_ACCURACY, w_gr);
		// if(y== check_y && x==check_x)	printf("weight[%d]: %.20f\n", i, weight[i]);
		sum_weight += weight[i];
		// if(y == check_y && x== check_x && enc->function_number == F_NUM)	printf("sum: %.20f | weight: %.20f\n", sum_weight, weight[i]);
	}
	if(sum_weight == 0){
		weight_coef = (double)weight_all;
	} else {
		weight_coef = (double)weight_all / sum_weight;
	}
	// if(y==check_y && x==check_x && enc->function_number == F_NUM)	printf("weight_coef: %.20f | sum_weight: %.20f\n", weight_coef, sum_weight);

	for(i=0; i<template_peak; i++){
		mask->class[peak] = cl;
		// if(y==check_y && x==check_x)	printf("weight[%d]:%f\n", i, weight[i] * weight_coef);
		mask->weight[peak] = (int)(weight[i] * weight_coef);
		if(mask->weight[peak] == 0)	continue;
		m_gr = enc->uquant[cl][u];
		m_prd = exam_array[y][x][i];
		// m_prd = enc->prd_class[y][x][cl];
		m_prd = CLIP(0, enc->maxprd, m_prd);
		#if CHECK_DEBUG
			if( y == check_y && x == check_x && enc->function_number == F_NUM)	printf("[temp_mask_parameter]\tm_prd[%d]: %d[%2d] | weight: %d | u: %d | gr: %d\n", peak, m_prd, cl, mask->weight[peak], u, m_gr);
		#endif

		mask->base[peak] = enc->bconv[m_prd];
		m_frac = enc->fconv[m_prd];
		mask->pm[peak] = enc->pmlist[m_gr] + m_frac;
		peak++;
	}

	return(peak);
}

int encode_w_gr_threshold(FILE *fp, ENCODER *enc, int flag)
{
	int gr, i, k, m, min_m, bits;
	cost_t cost, min_cost;
	PMODEL *pm;
	double p;
	pm = &enc->spm;

#if (!OPT_SIDEINFO)
	if (fp == NULL) return(0);
#endif
	/* Arithmetic */
	min_cost = INT_MAX;
	for (m = min_m = 0; m < 16; m++) {
		set_spmodel(pm, MAX_UPARA + 2, m);
		cost = 0.0;
		// for (cl = 0; cl < enc->num_class; cl++) {
			k = 0;
			for (gr = 1; gr < enc->num_group; gr++) {
				// i = enc->th[cl][gr - 1] - k;
				i = enc->w_gr[gr-1] - k;
				p = (double)pm->freq[i] / (pm->cumfreq[pm->size - k]);
				cost += -log(p);
				k += i;
				if (k > MAX_UPARA) break;
			}
		// }
		cost /= log(2.0);
		if (cost < min_cost) {
			min_cost = cost;
			min_m = m;
		}
	}
	set_spmodel(pm, MAX_UPARA + 2, min_m);
	p = log(pm->cumfreq[MAX_UPARA + 2]);
	if (fp == NULL) {
		/*if (flag == 1){
			for (i = 0; i < MAX_UPARA + 2; i++) {
				enc->th_cost[i] = (p - log(pm->freq[i])) / log(2.0);
			}
		}*/
		bits = (int)min_cost;
	} else {
		rc_encode(fp, enc->rc, min_m, 1, 16);
		// for (cl = 0; cl < enc->num_class; cl++) {
			k = 0;
			for (gr = 1; gr < enc->num_group; gr++) {
				// i = enc->th[cl][gr - 1] - k;
				i = enc->w_gr[gr-1] - k;
				rc_encode(fp, enc->rc, pm->cumfreq[i],  pm->freq[i],
					pm->cumfreq[pm->size - k]);
				k += i;
				if (k > MAX_UPARA) break;
			}

		/*if (enc->num_pmodel > 1) {
			for (gr = 0; gr < enc->num_group; gr++) {
				pm = enc->pmlist[gr];
				rc_encode(fp, enc->rc, pm->id, 1, enc->num_pmodel);
			}
		}*/
		bits = (int)enc->rc->code;
		enc->rc->code = 0;
	}
	return (bits);
}

#endif

void set_mask_parameter(ENCODER *enc,int y, int x, int u)
{
	int ty, tx, cl, i, peak, sample;
	int m_gr,m_prd,m_frac;
	int count_cl[enc->num_class];

	sample = win_sample[(int)enc->mask[y][x]];	//ウィンドウ内の画素数

	for (cl = 0; cl < enc->num_class; cl++) {
		count_cl[cl]  = 0;
	}

	for (i = 0; i < sample; i++){
		ty = y + mask_y[i];
		tx = x + mask_x[i];
		if(ty < 0) ty = 0;
		else if(ty >= enc->height) ty = enc->height - 1;
		if(tx < 0) tx = 0;
		else if(tx >= enc->width) tx = enc->width - 1;
		cl = enc->class[ty][tx];
		count_cl[cl]++;	//マスクにかかる領域にあるクラス毎の画素数
	}

	for(cl = peak = 0; cl < enc->num_class; cl++){
		if (count_cl[cl]!=0){
		#if TEMPLATE_MATCHING_ON
			if(cl == enc->temp_cl){
				for(m_gr=0; m_gr < enc->num_group; m_gr++){
					if(u < enc->w_gr[m_gr])	break;
				}
				if(m_gr >= enc->num_group)	m_gr = enc->num_group - 1;
				peak = temp_mask_parameter(enc, y, x, u, peak, cl, (count_cl[cl] << W_SHIFT) / sample, m_gr);
			} else {
		#endif
				mask->class[peak] = cl;		//マスクにかかる領域毎のクラスを保存
				mask->weight[peak] = ( (count_cl[cl] << W_SHIFT) / sample);	//	各ピーク毎の重み
				m_gr = enc->uquant[cl][u];
				m_prd = enc->prd_class[y][x][cl];
				m_prd = CLIP(0, enc->maxprd, m_prd);
				#if CHECK_DEBUG
					if( y == check_y && x == check_x && enc->function_number== F_NUM)	printf("[set_mask_parameter]\tm_prd[%d]: %d[%2d] | weight: %d | u: %d | gr: %d\n", 	peak, m_prd, cl, mask->weight[peak], u, m_gr);
				#endif

				mask->base[peak] = enc->bconv[m_prd];
				m_frac = enc->fconv[m_prd];
				mask->pm[peak] = enc->pmlist[m_gr] + m_frac;
				peak++;
		#if TEMPLATE_MATCHING_ON
			}
		#endif
		}
	}

	mask->num_peak = peak;	//ピークの数
}

int set_mask_parameter_optimize(ENCODER *enc,int y, int x, int u, int r_cl)
{
	int cl, peak;
	int m_gr,m_prd,m_frac,r_peak;
	r_peak = -1;

	for(cl = peak = 0; cl < enc->num_class; cl++){
		if (enc->weight[y][x][cl] != 0){
		#if TEMPLATE_MATCHING_ON
			if(cl == enc->temp_cl){
				for(m_gr=0; m_gr < enc->num_group; m_gr++){
					if(u < enc->w_gr[m_gr])	break;
				}
				peak = temp_mask_parameter(enc, y, x, u, peak, cl, enc->weight[y][x][cl], m_gr);
				if(cl == r_cl)	r_peak = peak - (enc->temp_peak_num - 1);	//マッチングコストが一番小さい事例のピーク番号
			} else {
		#endif
				mask->class[peak] = cl;
				mask->weight[peak] = enc->weight[y][x][cl];
				m_prd = enc->prd_class[y][x][cl];
				m_prd = CLIP(0, enc->maxprd, m_prd);
				mask->base[peak] = enc->bconv[m_prd];
				m_frac = enc->fconv[m_prd];
				m_gr = enc->uquant[cl][u];
				#if CHECK_DEBUG
					if( y == check_y && x == check_x && enc->function_number == F_NUM)	printf("[set_mask_parameter_opt] m_prd[%d]: %d[%2d] | weight: %d | gr: %d\n", 	peak, m_prd, cl, mask->weight[peak], m_gr);
				#endif
				if (cl == r_cl)	 r_peak = peak;	//当該ブロックのピーク番号
				mask->pm[peak] = enc->pmlist[m_gr] + m_frac;
				peak++;
		#if TEMPLATE_MATCHING_ON
			}
		#endif
		}
	}

	mask->num_peak = peak;	//ピークの数
	return(r_peak);
}

cost_t calc_cost2(ENCODER *enc, int tly, int tlx, int bry, int brx)
{
	cost_t cost;
	int x, y, u, cl, gr, e, base;
	int *upara_p, *encval_p;
	char *class_p, *group_p;
	double a;
	PMODEL *pm;

	a = 1.0 / log(2.0);

	if (bry > enc->height) bry = enc->height;
	if (tlx < 0) tlx = 0;
	if (tly < 0) tly = 0;
	if (brx > enc->width) brx = enc->width;
	cost = 0;
	for (y = tly; y < bry; y++) {
		class_p = &enc->class[y][tlx];
		group_p = &enc->group[y][tlx];
		upara_p = &enc->upara[y][tlx];
		encval_p = &enc->encval[y][tlx];
		for (x = tlx; x < brx; x++) {
			*upara_p++ = u = calc_uenc(enc, y, x);
			set_mask_parameter(enc,y,x,u);
			cl = *class_p++;
			*group_p++ = gr = enc->uquant[cl][u];
			e = *encval_p++;
			if (mask->num_peak == 1){
				base = mask->base[0];
				pm = mask->pm[0];
				enc->cost[y][x] = pm->cost[base + e] + pm->subcost[base];
			}else{
				set_pmodel_mult_cost(mask,enc->maxval+1,e);
				enc->cost[y][x] = a * (log(mask->cumfreq)-log(mask->freq)) ;
			}
			cost += enc->cost[y][x];
			if(y==check_y && x==check_x && enc->function_number == F_NUM)	printf("cost(%3d,%3d): %f\n", y, x, enc->cost[y][x]);
			// if(enc->function_number== F_NUM) printf("(%3d,%3d) cost: %d | cl: %d\n", y, x, (int)cost, cl );
		}
	}
	return (cost);
}


cost_t calc_cost(ENCODER *enc, int tly, int tlx, int bry, int brx)		//コストを算出
{
	cost_t cost;
	int x, y, u, cl, gr, prd, e, base, frac;
	int *upara_p, *prd_p, *encval_p;
	char *class_p, *group_p;
	PMODEL *pm;

# if OPTIMIZE_MASK_LOOP
	if(enc->optimize_loop == 2){
		return (calc_cost2(enc,tly,tlx,bry,brx));
	}
#endif

	if (bry > enc->height) bry = enc->height;
	if (brx > enc->width) brx = enc->width;
	if (tlx < 0) tlx = 0;
	if (tly < 0) tly = 0;
	cost = 0;
	for (y = tly; y < bry; y++) {
		class_p = &enc->class[y][tlx];
		group_p = &enc->group[y][tlx];
		upara_p = &enc->upara[y][tlx];
		encval_p = &enc->encval[y][tlx];
		prd_p = &enc->prd[y][tlx];
		for (x = tlx; x < brx; x++) {
			cl = *class_p++;
			*upara_p++ = u = calc_uenc(enc, y, x);
			*group_p++ = gr = enc->uquant[cl][u];
			e = *encval_p++;
			prd = *prd_p++;
			prd = CLIP(0, enc->maxprd, prd);
			base = enc->bconv[prd];
			frac = enc->fconv[prd];
			pm = enc->pmlist[gr] + frac;
			enc->cost[y][x] = pm->cost[base + e] + pm->subcost[base];
			if(y==check_y && x==check_x && enc->function_number == F_NUM)	printf("cost(%3d,%3d): %f\n", y, x, enc->cost[y][x]);
			cost += enc->cost[y][x];
		}
	}
	if(cost < 0) cost = INT_MAX;
	return (cost);
}

cost_t design_predictor(ENCODER *enc, int f_mmse)
{
	double **mat, *weight, w, e, d, pivot;
	int x, y, i, j, k, cl, gr, pivpos, *index, *roff_p, *org_p, *nzc_p;
	enc->function_number = 1;
	mat = (double **)alloc_2d_array(enc->prd_order, enc->prd_order + 1, sizeof(double));
	index = (int *)alloc_mem(sizeof(int) * enc->prd_order);
	weight = (double *)alloc_mem(sizeof(double) * enc->num_group);

	for (gr = 0; gr < enc->num_group; gr++) {
		if (f_mmse) {
			weight[gr] = 1.0;
		} else {
			weight[gr] = 1.0 / (enc->sigma[gr] * enc->sigma[gr]);
		}
	}
	for (cl = 0; cl < enc->num_class; cl++) {
		nzc_p = enc->nzconv[cl];
#if TEMPLATE_MATCHING_ON
		if(cl == 0){
			nzc_p[0] = -1;	//nzcの頭にフラグを入れる
			for(i=0; i < enc->prd_order; i++){
				enc->predictor[cl][i] = 0;
			}
			enc->predictor[cl][0] = TEMPLATE_FLAG;
			continue;
		}
#endif
		for (i = 0; i < enc->prd_order; i++) {
			for (j = 0; j <= enc->prd_order; j++) {
				mat[i][j] = 0.0;
			}
		}
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				if (enc->class[y][x] != cl) {
					x += BASE_BSIZE - 1;
					continue;
				}
				gr = enc->group[y][x];
				roff_p = enc->roff[y][x];
				org_p = &enc->org[y][x];
				for (i = 0; i < enc->prd_order; i++) {
					w = weight[gr] * org_p[roff_p[nzc_p[i]]];
					for (j = i; j < enc->prd_order; j++) {
						mat[i][j] += w * org_p[roff_p[nzc_p[j]]];
					}
					mat[i][enc->prd_order] += w * org_p[0];
				}
			}
		}
		// 掃き出し法
		for (i = 0; i < enc->prd_order; i++) {
			index[i] = i;
			for (j = 0; j < i; j++) {
				mat[i][j] = mat[j][i];
			}
		}
		for (i = 0; i < enc->prd_order; i++) {
			pivpos = i;
			pivot = fabs(mat[index[i]][i]);
			for (k = i + 1; k < enc->prd_order; k++) {
				if (fabs(mat[index[k]][i]) > pivot) {
					pivot = fabs(mat[index[k]][i]);
					pivpos = k;
				}
			}
			k = index[i];
			index[i] = index[pivpos];
			index[pivpos] = k;
			if (pivot > 1E-10) {
				d = mat[index[i]][i];
				for (j = i; j <= enc->prd_order; j++) {
					mat[index[i]][j] /= d;
				}
				for (k = 0; k < enc->prd_order; k++) {
					if (k == i) continue;
					d = mat[index[k]][i];
					for (j = i; j <= enc->prd_order; j++) {
						mat[index[k]][j] -= d * mat[index[i]][j];
					}
				}
			}
		}
		w = (1 << enc->coef_precision);		// 予測係数の精度
		e = 0.0;
		for (i = 0; i < enc->prd_order; i++) {
			if (fabs(mat[index[i]][i]) > 1E-10) {
				d = mat[index[i]][enc->prd_order] * w;
			} else {
				d = 0.0;
			}
			k = (int)d;
			if (k > d) k--;
			if (k < -enc->max_coef) {
				d = k = -enc->max_coef;
			} else if (k > enc->max_coef) {
				d = k = enc->max_coef;
			}
			enc->predictor[cl][nzc_p[i]] = k;
			d -= k;
			e += d;
			mat[index[i]][enc->prd_order] = d;
		}
		/* minimize mean rounding errors */
		k = (int)(e + 0.5);
		for (;k > 0; k--) {
			d = 0;
			for (j = i = 0; i < enc->prd_order; i++) {
				if (mat[index[i]][enc->prd_order] > d) {
					d = mat[index[i]][enc->prd_order];
					j = i;
				}
			}
			if (enc->predictor[cl][nzc_p[j]] < enc->max_coef) enc->predictor[cl][nzc_p[j]]++;
			mat[index[j]][enc->prd_order] = 0;
		}
	}//cl loop fin
	free(weight);
	free(index);
	free(mat);

	predict_region(enc, 0, 0, enc->height, enc->width);
	return (calc_cost(enc, 0, 0, enc->height, enc->width));
}


cost_t optimize_group(ENCODER *enc)
{
	cost_t cost, min_cost, **cbuf, *dpcost, *cbuf_p, *thc_p;
	int x, y, th1, th0, k, u, cl, gr, prd, e, base, frac;
	int **trellis;
	PMODEL *pm, **pm_p;
	enc->function_number = 2;
	trellis = (int **)alloc_2d_array(enc->num_group, MAX_UPARA + 2,
		sizeof(int));
	dpcost = (cost_t *)alloc_mem((MAX_UPARA + 2) * sizeof(cost_t));
	cbuf = (cost_t **)alloc_2d_array(enc->num_group, MAX_UPARA + 2,
		sizeof(cost_t));
	thc_p = enc->th_cost;
	for (k = 0; k < MAX_UPARA + 2; k++) trellis[0][k] = 0;
	/* Dynamic programming */
	for (cl = 0; cl < enc->num_class; cl++) {
		//if (enc->cl_hist[cl] == 0) continue;
		for (gr = 0; gr < enc->num_group; gr++) {
			cbuf_p = cbuf[gr];
			for (u = 0; u < MAX_UPARA + 2; u++) {
				cbuf_p[u] = 0;
			}
		}
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				if (enc->class[y][x] == cl) {
					u = enc->upara[y][x] + 1;
					e = enc->encval[y][x];
					prd = enc->prd[y][x];
					prd = CLIP(0, enc->maxprd, prd);
					base = enc->bconv[prd];
					frac = enc->fconv[prd];
					pm_p = enc->pmlist;
					for (gr = 0; gr < enc->num_group; gr++) {
						pm = (*pm_p++) + frac;
						cbuf[gr][u] += pm->cost[base + e] + pm->subcost[base];
					}
				}
			}
		}
		for (gr = 0; gr < enc->num_group; gr++) {
			cbuf_p = cbuf[gr];
			for (u = 1; u < MAX_UPARA + 2; u++) {
				cbuf_p[u] += cbuf_p[u - 1];
			}
		}
		cbuf_p = cbuf[0];
		for (u = 0; u < MAX_UPARA + 2; u++) {
			dpcost[u] = cbuf_p[u] + thc_p[u];
		}
		for (gr = 1; gr < enc->num_group - 1; gr++) {
			cbuf_p = cbuf[gr];
			/* minimize (cbuf_p[th1] - cbuf_p[th0] + dpcost[th0]) */
			for (th1 = MAX_UPARA + 1; th1 >= 0; th1--) {
				th0 = th1;
				min_cost = dpcost[th1] - cbuf_p[th1] + thc_p[0];
				for (k = 0; k < th1; k++) {
					cost = dpcost[k] - cbuf_p[k] + thc_p[th1 - k];
					if (cost < min_cost) {
						min_cost = cost;
						th0 = k;
					}
				}
				dpcost[th1] = min_cost + cbuf_p[th1];
				trellis[gr][th1] = th0;
			}
		}

		cbuf_p = cbuf[gr];
		/* minimize (cbuf_p[th1] - cbuf_p[th0] + dpcost[th0]) for last group */
		th1 = MAX_UPARA + 1;
		th0 = th1;
		min_cost = dpcost[th1] - cbuf_p[th1];
		for (k = 0; k < th1; k++) {
			cost = dpcost[k] - cbuf_p[k];
			if (cost < min_cost) {
				min_cost = cost;
				th0 = k;
			}
		}
		trellis[gr][th1] = th0;

		for (gr = enc->num_group - 1; gr > 0; gr--) {
			th1 = trellis[gr][th1];
			enc->th[cl][gr - 1] = th1;
		}
	}
	/* set context quantizer */
	for (cl = 0; cl < enc->num_class; cl++) {
		u = 0;
		for (gr = 0; gr < enc->num_group; gr++) {
			for (; u < enc->th[cl][gr]; u++) {
				enc->uquant[cl][u] = gr;
			}
		}
	}
	/* renew groups */
	cost = 0;
	pm_p = enc->pmlist;
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			cl = enc->class[y][x];
			u = enc->upara[y][x];
			enc->group[y][x] = gr = enc->uquant[cl][u];
			e = enc->encval[y][x];
			prd = enc->prd[y][x];
			prd = CLIP(0, enc->maxprd, prd);
			base = enc->bconv[prd];
			pm = pm_p[gr] + enc->fconv[prd];
			cost += pm->cost[base + e] + pm->subcost[base];
		}
	}
	/* optimize probability models */
	if (enc->optimize_loop > 1 && enc->num_pmodel > 1) {
		if (enc->num_pmodel > MAX_UPARA + 2) {
			free(cbuf);
			cbuf = (cost_t **)alloc_2d_array(enc->num_group, enc->num_pmodel,
				sizeof(cost_t));
		}
		for (gr = 0; gr < enc->num_group; gr++) {
			for (k = 0; k < enc->num_pmodel; k++) {
				cbuf[gr][k] = 0;
			}
		}
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				gr = enc->group[y][x];
				e = enc->encval[y][x];
				prd = enc->prd[y][x];
				prd = CLIP(0, enc->maxprd, prd);
				base = enc->bconv[prd];
				frac = enc->fconv[prd];
				for (k = 0; k < enc->num_pmodel; k++) {
					pm = enc->pmodels[gr][k] + frac;
					cbuf[gr][k] += pm->cost[base + e] + pm->subcost[base];
				}
			}
		}
		for (gr = 0; gr < enc->num_group; gr++) {
			pm = enc->pmodels[gr][0];
			cost = cbuf[gr][0];
			for (k = 1; k < enc->num_pmodel; k++) {
				if (cost > cbuf[gr][k]) {
					cost = cbuf[gr][k];
					pm = enc->pmodels[gr][k];
				}
			}
			pm_p[gr] = pm;
		}
		cost = 0.0;
		for (gr = 0; gr < enc->num_group; gr++) {
			cost += cbuf[gr][pm_p[gr]->id];
		}
	}
	free(cbuf);
	free(dpcost);
	free(trellis);
	return (cost);
}

void set_prdbuf(ENCODER *enc, int **prdbuf, int **errbuf,
				int tly, int tlx, int bufsize)
{
	int x, y, brx, bry, cl, k, l, prd, *prdbuf_p, *errbuf_p, *coef_p;
	int buf_ptr, org, *org_p, *roff_p, *nzc_p;

	brx = (tlx + bufsize < enc->width) ? (tlx + bufsize) : enc->width;
	bry = (tly + bufsize < enc->height) ? (tly + bufsize) : enc->height;
	for (cl = 0; cl < enc->num_class; cl++) {
		buf_ptr = bufsize * (tly % bufsize) + tlx % bufsize;
		nzc_p = enc->nzconv[cl];
		for (y = tly; y < bry; y++) {
			prdbuf_p = &prdbuf[cl][buf_ptr];
			errbuf_p = &errbuf[cl][buf_ptr];
			buf_ptr += bufsize;
			org_p = &enc->org[y][tlx];
			for (x = tlx; x < brx; x++) {
				if (cl == enc->class[y][x]) {
					*prdbuf_p++ = enc->prd[y][x];
					*errbuf_p++ = enc->err[y][x];
					org_p++;
				} else {
					coef_p = enc->predictor[cl];
					roff_p = enc->roff[y][x];
					prd = 0;
					for (k = 0; k < enc->num_nzcoef[cl]; k++) {
						if(nzc_p[k] == TEMPLATE_FLAG){
							prd = exam_array[y][x][0];
							break;
						} else {
							l = nzc_p[k];
							prd += org_p[roff_p[l]] * (coef_p[l]);
						}
					}
					org = *org_p++;
					*prdbuf_p++ = prd;
					prd = CLIP(0, enc->maxprd, prd);
					prd >>= (enc->coef_precision - 1);
					*errbuf_p++ = enc->econv[org][prd];
				}
			}
		}
	}
}

cost_t rap_area_cost(ENCODER *enc, int tly, int tlx, int bry, int brx)
{
	cost_t cost;
	int d;

	cost = 0;
			if (tly != 0){
				if (tlx != 0){//º¸¾å
					d = enc->mask[tly-1][tlx-1];
					cost += calc_cost2(enc, tly-d, tlx-d, tly, tlx);
				}

				d = enc->mask[tly-1][tlx];//¾å
				cost += calc_cost2(enc, tly-d, tlx, tly, brx);

				if (brx != enc->width){//±¦¾å
					d = enc->mask[tly-1][brx];
					cost += calc_cost2(enc, tly-d, brx, tly, brx+d);
				}
			}

			if (tlx != 0){//º¸
				d = enc->mask[tly][tlx-1];
				cost += calc_cost2(enc, tly, tlx-d, bry, tlx);
#if 1
				if (bry != enc->height){//º¸²¼
					d = enc->mask[bry][tlx-1];
					cost += calc_cost2(enc, bry, tlx-d, bry+d, tlx);
				}
#endif
			}

#if 1
			if (brx != enc->width){//±¦
				d = enc->mask[tly][brx];
				cost += calc_cost2(enc, tly, brx, bry, brx+d);
			}
#endif

#if 1
			if (bry != enc->height){//²¼
				d = enc->mask[bry][tlx];
				cost += calc_cost2(enc, bry, tlx, bry+d, brx);
#if 1
				if (brx != enc->width){//±¦²¼
					d = enc->mask[bry][brx];
					cost += calc_cost2(enc, bry, brx, bry+d, brx+d);
				}
#endif
			}
#endif
	return(cost);
}

int find_class(ENCODER *enc, int **prdbuf, int **errbuf,
			   int tly, int tlx, int bry, int brx, int bufsize, cost_t *err_cost)
{
	cost_t cost, min_cost;
	int x, y, bufptr, cl, min_cl;
	char *class_p;
	int *prd_p, *prdbuf_p, *err_p, *errbuf_p;
	min_cost = 1E8;
	min_cl = 0;

	for	(cl = 0; cl < enc->num_class; cl++) {
		bufptr = bufsize * (tly % bufsize) + tlx % bufsize;
		for (y = tly; y < bry; y++) {
			class_p = &enc->class[y][tlx];
			prd_p = &enc->prd[y][tlx];
			prdbuf_p = &prdbuf[cl][bufptr];
			err_p = &enc->err[y][tlx];
			errbuf_p = &errbuf[cl][bufptr];
			bufptr += bufsize;
			for (x = tlx; x < brx; x++) {
				*class_p++ = cl;
				*prd_p++ = *prdbuf_p++;
				*err_p++ = *errbuf_p++;
			}
		}
		cost = calc_cost(enc, tly, tlx, bry, brx);
		err_cost[cl] = cost;
		if( enc->optimize_loop == 2) {
#if OPTIMIZE_MASK_LOOP
			cost += rap_area_cost(enc, tly, tlx, bry, brx);
		  err_cost[cl] = cost;
#endif
			cost += enc->class_cost[enc->mtfbuf[cl]];
		}
		if (cost < min_cost) {
			min_cost = cost;
			min_cl = cl;
		}
	}
	bufptr = bufsize * (tly % bufsize) + tlx % bufsize;
	for (y = tly; y < bry; y++) {
		class_p = &enc->class[y][tlx];
		prd_p = &enc->prd[y][tlx];
		prdbuf_p = &prdbuf[min_cl][bufptr];
		err_p = &enc->err[y][tlx];
		errbuf_p = &errbuf[min_cl][bufptr];
		bufptr += bufsize;
		for (x = tlx; x < brx; x++) {
			*class_p++ = min_cl;
			*prd_p++ = *prdbuf_p++;
			*err_p++ = *errbuf_p++;
		}
	}
	return (min_cl);
}

cost_t vbs_class(ENCODER *enc, int **prdbuf, int **errbuf, int tly, int tlx,
				 int blksize, int width, int level, int *blk)
{
	int y, x, k, bry, brx, cl, bufsize, bufptr, ctx, s_blk;
	int mtf_save[MAX_CLASS];
	char **qtmap;
	cost_t cost1, cost2, qtcost, *err_cost;
	char *class_p;
	int *prd_p, *prdbuf_p, *err_p, *errbuf_p;

	s_blk = *blk;
	if (enc->quadtree_depth >= 0 && enc->optimize_loop > 1) {
		bufsize = MAX_BSIZE;
	} else {
		bufsize = BASE_BSIZE;
	}
	brx = (tlx + blksize < enc->width) ? (tlx + blksize) : enc->width;
	bry = (tly + blksize < enc->height) ? (tly + blksize) : enc->height;
	if (tlx >= brx || tly >= bry) return (0);
	for (k = 0; k < enc->num_class; k++) {
		mtf_save[k] = enc->mtfbuf[k];
	}
	mtf_classlabel(enc->class, enc->mtfbuf, tly, tlx,
		blksize, width, enc->num_class);
	err_cost = (cost_t *)alloc_mem(enc->num_class * sizeof(cost_t));
	cl = find_class(enc, prdbuf, errbuf, tly, tlx, bry, brx, bufsize, err_cost);
	qtcost = enc->class_cost[enc->mtfbuf[cl]];
	if (level > 0) {
		/* context for quad-tree flag */
		ctx = 0;
		qtmap = enc->qtmap[level - 1];
		y = (tly / MIN_BSIZE) >> level;
		x = (tlx / MIN_BSIZE) >> level;
		if (y > 0) {
			if (qtmap[y - 1][x] == 1) ctx++;
			if (brx < width && qtmap[y - 1][x + 1] == 1) ctx++;
		}
		if (x > 0 && qtmap[y][x - 1] == 1) ctx++;
		ctx = ((level - 1) * 4 + ctx) << 1;
		/* Quad-tree partitioning */
		cost1 = calc_cost(enc, tly, tlx, bry, brx)
			+ enc->class_cost[enc->mtfbuf[cl]] + enc->qtflag_cost[ctx];
#if OPTIMIZE_MASK_LOOP
		cost1 +=  rap_area_cost(enc, tly, tlx, bry, brx);
#endif
		blksize >>= 1;
		for (k = 0; k < enc->num_class; k++) {
			enc->mtfbuf[k] = mtf_save[k];
		}
		qtcost = enc->qtflag_cost[ctx + 1];
		qtcost += vbs_class(enc, prdbuf, errbuf, tly, tlx,
			blksize, width, level - 1, blk);
		qtcost += vbs_class(enc, prdbuf, errbuf, tly, tlx+blksize,
			blksize, width, level - 1, blk);
		qtcost += vbs_class(enc, prdbuf, errbuf, tly+blksize, tlx,
			blksize, width, level - 1, blk);
		qtcost += vbs_class(enc, prdbuf, errbuf, tly+blksize, tlx+blksize,
			blksize, brx, level - 1, blk);
		cost2 = calc_cost(enc, tly, tlx, bry, brx) + qtcost;
#if OPTIMIZE_MASK_LOOP
		cost2 += rap_area_cost(enc, tly, tlx, bry, brx);
#endif
		if (cost1 < cost2) {
			blksize <<= 1;
			for (k = 0; k < enc->num_class; k++) {
				enc->mtfbuf[k] = mtf_save[k];
			}
			mtf_classlabel(enc->class, enc->mtfbuf, tly, tlx,
				blksize, width, enc->num_class);
			qtcost = enc->class_cost[enc->mtfbuf[cl]]
			+ enc->qtflag_cost[ctx];
			bufptr = bufsize * (tly % bufsize) + tlx % bufsize;
			for (y = tly; y < bry; y++) {
				class_p = &enc->class[y][tlx];
				prd_p = &enc->prd[y][tlx];
				prdbuf_p = &prdbuf[cl][bufptr];
				err_p = &enc->err[y][tlx];
				errbuf_p = &errbuf[cl][bufptr];
				bufptr += bufsize;
				for (x = tlx; x < brx; x++) {
					*class_p++ = cl;
					*prd_p++ = *prdbuf_p++;
					*err_p++ = *errbuf_p++;
				}
			}
			tly = (tly / MIN_BSIZE) >> level;
			tlx = (tlx / MIN_BSIZE) >> level;
			bry = tly + 1;
			brx = tlx + 1;
			for (; level > 0; level--) {
				qtmap = enc->qtmap[level - 1];
				for (y = tly; y < bry; y++) {
					for (x = tlx; x < brx; x++) {
						qtmap[y][x] = 0;
					}
				}
				tly <<= 1;
				tlx <<= 1;
				bry <<= 1;
				brx <<= 1;
			}
#if AUTO_DEL_CL
			*blk = s_blk + 1;
			for (cl = 0; cl < enc->num_class; cl++) {
				enc->err_cost[s_blk][cl] = err_cost[cl];
			}
#endif
		} else {
			qtmap[y][x] = 1;
		}
	} else {
#if AUTO_DEL_CL
		*blk = s_blk + 1;
		for (cl = 0; cl < enc->num_class; cl++) {
			enc->err_cost[s_blk][cl] = err_cost[cl];
		}
#endif
	}
	free(err_cost);
	return (qtcost);
}

void count_cl(ENCODER *enc)
{
	int cl, y, x;

	for (cl = 0; cl < enc->num_class; cl++) enc->cl_hist[cl] = 0;
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			enc->cl_hist[(int)enc->class[y][x]]++;
		}
	}
}

cost_t optimize_class(ENCODER *enc)
{
	int y, x, i, blksize, level, blk = 0;
	int **prdbuf, **errbuf;
	enc->function_number = 3;
	if (enc->quadtree_depth >= 0 && enc->optimize_loop > 1) {
		level = enc->quadtree_depth;
		blksize = MAX_BSIZE;
	} else {
		level = 0;
		blksize = BASE_BSIZE;
	}
	for (i = 0; i < enc->num_class; i++) {
		enc->mtfbuf[i] = i;
	}
	prdbuf =(int **)alloc_2d_array(enc->num_class, blksize * blksize,
		sizeof(int));
	errbuf =(int **)alloc_2d_array(enc->num_class, blksize * blksize,
		sizeof(int));
	for (y = 0; y < enc->height; y += blksize) {
		for (x = 0; x < enc->width; x += blksize) {
			set_prdbuf(enc, prdbuf, errbuf, y, x, blksize);
			vbs_class(enc, prdbuf, errbuf, y, x,
				blksize, enc->width, level, &blk);
		}
	}
#if TEMPLATE_MATCHING_ON
	int cl;
	enc->temp_cl = -1;
	for(cl=0; cl<enc->num_class; cl++){
		if(enc->optimize_loop==1){
			if(enc->predictor[cl][0] == TEMPLATE_FLAG){
				enc->temp_cl = cl;
				break;
			} else {
				enc->temp_cl = -1;
			}
		} else if(enc->optimize_loop==2){
			if(enc->nzconv[cl][0] == -1){
				enc->temp_cl = cl;
				break;
			} else {
				enc->temp_cl = -1;
			}
		}
	}
	if(enc->optimize_loop == 2)	mask->temp_cl = enc->temp_cl;
#endif
#if CHECK_CLASS
	if(enc->optimize_loop == 1){
		for(y=0; y<enc->height; y += blksize){
			for(x=0; x<enc->width; x += blksize){
				printf("%2d ", enc->class[y][x]);
			}
			printf("\n");
		}
	}
#endif
	free(errbuf);
	free(prdbuf);
	count_cl(enc);
	return (calc_cost(enc, 0, 0, enc->height, enc->width));
}

void set_weight_flag(ENCODER *enc)
{
	int y, x, ty, tx, cl, i, sample;
	int count_cl[enc->num_class];

	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			sample = win_sample[(int)enc->mask[y][x]];	//ウィンドウ内の画素数(enc->maskの値に依存(1,9,25,49,81))
			for (cl = 0; cl < enc->num_class; cl++) {
				count_cl[cl]  = 0;
			}

			for (i = 0; i < sample; i++){
				ty = y + mask_y[i];
				tx = x + mask_x[i];	//マスクにかかる位置の画素を算出
				if(ty < 0) ty = 0;
				else if(ty >= enc->height) ty = enc->height - 1;
				if(tx < 0) tx = 0;
				else if(tx >= enc->width) tx = enc->width - 1;
				cl = enc->class[ty][tx];
				count_cl[cl]++;		//マスク内のクラスのヒストグラム
			}

			for(cl = 0; cl < enc->num_class; cl++){
				enc->weight[y][x][cl] = 0;
				if (count_cl[cl] != 0){
					enc->weight[y][x][cl] =( (count_cl[cl] << W_SHIFT) / sample);//確率モデルの重み
					//if(y==100&&x==100)printf("count_cl[cl] = %d sample = %d weight[y][x][%d] = %d\n",count_cl[cl],sample,cl,enc->weight[y][x][cl]);
				}
			}
		}
	}
}



#if AUTO_PRD_ORDER
/*****************************************************************************
enc->nzconv		: Non-Zero Coef to The Order (Index) table.
enc->num_nzcoef 	: Number of Non-Zero Coef
*****************************************************************************/
void set_prd_pels(ENCODER *enc)
{
	int cl, k, i;

	for (cl = 0; cl < enc->num_class; cl++) {
		k = 0;
		if(enc->nzconv[cl][0] == -1){
			/*for(i = 0; i < enc->max_prd_order; i++){
				enc->nzconv[cl][i] = TEMPLATE_FLAG;
			}*/
			k=-1;
		} else {
			for (i = 0; i < enc->max_prd_order; i++) {
				if (enc->predictor[cl][i] != 0) {
					enc->nzconv[cl][k] = i;
					k++;
				}
			}

			for (i = k; i < enc->max_prd_order; i++) {
				enc->nzconv[cl][i] = (int)1E5;			//for bugfix
			}
		}
		enc->num_nzcoef[cl] = k;
	}
}

/*****************************************************************************
Coef Modification
*****************************************************************************/
#if OPTIMIZE_MASK_LOOP
void optimize_coef(ENCODER *enc, int cl, int pos, int *num_eff)
{					//if AUTO_PRD_ORDER 1 / OPTIMIZE_MASK_LOOP 1
#define S_RANGE 2		//Search Range  ex. 2 -> +-1
#define SUBS_RANGE 2	//Search Range of One Coef Modification
	cost_t *cbuf, *cbuf_p, *cbuf_p2, c_base, *coef_cost_p1, *coef_cost_p2;
	double a;
	int i, j, k, l, x, y, df1, df2, df_f, *org_p, base, *roff_p;
	int coef_abs_pos, coef_abs, coef_pos;
	int num_swap_search, *swap_search;
	int prd, maxprd, *coef_p, prd_c;
	int u, peak, m_gr;
	char *class_p, onoff;
	img_t *bconv_p, *fconv_p;
	PMODEL *pm, *pm_p;
	int diffy, diffx, *diff, *diff_p, *diff_p2, sgn;

	num_swap_search = enc->max_prd_order - enc->num_nzcoef[cl];
	cbuf = (cost_t *)alloc_mem((SUBS_RANGE + (enc->max_prd_order * S_RANGE) + num_swap_search) * sizeof(cost_t));
	diff = (int *) alloc_mem ((SUBS_RANGE + (enc->max_prd_order * S_RANGE)+ num_swap_search) * sizeof(int));
	//table of exchange coef
	swap_search = (int *)alloc_mem(num_swap_search * sizeof(int));
	// nzc_p = enc->nzconv[cl];
	coef_p = enc->predictor[cl];
	coef_pos = coef_p[pos];
	coef_abs_pos = (coef_pos < 0) ? -coef_pos : coef_pos;
	cbuf_p = cbuf;
	cbuf_p2 = &cbuf[SUBS_RANGE + (enc->max_prd_order * S_RANGE)];
	diff_p = diff;
	diff_p2 = &diff[SUBS_RANGE + (enc->max_prd_order * S_RANGE)];
	i = enc->ord2mhd[pos];
	coef_cost_p1 = enc->coef_cost[enc->zero_m[i]][enc->coef_m[i]];
	for (i = 0; i < (SUBS_RANGE >> 1); i++) {
		df1 = coef_pos - (i + 1);
		for (j = 0; j < 2; j++) {
			y = df1;
			sgn = 0;
			if (y < 0) {
				y = -y;
				sgn = 1;
			}
			if (y > enc->max_coef) y = enc->max_coef;
			*cbuf_p++ = coef_cost_p1[y];
			*diff_p++ = (sgn) ? -y - coef_pos : y - coef_pos;
			df1 += (i + 1) << 1;
		}
	}
	k = 0;
	for (i = 0; i < enc->max_prd_order; i++) {
		onoff = (coef_p[i] == 0) ? 1 : 0;
		j = enc->ord2mhd[i];
		coef_cost_p2 = enc->coef_cost[enc->zero_m[j]][enc->coef_m[j]];
		df2 = coef_p[i];
		coef_abs = (df2 < 0) ? -df2 : df2;
		c_base = coef_cost_p2[coef_abs];
		for (l = 0; l < (S_RANGE >> 1); l++) {
			df1 = coef_pos - (l + 1);
			df2 = coef_p[i] + (l + 1);
			for (j = 0; j < 2; j++) {		//change "+ or -" modification
				y = df1;
				x = df2;
				sgn = 0;
				if (y < 0) {
					y = -y;
					sgn = 1;
				}
				diffy = y - enc->max_coef;
				if (diffy < 0) diffy = 0;
				if (y > enc->max_coef) y = enc->max_coef;
				if (x < 0) x = -x;
				diffx = x - enc->max_coef;
				if (diffx < 0) diffx = 0;
				if (x > enc->max_coef) x = enc->max_coef;
				if (diffy > 0 || diffx > 0) {
					if (diffy > diffx) {
						x += (j) ? diffy - diffx : diffx - diffy;
						if (x < 0) {
							x = -x;
						}
					} else {
						y += (j) ? diffy - diffx : diffx - diffy;
						if (y < 0) {
							y = -y;
							sgn = (sgn) ? 0 : 1;
						}
					}
				}
				*cbuf_p++ = coef_cost_p1[y] + coef_cost_p2[x] - c_base;
				*diff_p++ = (sgn) ? -y - coef_pos : y - coef_pos;
				df1 += (l + 1) << 1;
				df2 -= (l + 1) << 1;
			}
		}
		if (onoff == 1) {
			*cbuf_p2++ = coef_cost_p1[0] + coef_cost_p2[coef_abs_pos] - coef_cost_p2[0];
			*diff_p2++ = -coef_pos;
			swap_search[k++] = i;
		}
	}
	for (i = 0; i < S_RANGE; i++) {
		cbuf[SUBS_RANGE + pos * S_RANGE + i] = coef_cost_p1[coef_abs_pos];
		diff[SUBS_RANGE + pos * S_RANGE + i] = 0;
	}
	// before here -> the calculation of side cost
	bconv_p = enc->bconv;
	fconv_p = enc->fconv;
	maxprd = enc->maxprd;
	// shift = enc->coef_precision - 1;
  a = 1.0 / log(2.0);
	for (y = 0; y < enc->height; y++) {
		// class_p = enc->class[y];
		for (x = 0; x < enc->width; x++) {
			// if (cl != *class_p++) continue;
	  	if (enc->weight[y][x][cl] == 0) continue;
			u = enc->upara[y][x];
			peak = set_mask_parameter_optimize(enc, y, x, u, cl);
			roff_p = enc->roff[y][x];
			prd = enc->prd_class[y][x][cl];
			org_p = &enc->org[y][x];
			pm_p = enc->pmlist[(int)enc->group[y][x]];
			df1 = org_p[roff_p[pos]];
			cbuf_p = cbuf;
			cbuf_p2 = &cbuf[SUBS_RANGE + (enc->max_prd_order * S_RANGE)];
			diff_p = diff;
			diff_p2 = &diff[SUBS_RANGE + (enc->max_prd_order * S_RANGE)];
			df_f = df1;
			//Only One Coef Modification
			for (i = 0; i < (SUBS_RANGE >> 1); i++) {
				for (j = 0; j < 2; j++) {	//change "+ or -" modification
					prd_c = prd + df_f * (*diff_p++);
					prd_c = CLIP(0, maxprd, prd_c);
					if (mask->num_peak == 1){
						base = bconv_p[prd_c];
						pm = pm_p + fconv_p[prd_c];
						*cbuf_p++ += pm->cost[*org_p + base]
							+ pm->subcost[base];
					}else{
						mask->base[peak] = bconv_p[prd_c];
						m_gr = enc->uquant[cl][u];
						mask->pm[peak] = enc->pmlist[m_gr] + fconv_p[prd_c];
						set_pmodel_mult_cost(mask,enc->maxval+1,*org_p);
						*cbuf_p++ += a * (log(mask->cumfreq)-log(mask->freq));
					}
				}
			}
			//Modification and Search The Other Coef
			for (i = 0; i < enc->max_prd_order; i++) {
				df2 = org_p[roff_p[i]];
				df_f = df1 - df2;
				for (l = 0; l < (S_RANGE >> 1); l++) {
					for (j = 0; j < 2; j++) {	//change "+ or -"
						prd_c = prd + df_f * (*diff_p++);
						prd_c = CLIP(0, maxprd, prd_c);
						if (mask->num_peak == 1){
							base = bconv_p[prd_c];
							pm = pm_p + fconv_p[prd_c];
							*cbuf_p++ += pm->cost[*org_p + base]
								+ pm->subcost[base];
						}else{
							mask->base[peak] = bconv_p[prd_c];
							m_gr = enc->uquant[cl][u];
							mask->pm[peak] = enc->pmlist[m_gr] + fconv_p[prd_c];
							set_pmodel_mult_cost(mask,enc->maxval+1,*org_p);
							*cbuf_p++ += a * (log(mask->cumfreq)-log(mask->freq));
						}
					}
				}
			}
			//exchange two coefs
			for (j = 0; j < num_swap_search; j++) {
				k = swap_search[j];
				prd_c = prd + (df1 - org_p[roff_p[k]]) * (*diff_p2++);
				prd_c = CLIP(0, maxprd, prd_c);
				if (mask->num_peak == 1){
					base = bconv_p[prd_c];
					pm = pm_p + fconv_p[prd_c];
					*cbuf_p2++ += pm->cost[*org_p + base]
						+ pm->subcost[base];
				}else{
					mask->base[peak] = bconv_p[prd_c];
					m_gr = enc->uquant[cl][u];
					mask->pm[peak] = enc->pmlist[m_gr] + fconv_p[prd_c];
					set_pmodel_mult_cost(mask,enc->maxval+1,*org_p);
					*cbuf_p2++ += a * (log(mask->cumfreq)-log(mask->freq));
				}
			}
		}
	}

	j = SUBS_RANGE + pos * S_RANGE;
	for (i = 0; i < SUBS_RANGE + (enc->max_prd_order * S_RANGE) + num_swap_search; i++) {
		if (cbuf[i] < cbuf[j]) {
			j = i;
		}
	}
	free(cbuf);
	if (j == SUBS_RANGE + pos * S_RANGE) {
		free(swap_search);
		free(diff);
		return;
	}
	if (j < SUBS_RANGE) {
		for (y = 0; y < enc->height; y++) {
			class_p = enc->class[y];
			for (x = 0; x < enc->width; x++) {
				if (cl == enc->temp_cl ){
					enc->prd_class[y][x][cl] = exam_array[y][x][0];
				} else if (cl == *class_p++) {
					org_p = &enc->org[y][x];
					roff_p = enc->roff[y][x];
					enc->prd[y][x] += org_p[roff_p[pos]] * diff[j];
					enc->prd_class[y][x][cl] += org_p[roff_p[pos]] * diff[j];
				}
			}
		}
		coef_p[pos] += diff[j];
	} else {
		if (j < SUBS_RANGE + (enc->max_prd_order * S_RANGE)) {
			i = j - SUBS_RANGE;
			i /= S_RANGE;
		} else {
			i = j - SUBS_RANGE - (enc->max_prd_order * S_RANGE);
			i = swap_search[i];
		}
		for (y = 0; y < enc->height; y++) {
			class_p = enc->class[y];
			for (x = 0; x < enc->width; x++) {
			#if TEMPLATE_MATCHING_ON
				if (cl == enc->temp_cl){
					enc->prd_class[y][x][cl] = exam_array[y][x][0];
				} else
			#endif
				if (cl == *class_p++) {
					org_p = &enc->org[y][x];
					roff_p = enc->roff[y][x];
					enc->prd[y][x] += (org_p[roff_p[pos]] - org_p[roff_p[i]]) * diff[j];
					enc->prd_class[y][x][cl] += (org_p[roff_p[pos]] - org_p[roff_p[i]]) * diff[j];
				}
			}
		}
		coef_p[pos] += diff[j];
		coef_p[i] -= diff[j];
	}
	free(swap_search);
	free(diff);
	if (diff[j] != 0) (*num_eff)++;
}



#else
//if AUTO_PRD_ORDER	1 / OPTIMIZE_MASK_LOOP	0
void optimize_coef(ENCODER *enc, int cl, int pos, int *num_eff)
{
#define S_RANGE 2		//Search Range  ex. 2 -> +-1
#define SUBS_RANGE 2	//Search Range of One Coef Modification
	cost_t *cbuf, *cbuf_p, *cbuf_p2, c_base, *coef_cost_p1, *coef_cost_p2;
	int i, j, k, l, x, y, df1, df2, df_f, *org_p, base, *roff_p;
	int coef_abs_pos, coef_abs, coef_pos;
	int num_swap_search, *swap_search;
	int prd, shift, maxprd, *coef_p, *nzc_p, prd_c;
	char *class_p, onoff;
	img_t *bconv_p, *fconv_p;
	PMODEL *pm, *pm_p;
	int diffy, diffx, *diff, *diff_p, *diff_p2, sgn;

	num_swap_search = enc->max_prd_order - enc->num_nzcoef[cl];
	cbuf = (cost_t *)alloc_mem((SUBS_RANGE + (enc->max_prd_order * S_RANGE) + num_swap_search) * sizeof(cost_t));
	diff = (int *) alloc_mem ((SUBS_RANGE + (enc->max_prd_order * S_RANGE)+ num_swap_search) * sizeof(int));
	//table of exchange coef
	swap_search = (int *)alloc_mem(num_swap_search * sizeof(int));
	nzc_p = enc->nzconv[cl];
	coef_p = enc->predictor[cl];
	coef_pos = coef_p[pos];
	coef_abs_pos = (coef_pos < 0) ? -coef_pos : coef_pos;
	cbuf_p = cbuf;
	cbuf_p2 = &cbuf[SUBS_RANGE + (enc->max_prd_order * S_RANGE)];
	diff_p = diff;
	diff_p2 = &diff[SUBS_RANGE + (enc->max_prd_order * S_RANGE)];
	i = enc->ord2mhd[pos];
	coef_cost_p1 = enc->coef_cost[enc->zero_m[i]][enc->coef_m[i]];
	for (i = 0; i < (SUBS_RANGE >> 1); i++) {
		df1 = coef_pos - (i + 1);
		for (j = 0; j < 2; j++) {
			y = df1;
			sgn = 0;
			if (y < 0) {
				y = -y;
				sgn = 1;
			}
			if (y > enc->max_coef) y = enc->max_coef;
			*cbuf_p++ = coef_cost_p1[y];
			*diff_p++ = (sgn) ? -y - coef_pos : y - coef_pos;
			df1 += (i + 1) << 1;
		}
	}
	k = 0;
	for (i = 0; i < enc->max_prd_order; i++) {
		onoff = (coef_p[i] == 0) ? 1 : 0;
		j = enc->ord2mhd[i];
		coef_cost_p2 = enc->coef_cost[enc->zero_m[j]][enc->coef_m[j]];
		df2 = coef_p[i];
		coef_abs = (df2 < 0) ? -df2 : df2;
		c_base = coef_cost_p2[coef_abs];
		for (l = 0; l < (S_RANGE >> 1); l++) {
			df1 = coef_pos - (l + 1);
			df2 = coef_p[i] + (l + 1);
			for (j = 0; j < 2; j++) {		//change "+ or -" modification
				y = df1;
				x = df2;
				sgn = 0;
				if (y < 0) {
					y = -y;
					sgn = 1;
				}
				diffy = y - enc->max_coef;
				if (diffy < 0) diffy = 0;
				if (y > enc->max_coef) y = enc->max_coef;
				if (x < 0) x = -x;
				diffx = x - enc->max_coef;
				if (diffx < 0) diffx = 0;
				if (x > enc->max_coef) x = enc->max_coef;
				if (diffy > 0 || diffx > 0) {
					if (diffy > diffx) {
						x += (j) ? diffy - diffx : diffx - diffy;
						if (x < 0) {
							x = -x;
						}
					} else {
						y += (j) ? diffy - diffx : diffx - diffy;
						if (y < 0) {
							y = -y;
							sgn = (sgn) ? 0 : 1;
						}
					}
				}
				*cbuf_p++ = coef_cost_p1[y] + coef_cost_p2[x] - c_base;
				*diff_p++ = (sgn) ? -y - coef_pos : y - coef_pos;
				df1 += (l + 1) << 1;
				df2 -= (l + 1) << 1;
			}
		}
		if (onoff == 1) {
			*cbuf_p2++ = coef_cost_p1[0] + coef_cost_p2[coef_abs_pos] - coef_cost_p2[0];
			*diff_p2++ = -coef_pos;
			swap_search[k++] = i;
		}
	}
	for (i = 0; i < S_RANGE; i++) {
		cbuf[SUBS_RANGE + pos * S_RANGE + i] = coef_cost_p1[coef_abs_pos];
		diff[SUBS_RANGE + pos * S_RANGE + i] = 0;
	}
	// before here -> the calculation of side cost
	bconv_p = enc->bconv;
	fconv_p = enc->fconv;
	maxprd = enc->maxprd;
	shift = enc->coef_precision - 1;
	for (y = 0; y < enc->height; y++) {
		class_p = enc->class[y];
		for (x = 0; x < enc->width; x++) {
			if (cl != *class_p++) continue;
			roff_p = enc->roff[y][x];
			prd = enc->prd[y][x];
			org_p = &enc->org[y][x];
			pm_p = enc->pmlist[(int)enc->group[y][x]];
			df1 = org_p[roff_p[pos]];
			cbuf_p = cbuf;
			cbuf_p2 = &cbuf[SUBS_RANGE + (enc->max_prd_order * S_RANGE)];
			diff_p = diff;
			diff_p2 = &diff[SUBS_RANGE + (enc->max_prd_order * S_RANGE)];
			df_f = df1;
			//Only One Coef Modification
			for (i = 0; i < (SUBS_RANGE >> 1); i++) {
				for (j = 0; j < 2; j++) {	//change "+ or -" modification
					prd_c = prd + df_f * (*diff_p++);
					prd_c = CLIP(0, maxprd, prd_c);
					base = bconv_p[prd_c];
					pm = pm_p + fconv_p[prd_c];
					*cbuf_p++ += pm->cost[*org_p + base]
					+ pm->subcost[base];
				}
			}
			//Modification and Search The Other Coef
			for (i = 0; i < enc->max_prd_order; i++) {
				df2 = org_p[roff_p[i]];
				df_f = df1 - df2;
				for (l = 0; l < (S_RANGE >> 1); l++) {
					for (j = 0; j < 2; j++) {	//change "+ or -"
						prd_c = prd + df_f * (*diff_p++);
						prd_c = CLIP(0, maxprd, prd_c);
						base = bconv_p[prd_c];
						pm = pm_p + fconv_p[prd_c];
						*cbuf_p++ += pm->cost[*org_p + base]
						+ pm->subcost[base];
					}
				}
			}
			//exchange two coefs
			for (j = 0; j < num_swap_search; j++) {
				k = swap_search[j];
				prd_c = prd + (df1 - org_p[roff_p[k]]) * (*diff_p2++);
				prd_c = CLIP(0, maxprd, prd_c);
				base = bconv_p[prd_c];
				pm = pm_p + fconv_p[prd_c];
				*cbuf_p2++ += pm->cost[*org_p + base]
				+ pm->subcost[base];
			}
		}
	}
	j = SUBS_RANGE + pos * S_RANGE;
	for (i = 0; i < SUBS_RANGE + (enc->max_prd_order * S_RANGE) + num_swap_search; i++) {
		if (cbuf[i] < cbuf[j]) {
			j = i;
		}
	}
	free(cbuf);
	if (j == SUBS_RANGE + pos * S_RANGE) {
		free(swap_search);
		free(diff);
		return;
	}
	if (j < SUBS_RANGE) {
		for (y = 0; y < enc->height; y++) {
			class_p = enc->class[y];
			for (x = 0; x < enc->width; x++) {
				if (cl == *class_p++) {
					org_p = &enc->org[y][x];
					roff_p = enc->roff[y][x];
					enc->prd[y][x] += org_p[roff_p[pos]] * diff[j];
				}
			}
		}
		coef_p[pos] += diff[j];
	} else {
		if (j < SUBS_RANGE + (enc->max_prd_order * S_RANGE)) {
			i = j - SUBS_RANGE;
			i /= S_RANGE;
		} else {
			i = j - SUBS_RANGE - (enc->max_prd_order * S_RANGE);
			i = swap_search[i];
		}
		for (y = 0; y < enc->height; y++) {
			class_p = enc->class[y];
			for (x = 0; x < enc->width; x++) {
				if (cl == *class_p++) {
					org_p = &enc->org[y][x];
					roff_p = enc->roff[y][x];
					enc->prd[y][x] += (org_p[roff_p[pos]] - org_p[roff_p[i]]) * diff[j];
				}
			}
		}
		coef_p[pos] += diff[j];
		coef_p[i] -= diff[j];
	}
	free(swap_search);
	free(diff);
	if (diff[j] != 0) (*num_eff)++;
}

#endif

cost_t optimize_predictor(ENCODER *enc)	//when AUTO_PRD_ORDER 1
{
	int cl, pos, k, num_eff, y, x;
#ifndef RAND_MAX
#  define RAND_MAX 32767
#endif
	enc->function_number = 4;
	for (cl = 0; cl < enc->num_class; cl++) {
		num_eff = 0;
		if (enc->cl_hist[cl] == 0) continue;

		if(enc->num_nzcoef[cl] != -1){
			for (k = 0; k < enc->num_search[cl]; k++) {
				if (enc->num_nzcoef[cl] == 0) continue;
				pos = (int)(((double)rand() * enc->num_nzcoef[cl]) / (RAND_MAX+1.0));	//	いじる係数をrand関数で選択
				pos = enc->nzconv[cl][pos];
				optimize_coef(enc, cl, pos, &num_eff);	//予測係数の最適化
				set_prd_pels(enc);
			}
		} else {
			for(y=0; y<enc->height; y++){
				for(x=0; x<enc->width; x++){
					if(enc->class[y][x] == cl){
						enc->prd_class[y][x][cl] = exam_array[y][x][0];
					}
				}
			}
			set_prd_pels(enc);
		}
		enc->num_search[cl] = num_eff + 3;

#if CHECK_PREDICTOR
		printf("[%2d] ", cl);
		if(enc->num_nzcoef[cl] == -1){
			for(k=0; k<5; k++){
				printf("%2d," ,enc->predictor[cl][enc->nzconv[cl][k]]);
			}
		} else {
			for(k=0; k<enc->num_nzcoef[cl]; k++){
				printf("%2d," ,enc->predictor[cl][enc->nzconv[cl][k]]);
			}
		}
		printf("\n");
#endif

	}
	save_prediction_value(enc);
	predict_region(enc, 0, 0, enc->height, enc->width);
	return (calc_cost(enc, 0, 0, enc->height, enc->width));
}

#else 	//AUTO_PRD_ORDER	0

#if OPTIMIZE_MASK_LOOP
void optimize_coef(ENCODER *enc, int cl, int pos1, int pos2)
{		//if AUTO_PRD_ORDER	0 / OPTIMIZE_MASK_LOOP 1
#define SEARCH_RANGE 11
#define SUBSEARCH_RANGE 3
	cost_t cbuf[SEARCH_RANGE * SUBSEARCH_RANGE], *cbuf_p;
	double a;
	int i, j, k, x, y, df1, df2, base;
	int prd, prd_f, shift, maxprd, *coef_p, *roff_p, *org_p;
	int u, peak, m_gr;
	char *class_p;
	img_t *bconv_p, *fconv_p;
	PMODEL *pm, *pm_p;

	cbuf_p = cbuf;
	coef_p = enc->predictor[cl];
	k = 0;
	for (i = 0; i < SEARCH_RANGE; i++) {
		y = coef_p[pos1] + i - (SEARCH_RANGE >> 1);
		if (y < 0) y = -y;
		if (y > enc->max_coef) y = enc->max_coef;
		for (j = 0; j < SUBSEARCH_RANGE; j++) {
			x = coef_p[pos2] - (i - (SEARCH_RANGE >> 1))
				- (j - (SUBSEARCH_RANGE >> 1));
			if (x < 0) x = -x;
			if (x > enc->max_coef) x = enc->max_coef;

			cbuf_p[k++] = enc->coef_cost[enc->coef_m[pos1]][y]
				+ enc->coef_cost[enc->coef_m[pos2]][x];

		}
	}
	bconv_p = enc->bconv;
	fconv_p = enc->fconv;
	maxprd = enc->maxprd;
	shift = enc->coef_precision - 1;
  a = 1.0 / log(2.0);
	for (y = 0; y < enc->height; y++) {
//		class_p = enc->class[y];
		for (x = 0; x < enc->width; x++) {
			//		if (cl != *class_p++) continue;
	  	if (enc->weight[y][x][cl] == 0) continue;
			u = enc->upara[y][x];
			peak = set_mask_parameter_optimize(enc, y, x, u, cl);
			roff_p = enc->roff[y][x];
			prd = enc->prd_class[y][x][cl];
			org_p = &enc->org[y][x];
			pm_p = enc->pmlist[(int)enc->group[y][x]];
			df1 = org_p[roff_p[pos1]];
			df2 = org_p[roff_p[pos2]];
			prd_f = prd - (df1 - df2) * (SEARCH_RANGE >> 1)
				+ df2 * (SUBSEARCH_RANGE >> 1);
			cbuf_p = cbuf;
			for (i = 0; i < SEARCH_RANGE; i++) {
				for (j = 0; j < SUBSEARCH_RANGE; j++) {
					prd = prd_f;
					prd = CLIP(0, maxprd, prd);
					if (mask->num_peak == 1){
						base = bconv_p[prd];
						pm = pm_p + fconv_p[prd];
						(*cbuf_p++) += pm->cost[*org_p + base]
							+ pm->subcost[base];
					}else{
						mask->base[peak] = bconv_p[prd];
						m_gr = enc->uquant[cl][u];
						mask->pm[peak] = enc->pmlist[m_gr] + fconv_p[prd];
						set_pmodel_mult_cost(mask,enc->maxval+1,*org_p);
						(*cbuf_p++) += a * (log(mask->cumfreq)-log(mask->freq));
					}
					prd_f -= df2;
				}
				prd_f += df1 + df2 * (SUBSEARCH_RANGE - 1);
			}
		}
	}
	//						printf("b");

	cbuf_p = cbuf;
	j = (SEARCH_RANGE * SUBSEARCH_RANGE) >> 1;
	for (i = 0; i < SEARCH_RANGE * SUBSEARCH_RANGE; i++) {
		if (cbuf_p[i] < cbuf_p[j]) {
			j = i;
		}
	}
	i = (j / SUBSEARCH_RANGE) - (SEARCH_RANGE >> 1);
	j = (j % SUBSEARCH_RANGE) - (SUBSEARCH_RANGE >> 1);
	y = coef_p[pos1] + i;
	x = coef_p[pos2] - i - j;
	if (y < -enc->max_coef) y = -enc->max_coef;
	else if (y > enc->max_coef) y = enc->max_coef;
	if (x < -enc->max_coef) x = -enc->max_coef;
	else if (x > enc->max_coef) x = enc->max_coef;
	i = y - coef_p[pos1];
	j = x - coef_p[pos2];
	if (i != 0 || j != 0) {
		for (y = 0; y < enc->height; y++) {
			class_p = enc->class[y];
			for (x = 0; x < enc->width; x++) {
				if (cl == enc->temp_cl){
					enc->prd_class[y][x][cl] = exam_array[y][x][0];
				} else if (cl == *class_p++) {
					roff_p = enc->roff[y][x];
					org_p = &enc->org[y][x];
//					enc->prd[y][x] += org_p[roff_p[pos1]] * i
//						+ org_p[roff_p[pos2]] * j;
					enc->prd_class[y][x][cl] += org_p[roff_p[pos1]] * i
						+ org_p[roff_p[pos2]] * j;
				}
			}
		}
		coef_p[pos1] += i;
		coef_p[pos2] += j;
	}
}

#else

void optimize_coef(ENCODER *enc, int cl, int pos1, int pos2)
{				//if AUTO_PRD_ORDER	0 / OPTIMIZE_MASK_LOOP 0
#define SEARCH_RANGE 11
#define SUBSEARCH_RANGE 3
	cost_t cbuf[SEARCH_RANGE * SUBSEARCH_RANGE], *cbuf_p;
	int i, j, k, x, y, df1, df2, base;
	int prd, prd_f, shift, maxprd, *coef_p, *roff_p, *org_p;
	char *class_p;
	img_t *bconv_p, *fconv_p;
	PMODEL *pm, *pm_p;

	cbuf_p = cbuf;
	coef_p = enc->predictor[cl];
	k = 0;
	for (i = 0; i < SEARCH_RANGE; i++) {
		y = coef_p[pos1] + i - (SEARCH_RANGE >> 1);
		if (y < 0) y = -y;
		if (y > enc->max_coef) y = enc->max_coef;
		for (j = 0; j < SUBSEARCH_RANGE; j++) {
			x = coef_p[pos2] - (i - (SEARCH_RANGE >> 1))
				- (j - (SUBSEARCH_RANGE >> 1));
			if (x < 0) x = -x;
			if (x > enc->max_coef) x = enc->max_coef;

			cbuf_p[k++] = enc->coef_cost[enc->coef_m[pos1]][y]
			+ enc->coef_cost[enc->coef_m[pos2]][x];

		}
	}
	bconv_p = enc->bconv;
	fconv_p = enc->fconv;
	maxprd = enc->maxprd;
	shift = enc->coef_precision - 1;
	for (y = 0; y < enc->height; y++) {
		class_p = enc->class[y];
		for (x = 0; x < enc->width; x++) {
			if (cl != *class_p++) continue;
			roff_p = enc->roff[y][x];
			prd = enc->prd[y][x];
			org_p = &enc->org[y][x];
			pm_p = enc->pmlist[(int)enc->group[y][x]];
			df1 = org_p[roff_p[pos1]];
			df2 = org_p[roff_p[pos2]];
			prd_f = prd - (df1 - df2) * (SEARCH_RANGE >> 1)
				+ df2 * (SUBSEARCH_RANGE >> 1);
			cbuf_p = cbuf;
			for (i = 0; i < SEARCH_RANGE; i++) {
				for (j = 0; j < SUBSEARCH_RANGE; j++) {
					prd = prd_f;
					prd = CLIP(0, maxprd, prd);
					base = bconv_p[prd];
					pm = pm_p + fconv_p[prd];
					(*cbuf_p++) += pm->cost[*org_p + base]
					+ pm->subcost[base];
					prd_f -= df2;
				}
				prd_f += df1 + df2 * (SUBSEARCH_RANGE - 1);
			}
		}
	}
	cbuf_p = cbuf;
	j = (SEARCH_RANGE * SUBSEARCH_RANGE) >> 1;
	for (i = 0; i < SEARCH_RANGE * SUBSEARCH_RANGE; i++) {
		if (cbuf_p[i] < cbuf_p[j]) {
			j = i;
		}
	}
	i = (j / SUBSEARCH_RANGE) - (SEARCH_RANGE >> 1);
	j = (j % SUBSEARCH_RANGE) - (SUBSEARCH_RANGE >> 1);
	y = coef_p[pos1] + i;
	x = coef_p[pos2] - i - j;
	if (y < -enc->max_coef) y = -enc->max_coef;
	else if (y > enc->max_coef) y = enc->max_coef;
	if (x < -enc->max_coef) x = -enc->max_coef;
	else if (x > enc->max_coef) x = enc->max_coef;
	i = y - coef_p[pos1];
	j = x - coef_p[pos2];
	if (i != 0 || j != 0) {
		for (y = 0; y < enc->height; y++) {
			class_p = enc->class[y];
			for (x = 0; x < enc->width; x++) {
				if (cl == *class_p++) {
					roff_p = enc->roff[y][x];
					org_p = &enc->org[y][x];
					enc->prd[y][x] += org_p[roff_p[pos1]] * i
						+ org_p[roff_p[pos2]] * j;
				}
			}
		}
		coef_p[pos1] += i;
		coef_p[pos2] += j;
	}
}
#endif

cost_t optimize_predictor(ENCODER *enc)	//when AUTO_PRD_ORDER 0
{
	int cl, k, pos1, pos2;
#ifndef RAND_MAX
#  define RAND_MAX 32767
#endif

	for (cl = 0; cl < enc->num_class; cl++) {
		//if (enc->cl_hist[cl] == 0) continue;
		for (k = 0; k < enc->max_prd_order; k++) {
retry:
			pos1 = (int)(((double)rand() * enc->max_prd_order) / (RAND_MAX+1.0));
			pos2 = (int)(((double)rand() * enc->max_prd_order) / (RAND_MAX+1.0));
			if (pos1 == pos2) goto retry;
			optimize_coef(enc, cl, pos1, pos2);
		}
	}
	save_prediction_value(enc);
	predict_region(enc, 0, 0, enc->height, enc->width);
	return (calc_cost(enc, 0, 0, enc->height, enc->width));
}
#endif

#if MULT_PEAK_MODE
cost_t optimize_mask(ENCODER *enc)
{
	cost_t cost, min_cost;
	int y, x, ty, tx, min_m, m, i, j;
	enc->function_number = 7;
//	min_cost = 1E8;

	for (y = 0; y < enc->height; y += WIN_BSIZE) {
		for (x = 0; x < enc->width; x += WIN_BSIZE) {
			if ((y + WIN_BSIZE) > enc->height) {
				ty = enc->height % WIN_BSIZE;
			} else {
				ty = WIN_BSIZE;
			}
			if ((x + WIN_BSIZE) > enc->width) {
				tx = enc->width % WIN_BSIZE;
			} else {
				tx = WIN_BSIZE;
			}
			min_m = 0;
			min_cost = 1E8;
			for (m = 0; m < NUM_MASK; m++) {
				for (i = 0; i < ty; i++) {
					for (j = 0; j < tx; j++) {
						enc->mask[y + i][x + j] = m;
					}
				}
				cost = calc_cost2(enc, y, x, y + ty, x + tx);
				if (cost < min_cost) {
					min_m = m;
					min_cost = cost;
				}
			}

			for(i = 0; i < ty; i++){
				for(j = 0; j < tx; j++){
					enc->mask[y + i][x + j] = min_m;
				}
			}
//			printf("min_m = %d\n",min_m);
		}
	}
	cost = calc_cost2(enc, 0, 0, enc->height, enc->width);
	return(cost);
}

void set_mask_parameter_optimize2(ENCODER *enc,int y, int x, int u)
{
	int cl, peak, gr;

	for(cl = peak = 0; cl < enc->num_class; cl++){
		if (enc->weight[y][x][cl] != 0){
		#if TEMPLATE_MATCHING_ON
			if(cl == enc->temp_cl){
				for(gr=0; gr<enc->num_group; gr++){
					if(u < enc->w_gr[gr])	break;
				}
				peak = temp_mask_parameter(enc, y, x, u, peak, cl, enc->weight[y][x][cl], gr);
			} else {
		#endif
				mask->class[peak] = cl;
				mask->weight[peak] = enc->weight[y][x][cl];
				peak++;
		#if TEMPLATE_MATCHING_ON
			}
		#endif
		}
	}

	mask->num_peak = peak;	//ピークの数
}

#if TEMPLATE_MATCHING_ON
int set_mask_parameter_optimize_temp(ENCODER *enc,int y, int x, int u, int r_cl, int w_gr)
{
	int cl, peak;
	int m_gr,m_prd,m_frac,r_peak;
	r_peak = -1;

	for(cl = peak = 0; cl < enc->num_class; cl++){
		if (enc->weight[y][x][cl] != 0){
			if(cl == enc->temp_cl){
				peak = temp_mask_parameter(enc, y, x, u, peak, cl, enc->weight[y][x][cl], w_gr);
				if(cl == r_cl)	r_peak = peak - (enc->temp_peak_num - 1);	//マッチングコストが一番小さい事例のピーク番号
			} else {
				mask->class[peak] = cl;
				mask->weight[peak] = enc->weight[y][x][cl];
				m_prd = enc->prd_class[y][x][cl];
				m_prd = CLIP(0, enc->maxprd, m_prd);
				mask->base[peak] = enc->bconv[m_prd];
				m_frac = enc->fconv[m_prd];
				m_gr = enc->uquant[cl][u];
				#if CHECK_DEBUG
					if( y == check_y && x == check_x && enc->function_number == F_NUM)	printf("[set_mask_parameter_opt] m_prd[%d]: %d[%2d] | weight: %d | gr: %d\n", 	peak, m_prd, cl, mask->weight[peak], m_gr);
				#endif
				if (cl == r_cl)	 r_peak = peak;	//当該ブロックのピーク番号
				mask->pm[peak] = enc->pmlist[m_gr] + m_frac;
				peak++;
			}
		}
	}

	mask->num_peak = peak;	//ピークの数
	return(r_peak);
}

void make_th(ENCODER *enc){
	int gr=0, y, x, u, peak, e, prd, th1, th0, frac, **trellis, cl=0, base, k;
	double a;
	cost_t cost, **cbuf, *cbuf_p, *dpcost, *thc_p, min_cost;
	PMODEL *pm;

	trellis = (int **)alloc_2d_array(enc->num_group, MAX_UPARA + 2,
		sizeof(int));
	dpcost = (cost_t *)alloc_mem((MAX_UPARA + 2) * sizeof(cost_t));
	cbuf = (cost_t **)alloc_2d_array(enc->num_group, MAX_UPARA + 2,
		sizeof(cost_t));
	thc_p = enc->th_cost;	//閾値の符号化に必要な符号量
	a = 1.0 / log(2.0);

	for (gr = 0; gr < enc->num_group; gr++) {
		cbuf_p = cbuf[gr];
		for (u = 0; u < MAX_UPARA + 2; u++) {
			cbuf_p[u] = 0;
		}
	}
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			// if (enc->weight[y][x][cl] == 0) continue;
			if (enc->class[y][x] != enc->temp_cl) continue;
			cl = enc->class[y][x];
			u = enc->upara[y][x] + 1;
			// peak = set_mask_parameter_optimize(enc, y, x, u-1,cl);
			e = enc->encval[y][x];
			prd = enc->prd_class[y][x][cl];
			prd = CLIP(0, enc->maxprd, prd);
			frac = enc->fconv[prd];
			for (gr = 0; gr < enc->num_group; gr++) {
				peak = set_mask_parameter_optimize_temp(enc, y, x, u-1,cl, gr);
				// mask->pm[peak] = enc->pmlist[gr] + frac;
				if (mask->num_peak == 1){
					base = mask->base[0];
					pm = mask->pm[0];
					cbuf[gr][u] += pm->cost[base + e] + pm->subcost[base];
				}else{
					set_pmodel_mult_cost(mask,enc->maxval+1,e);
					cbuf[gr][u] += a * (log(mask->cumfreq)-log(mask->freq));
				}
			}
		}
	}

	for (gr = 0; gr < enc->num_group; gr++) {
		cbuf_p = cbuf[gr];
		for (u = 1; u < MAX_UPARA + 2; u++) {
			cbuf_p[u] += cbuf_p[u - 1];	//cbuf_p:ある特徴量までの特徴量の総和が格納
		}
	}
	cbuf_p = cbuf[0];
	for (u = 0; u < MAX_UPARA + 2; u++) {
		dpcost[u] = cbuf_p[u] + thc_p[u];
	}
	for (gr = 1; gr < enc->num_group - 1; gr++) {
		cbuf_p = cbuf[gr];
		/* minimize (cbuf_p[th1] - cbuf_p[th0] + dpcost[th0]) */
		for (th1 = MAX_UPARA + 1; th1 >= 0; th1--) {
			th0 = th1;
			min_cost = dpcost[th1] - cbuf_p[th1] + thc_p[0];
			for (k = 0; k < th1; k++) {
				cost = dpcost[k] - cbuf_p[k] + thc_p[th1 - k];
				if (cost < min_cost) {
					min_cost = cost;
					th0 = k;
				}
			}
			dpcost[th1] = min_cost + cbuf_p[th1];
			trellis[gr][th1] = th0;
		}
	}
	cbuf_p = cbuf[gr];
	/* minimize (cbuf_p[th1] - cbuf_p[th0] + dpcost[th0]) for last group */
	th1 = MAX_UPARA + 1;
	th0 = th1;
	min_cost = dpcost[th1] - cbuf_p[th1];
	for (k = 0; k < th1; k++) {
		cost = dpcost[k] - cbuf_p[k];
		if (cost < min_cost) {
			min_cost = cost;
			th0 = k;
		}
	}
	trellis[gr][th1] = th0;
	for (gr = enc->num_group - 1; gr > 0; gr--) {
		th1 = trellis[gr][th1];
		// enc->th[cl][gr - 1] = th1;
		enc->w_gr[gr - 1] = th1;
	}
	// enc->w_gr[0] = 0;
	// enc->w_gr[enc->num_group-1] = MAX_UPARA+1;
	return;
}

cost_t optimize_template(ENCODER *enc){
	char temp_num;
	int min_temp_num=0, i;
	cost_t cost, min_cost = INT_MAX;
	enc->function_number = 6;
	if(enc->temp_cl == -1)	return(calc_cost2(enc, 0, 0, enc->height, enc->width));

// コンテクスト毎の重み
	make_th(enc);
	/* for(i=0; i<enc->num_group; i++){
		printf("%d, %d\n", i, enc->w_gr[i]);
	} */

// ピークの数の最適化
	if(enc->cl_hist[enc->temp_cl] == 0)	return(calc_cost2(enc, 0, 0, enc->height, enc->width));

	for(temp_num = 1; temp_num<=TEMPLATE_CLASS_NUM; temp_num++){
		enc->temp_peak_num = temp_num;
		// for(gr=0; gr<enc->num_group; gr++){
			// enc->w_gr = gr;
			cost = calc_cost2(enc, 0, 0, enc->height, enc->width);
			// printf("%d,%f\n", temp_num, cost);
			if(cost < min_cost){
				min_cost = cost;
				min_temp_num = temp_num;
				// min_gr = gr;
			}
		// }
	}
	if(min_temp_num <= 0 || min_temp_num > TEMPLATE_CLASS_NUM)	min_temp_num = TEMPLATE_CLASS_NUM;
	enc->temp_peak_num = min_temp_num;

	return(min_cost);
}
#endif

cost_t optimize_group_mult(ENCODER *enc)
{
	cost_t cost, min_cost, **cbuf, *dpcost, *cbuf_p, *thc_p;
	double a;
	int x, y, th1, th0, k, u, cl, gr, prd, e, base, frac, peak,m_gr,count;
	int **trellis;
	enc->function_number = 5;
	PMODEL *pm, **pm_p;

	trellis = (int **)alloc_2d_array(enc->num_group, MAX_UPARA + 2,
		sizeof(int));
	dpcost = (cost_t *)alloc_mem((MAX_UPARA + 2) * sizeof(cost_t));
	cbuf = (cost_t **)alloc_2d_array(enc->num_group, MAX_UPARA + 2,
		sizeof(cost_t));
	thc_p = enc->th_cost;
	a = 1.0 / log(2.0);

	for (k = 0; k < MAX_UPARA + 2; k++) trellis[0][k] = 0;
	/* Dynamic programming */
	for (cl = 0; cl < enc->num_class; cl++) {
		//if (enc->cl_hist[cl] == 0) continue;
		for (gr = 0; gr < enc->num_group; gr++) {
			cbuf_p = cbuf[gr];
			for (u = 0; u < MAX_UPARA + 2; u++) {
				cbuf_p[u] = 0;
			}
		}
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				if (enc->weight[y][x][cl] == 0) continue;
				u = enc->upara[y][x] + 1;
				peak = set_mask_parameter_optimize(enc, y, x, u-1,cl);
				e = enc->encval[y][x];
				prd = enc->prd_class[y][x][cl];
				prd = CLIP(0, enc->maxprd, prd);
				frac = enc->fconv[prd];
				for (gr = 0; gr < enc->num_group; gr++) {
					mask->pm[peak] = enc->pmlist[gr] + frac;
					if (mask->num_peak == 1){
						base = mask->base[0];
						pm = mask->pm[0];
						cbuf[gr][u] += pm->cost[base + e] + pm->subcost[base];
					}else{
						set_pmodel_mult_cost(mask,enc->maxval+1,e);
						cbuf[gr][u] += a * (log(mask->cumfreq)-log(mask->freq));
					}
				}
			}
		}

		for (gr = 0; gr < enc->num_group; gr++) {
			cbuf_p = cbuf[gr];
			for (u = 1; u < MAX_UPARA + 2; u++) {
				cbuf_p[u] += cbuf_p[u - 1];
			}
		}
		cbuf_p = cbuf[0];
		for (u = 0; u < MAX_UPARA + 2; u++) {
			dpcost[u] = cbuf_p[u] + thc_p[u];
		}
		for (gr = 1; gr < enc->num_group - 1; gr++) {
			cbuf_p = cbuf[gr];
			/* minimize (cbuf_p[th1] - cbuf_p[th0] + dpcost[th0]) */
			for (th1 = MAX_UPARA + 1; th1 >= 0; th1--) {
				th0 = th1;
				min_cost = dpcost[th1] - cbuf_p[th1] + thc_p[0];
				for (k = 0; k < th1; k++) {
					cost = dpcost[k] - cbuf_p[k] + thc_p[th1 - k];
					if (cost < min_cost) {
						min_cost = cost;
						th0 = k;
					}
				}
				dpcost[th1] = min_cost + cbuf_p[th1];
				trellis[gr][th1] = th0;
			}
		}

		cbuf_p = cbuf[gr];
		/* minimize (cbuf_p[th1] - cbuf_p[th0] + dpcost[th0]) for last group */
		th1 = MAX_UPARA + 1;
		th0 = th1;
		min_cost = dpcost[th1] - cbuf_p[th1];
		for (k = 0; k < th1; k++) {
			cost = dpcost[k] - cbuf_p[k];
			if (cost < min_cost) {
				min_cost = cost;
				th0 = k;
			}
		}
		trellis[gr][th1] = th0;

		for (gr = enc->num_group - 1; gr > 0; gr--) {
			th1 = trellis[gr][th1];
			enc->th[cl][gr - 1] = th1;
		}
	/*** set***/
		u = 0;
		for (gr = 0; gr < enc->num_group; gr++) {
			for (; u < enc->th[cl][gr]; u++) {
				enc->uquant[cl][u] = gr;
			}
		}
	/*** *****/
	}

	/* set context quantizer */
	for (cl = 0; cl < enc->num_class; cl++) {
		u = 0;
		for (gr = 0; gr < enc->num_group; gr++) {
			for (; u < enc->th[cl][gr]; u++) {
				enc->uquant[cl][u] = gr;
			}
		}
	}

	/* renew groups */
	cost = 0;
	pm_p = enc->pmlist;
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			cl = enc->class[y][x];
			u = enc->upara[y][x];
			set_mask_parameter(enc,y,x,u);
			enc->group[y][x] = gr = enc->uquant[cl][u];
			e = enc->encval[y][x];
			if (mask->num_peak == 1){
				base = mask->base[0];
				pm = mask->pm[0];
				cost += pm->cost[base + e] + pm->subcost[base];
			}else{
				set_pmodel_mult_cost(mask,enc->maxval+1,e);
				cost += a * (log(mask->cumfreq)-log(mask->freq));
			}
		}
	}
printf ("[op_group]-> %d" ,(int)cost);	//しきい値毎に分散を最適化した時のコスト算出

	/* optimize probability models */
	if (enc->optimize_loop > 1 && enc->num_pmodel > 1) {
		if (enc->num_pmodel > MAX_UPARA + 2) {
			free(cbuf);
			cbuf = (cost_t **)alloc_2d_array(enc->num_group, enc->num_pmodel,
				sizeof(cost_t));
		}

		for (gr = 0; gr < enc->num_group; gr++) {
			for (k = 0; k < enc->num_pmodel; k++) {
				cbuf[gr][k] = 0;
			}
		}
		for (gr = 0; gr < enc->num_group; gr++){
			for (y = 0; y < enc->height; y++) {
				for (x = 0; x < enc->width; x++) {
					count = 0;//init
					e = enc->encval[y][x];
					u = enc->upara[y][x];
					set_mask_parameter_optimize2(enc, y, x, u);
					for (k = 0; k < enc->num_pmodel; k++) {
						for(peak = 0; peak < mask->num_peak; peak++){
							cl = mask->class[peak];
							prd = enc->prd_class[y][x][cl];
							prd = CLIP(0, enc->maxprd, prd);
							mask->base[peak] = enc->bconv[prd];
							frac = enc->fconv[prd];
							m_gr = enc->uquant[cl][u];
							if (m_gr == gr)	 {
								mask->pm[peak] =  enc->pmodels[gr][k] + frac;
								count++;
							}else{
								mask->pm[peak] = enc->pmlist[m_gr] + frac;
							}
						}
						if (count == 0) break;
						if (mask->num_peak == 1){
							base = mask->base[0];
							pm = mask->pm[0];
							cbuf[gr][k] += pm->cost[base + e] + pm->subcost[base];
						}else{
							set_pmodel_mult_cost(mask,enc->maxval+1,e);
							cbuf[gr][k] += a * (log(mask->cumfreq)-log(mask->freq));
						}
					}
				}
			}
			pm = enc->pmodels[gr][0];
			cost = cbuf[gr][0];
			for (k = 1; k < enc->num_pmodel; k++) {
				if (cost > cbuf[gr][k]) {
					cost = cbuf[gr][k];
					pm = enc->pmodels[gr][k];
				}
			}
			pm_p[gr] = pm;
		}
/*
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				gr = enc->group[y][x];
				e = enc->encval[y][x];
				prd = enc->prd[y][x];
				prd = CLIP(0, enc->maxprd, prd);
				base = enc->bconv[prd];
				frac = enc->fconv[prd];
				for (k = 0; k < enc->num_pmodel; k++) {
					pm = enc->pmodels[gr][k] + frac;
					cbuf[gr][k] += pm->cost[base + e] + pm->subcost[base];
				}
			}
		}
		for (gr = 0; gr < enc->num_group; gr++) {
			pm = enc->pmodels[gr][0];
			cost = cbuf[gr][0];
			for (k = 1; k < enc->num_pmodel; k++) {
				if (cost > cbuf[gr][k]) {
					cost = cbuf[gr][k];
					pm = enc->pmodels[gr][k];
				}
			}
			pm_p[gr] = pm;
		}
*/
		cost = 0.0;
/*
		for (gr = 0; gr < enc->num_group; gr++) {
			cost += cbuf[gr][pm_p[gr]->id];
		}
	}
*/
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				cl = enc->class[y][x];
				u = enc->upara[y][x];
				set_mask_parameter(enc,y,x,u);
				enc->group[y][x] = gr = enc->uquant[cl][u];
				e = enc->encval[y][x];
				if (mask->num_peak == 1){
					base = mask->base[0];
					pm = mask->pm[0];
					cost += pm->cost[base + e] + pm->subcost[base];
				}else{
					set_pmodel_mult_cost(mask,enc->maxval+1,e);
					cost += a * (log(mask->cumfreq)-log(mask->freq));
				}
			}
		}
	}
	printf ("[op_c]->");

	free(cbuf);
	free(dpcost);
	free(trellis);
	return (cost);
}

#endif

int putbits(FILE *fp, int n, uint x)
{
	static int bitpos = 8;
	static uint bitbuf = 0;
	int bits;

	bits = n;
	if (bits <= 0) return (0);
	while (n >= bitpos) {
		n -= bitpos;
		if (n < 32) {
			bitbuf |= ((x >> n) & (0xff >> (8 - bitpos)));
		}
		putc(bitbuf, fp);
		bitbuf = 0;
		bitpos = 8;
	}
	bitpos -= n;
	bitbuf |= ((x & (0xff >> (8 - n))) << bitpos);
	return (bits);
}

void remove_emptyclass(ENCODER *enc)
{
	int cl, i, k, x, y;

	for (cl = 0; cl < enc->num_class; cl++) {
		enc->mtfbuf[cl] = 0;
	}
	for (y = 0; y < enc->height; y += MIN_BSIZE) {
		for (x = 0; x < enc->width; x += MIN_BSIZE) {
			cl = enc->class[y][x];
			enc->mtfbuf[cl]++;
		}
	}
	for (i = cl = 0; i < enc->num_class; i++) {
		if (enc->mtfbuf[i] == 0) {
			enc->mtfbuf[i] = -1;
		} else {
			enc->mtfbuf[i] = cl++;	//クラス番号の要素に新しいクラス番号が入る
		}
	}
	if (cl == enc->num_class) return;	/* no empty class */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			i = enc->class[y][x];
			enc->class[y][x] = enc->mtfbuf[i];
		}
	}
	for (i = cl = 0; i < enc->num_class; i++) {
		if (enc->mtfbuf[i] < 0) continue;
		if (cl != i) {
			for (k = 0; k < enc->max_prd_order; k++) {
				enc->predictor[cl][k] = enc->predictor[i][k];
			}
			for (k = 0; k < enc->num_group - 1; k++) {
				enc->th[cl][k] = enc->th[i][k];
			}
#if MULT_PEAK_MODE
			for (k = 0; k < MAX_UPARA+1; k++) {
				enc->uquant[cl][k] = enc->uquant[i][k];
			}
			for (y = 0; y < enc->height; y++) {
				for (x = 0; x < enc->width; x++) {
					enc->prd_class[y][x][cl] = enc->prd_class[y][x][i];
				}
			}
#endif
		}
		cl++;
	}
	printf("M = %d\n", cl);
	enc->num_class = cl;
	predict_region(enc, 0, 0, enc->height, enc->width);
}


int write_header(ENCODER *enc, FILE *fp)
{
	int bits;

	bits = putbits(fp, 16, MAGIC_NUMBER);
	bits += putbits(fp, 8, VERSION);
	bits += putbits(fp, 16, enc->width);
	bits += putbits(fp, 16, enc->height);
	bits += putbits(fp, 16, enc->maxval);
	bits += putbits(fp, 6, enc->num_class);
	bits += putbits(fp, 6, enc->num_group);
	bits += putbits(fp, 7, enc->max_prd_order);
	bits += putbits(fp, 6, enc->num_pmodel - 1);
	bits += putbits(fp, 4, enc->coef_precision - 1);
	bits += putbits(fp, 3, enc->pm_accuracy);
	bits += putbits(fp, 1, (enc->quadtree_depth < 0)? 0 : 1);

#if TEMPLATE_MATCHING_ON
	// bits += putbits(fp, 6, enc->temp_cl);
	// printf("TEMP_CL : %d | ", enc->temp_cl);
	bits += putbits(fp, 6, enc->temp_peak_num);
	printf("TEMP_PEAK_NUM: %d\n", enc->temp_peak_num);
#endif

	return (bits);
}

void set_qtindex(ENCODER *enc, int *index, uint *hist, int *numidx,
				 int tly, int tlx, int blksize, int width, int level)
{
	int i, cl, x, y, ctx;
	char **qtmap;

	if (tly >= enc->height || tlx >= enc->width) return;
	if (level > 0) {
		/* context modeling for quad-tree flag */
		ctx = 0;
		qtmap = enc->qtmap[level - 1];
		y = ((tly + MIN_BSIZE - 1) / MIN_BSIZE) >> level;
		x = ((tlx + MIN_BSIZE - 1) / MIN_BSIZE) >> level;
		if (y > 0) {
			if (qtmap[y - 1][x] == 1) ctx++;
			if (tlx + blksize < width && qtmap[y - 1][x + 1] == 1) ctx++;
		}
		if (x > 0 && qtmap[y][x - 1] == 1) ctx++;
		ctx = ((level - 1) * 4 + ctx) << 1;
		if (qtmap[y][x] == 1) {
			ctx++;
			index[(*numidx)++] = -(ctx + 1);
			enc->qtctx[ctx]++;
			blksize >>= 1;
			set_qtindex(enc, index, hist, numidx, tly, tlx,
				blksize, width, level - 1);
			set_qtindex(enc, index, hist, numidx, tly, tlx + blksize,
				blksize, width, level - 1);
			set_qtindex(enc, index, hist, numidx, tly + blksize, tlx,
				blksize, width, level - 1);
			width = tlx + blksize * 2;
			if (width >= enc->width) width = enc->width;
			set_qtindex(enc, index, hist, numidx, tly + blksize, tlx + blksize,
				blksize, width, level - 1);
			return;
		} else {
			index[(*numidx)++] = -(ctx + 1);
			enc->qtctx[ctx]++;
		}
	}
	cl = enc->class[tly][tlx];
	mtf_classlabel(enc->class, enc->mtfbuf, tly, tlx, blksize, width, enc->num_class);
	i = enc->mtfbuf[cl];
	index[(*numidx)++] = i;
	hist[i]++;
	return;
}

int encode_class(FILE *fp, ENCODER *enc, int flag)
{
	int i, j, k, numidx, blksize, level, x, y, ctx, bits, *index;
	uint *hist;
	cost_t cost;
	PMODEL *pm;
	double p, c;
	int qtree_code[QUADTREE_DEPTH << 2], mtf_code[MAX_CLASS];
	cost_t qtflag_cost[QUADTREE_DEPTH << 3], *class_cost;

#if (!OPT_SIDEINFO)
	if (fp == NULL) return(0);
#endif
	if (enc->quadtree_depth >= 0 && enc->optimize_loop > 1) {
		level = enc->quadtree_depth;
		blksize = MAX_BSIZE;
		numidx = 0;
		y = (enc->height + MIN_BSIZE - 1) / MIN_BSIZE;
		x = (enc->width + MIN_BSIZE - 1) / MIN_BSIZE;
		for (k = 0; k <= level; k++) {
			numidx += y * x;
			y = (y + 1) >> 1;
			x = (x + 1) >> 1;
		}
		for (k = 0; k < QUADTREE_DEPTH << 3; k++) {
			enc->qtctx[k] = 0;
		}
	} else {
		level = 0;
		blksize = BASE_BSIZE;
		numidx = ((enc->height + (BASE_BSIZE - 1)) / BASE_BSIZE)
			* ((enc->width + (BASE_BSIZE - 1)) / BASE_BSIZE);
	}
	hist = (uint *)alloc_mem(enc->num_class * sizeof(uint));
	index = (int *)alloc_mem(numidx * sizeof(int));
	for (i = 0; i < enc->num_class; i++) {
		hist[i] = 0;
		enc->mtfbuf[i] = i;
	}
	numidx = 0;
	for (y = 0; y < enc->height; y += blksize) {
		for (x = 0; x < enc->width; x += blksize) {
			set_qtindex(enc, index, hist, &numidx, y, x,
				blksize, enc->width, level);
		}
	}
	bits = 0;
	/* Arithmetic */
	class_cost = (cost_t *)alloc_mem(enc->num_class * sizeof(cost_t));
	/* context modeling for quad-tree flag */
	if (level > 0) {
		for (ctx = 0; ctx < QUADTREE_DEPTH << 2; ctx++) {
			cost = INT_MAX;
			for (i = k = 0; i < 7; i++) {
				p = qtree_prob[i];
				c = -log(p) * (cost_t)enc->qtctx[(ctx << 1) + 1]
				-log(1.0 - p) * (cost_t)enc->qtctx[ctx << 1];
				if (c < cost) {
					k = i;
					cost = c;
				}
			}
			p = qtree_prob[qtree_code[ctx] = k];
			qtflag_cost[(ctx << 1) + 1] = -log(p) / log(2.0);
			qtflag_cost[ctx << 1] = -log(1.0 - p) / log(2.0);
		}
	}
	/* quantization of log-transformed probability */
	c = 0.0;
	for (i = 0; i < enc->num_class; i++) {
		c += (double)hist[i];
	}
	for (i = 0; i < enc->num_class; i++) {
		p = (double)hist[i] / c;
		if (p > 0.0) {
			mtf_code[i] = (int)(-log(p) / log(2.0) * (PMCLASS_LEVEL / PMCLASS_MAX));
			if (mtf_code[i] >= PMCLASS_LEVEL) {
				mtf_code[i] = PMCLASS_LEVEL - 1;
			}
		} else {
			mtf_code[i] = PMCLASS_LEVEL - 1;
		}
		p = exp(-log(2.0) * ((double)mtf_code[i] + 0.5)
			* PMCLASS_MAX / PMCLASS_LEVEL);
		class_cost[i] = -log(p) / log(2.0);
		hist[i] = (uint)(p * (1 << 10));
		if (hist[i] <= 0) hist[i] = 1;
	}
	if (fp == NULL) {
		cost = 0.0;
		for (j = 0; j < numidx; j++) {
			i = index[j];
			if (i < 0) {
				ctx = -(i + 1);
				cost += qtflag_cost[ctx];
			} else {
				cost += class_cost[i];
			}
		}
		if (flag == 1) {
			for (i = 0; i < (QUADTREE_DEPTH << 3); i++)
				enc->qtflag_cost[i] = qtflag_cost[i];
			for (i = 0; i < enc->num_class; i++)
				enc->class_cost[i] = class_cost[i];
		}
		bits = (int)cost;
	} else {	/* actually encode */
		PMODEL cpm[1];
		/* additional info. */
		pm = &enc->spm;
		if (level > 0) {
			set_spmodel(pm, 7, -1);
			for (ctx = 0; ctx < QUADTREE_DEPTH << 2; ctx++) {
				i = qtree_code[ctx];
				rc_encode(fp, enc->rc, pm->cumfreq[i], pm->freq[i],
					pm->cumfreq[pm->size]);
			}
		}
		set_spmodel(pm, PMCLASS_LEVEL, -1);
		for (i = 0; i < enc->num_class; i++) {
			j = mtf_code[i];
			rc_encode(fp, enc->rc, pm->cumfreq[j], pm->freq[j],
				pm->cumfreq[pm->size]);
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
		cpm->size = enc->num_class;
		cpm->freq = (uint *)alloc_mem((cpm->size * 2 + 1) * sizeof(uint));
		cpm->cumfreq = &(cpm->freq[cpm->size]);
		cpm->cumfreq[0] = 0;
		for (i = 0; i < enc->num_class; i++) {
			cpm->freq[i] = hist[i];
			cpm->cumfreq[i + 1] = cpm->cumfreq[i] + cpm->freq[i];
		}
		for (j = 0; j < numidx; j++) {
			i = index[j];
			if (i < 0) {
				i = -(i + 1);
				ctx = i & (~1);
				rc_encode(fp, enc->rc,
					pm->cumfreq[i] - pm->cumfreq[ctx], pm->freq[i],
					pm->cumfreq[ctx + 2] - pm->cumfreq[ctx]);
			} else {
				rc_encode(fp, enc->rc, cpm->cumfreq[i], cpm->freq[i],
					cpm->cumfreq[cpm->size]);
			}
		}
		bits += (int)enc->rc->code;
		enc->rc->code = 0;
	}
	free(index);
	free(hist);
	return (bits);
}

#if AUTO_PRD_ORDER

int encode_predictor(FILE *fp, ENCODER *enc, int flag)	//when AUTO_PRD_ORDER 1
{
	int cl, coef, sgn, k, m, min_m, bits, d;
	cost_t cost, min_cost, t_cost;
	PMODEL *pm;
	uint cumb;
	int zrfreq, nzfreq;

#if (!OPT_SIDEINFO)
	if (fp == NULL) return(0);
#endif
	t_cost = 0.0;
	for (d = 0; d < enc->prd_mhd; d++) {
		min_cost = INT_MAX;
		for (m = min_m = 0; m < 8; m++) {
			cost = 0.0;
			for (k = d * (d + 1); k < (d + 1) * (d + 2); k++) {
				for (cl = 0; cl < enc->num_class; cl++) {
				#if TEMPLATE_MATCHING_ON
					if(cl == enc->temp_cl) continue;
				#endif
					coef = enc->predictor[cl][k];
					if (coef < 0) coef = -coef;
					cost += enc->coef_cost[enc->zero_m[d]][m][coef];
				}
			}
			if (cost < min_cost) {
				min_cost = cost;
				min_m = m;
			}
		}
		if (flag) enc->coef_m[d] = min_m;
		min_cost = INT_MAX;
		for (m = min_m = 0; m < NUM_ZMODEL; m++) {
			cost = 0.0;
			for (k = d * (d + 1); k < (d + 1) * (d + 2); k++) {
				for (cl = 0; cl < enc->num_class; cl++) {
				#if TEMPLATE_MATCHING_ON
					if(cl == enc->temp_cl) continue;
				#endif
					coef = enc->predictor[cl][k];
					if (coef < 0) coef = -coef;
					cost += enc->coef_cost[m][enc->coef_m[d]][coef];
				}
			}
			if (cost < min_cost) {
				min_cost = cost;
				min_m = m;
			}
		}
		if (flag) enc->zero_m[d] = min_m;
		t_cost += min_cost;
	}
	bits = (int)t_cost;
#if CHECK_PREDICTOR
	if(fp != NULL){
		for(cl=0; cl<enc->num_class; cl++){
			if(cl == enc->temp_cl){
				printf("[%2d] TEMPLATE_CLASS\n", cl);
			} else {
				printf("[%2d]", cl);
				for(k=0; k<enc->max_prd_order; k++){
					printf("%d ", enc->predictor[cl][k]);
				}
				printf("\n");
			}
		}
	}
#endif
	/* Arithmetic */
	if (fp != NULL) {
		bits = 0;
		// PMODEL *pm;
		// uint cumb;
		// int zrfreq, nzfreq;
		pm = &enc->spm;
		for (d = 0; d < enc->prd_mhd; d++) {
			rc_encode(fp, enc->rc, enc->zero_m[d], 1, NUM_ZMODEL);
			rc_encode(fp, enc->rc, enc->coef_m[d], 1, 8);
			nzfreq = enc->zero_fr[enc->zero_m[d]];
			zrfreq = TOT_ZEROFR - nzfreq;
			set_spmodel(pm, enc->max_coef + 1, enc->coef_m[d]);
			cumb = pm->freq[0];
			for (k = d * (d + 1); k < (d + 1) * (d + 2); k++) {
				for (cl = 0; cl < enc->num_class; cl++) {
				#if TEMPLATE_MATCHING_ON
					if(cl == enc->temp_cl)continue;	//clがテンプレートマッチングをするクラスならcontinue
				#endif
					coef = enc->predictor[cl][k];
					if (coef == 0) {
						rc_encode(fp, enc->rc, 0, zrfreq, TOT_ZEROFR);
					} else {
						rc_encode(fp, enc->rc, zrfreq, nzfreq, TOT_ZEROFR);
						sgn = (coef < 0)? 1 : 0;
						if (coef < 0) coef = -coef;
						rc_encode(fp, enc->rc, pm->cumfreq[coef] - cumb,  pm->freq[coef],
							pm->cumfreq[pm->size] - cumb);
						rc_encode(fp, enc->rc, sgn, 1, 2);
					}
				}
			}
		}
		bits = (int)enc->rc->code;
		enc->rc->code = 0;
	}
	return (bits);
}

#else

/* change pmodel each position */
int encode_predictor(FILE *fp, ENCODER *enc, int flag)	//when AUTO_PRD_ORDER 0
{
	int cl, coef, sgn, k, m, min_m, bits;
	cost_t cost, min_cost, t_cost;
	PMODEL *pm;
	pm = &enc->spm;

#if (!OPT_SIDEINFO)
	if (fp == NULL) return(0);
#endif
	t_cost = 0.0;
	for (k = 0; k < enc->prd_order; k++) {
		min_cost = INT_MAX;
		for (m = min_m = 0; m < 16; m++) {
			cost = 0.0;
			for (cl = 0; cl < enc->num_class; cl++) {
				coef = enc->predictor[cl][k];
				if (coef < 0) coef = -coef;
				cost += enc->coef_cost[m][coef];
			}
			if (cost < min_cost) {
				min_cost = cost;
				min_m = m;
			}
		}
		t_cost += min_cost;
		if (flag) enc->coef_m[k] = min_m;
	}
	bits = (int)t_cost;
	if (fp != NULL) {
		bits = 0;
		for (k = 0; k < enc->prd_order; k++) {
			set_spmodel(pm, enc->max_coef + 1, enc->coef_m[k]);
			rc_encode(fp, enc->rc, enc->coef_m[k], 1, 16);
			for (cl = 0; cl < enc->num_class; cl++) {
				coef = enc->predictor[cl][k];
				sgn = (coef < 0)? 1 : 0;
				if (coef < 0) coef = -coef;
				rc_encode(fp, enc->rc, pm->cumfreq[coef],  pm->freq[coef],
					pm->cumfreq[pm->size]);
				if (coef > 0) {
					rc_encode(fp, enc->rc, sgn, 1, 2);
				}
			}
		}
		bits = (int)enc->rc->code;
		enc->rc->code = 0;
	}
	return (bits);
}

#endif

int encode_threshold(FILE *fp, ENCODER *enc, int flag)
{
	int cl, gr, i, k, m, min_m, bits;
	cost_t cost, min_cost;
	PMODEL *pm;
	double p;
	pm = &enc->spm;

#if (!OPT_SIDEINFO)
	if (fp == NULL) return(0);
#endif
	/* Arithmetic */
	min_cost = INT_MAX;
	for (m = min_m = 0; m < 16; m++) {
		set_spmodel(pm, MAX_UPARA + 2, m);
		cost = 0.0;
		for (cl = 0; cl < enc->num_class; cl++) {
			k = 0;
			for (gr = 1; gr < enc->num_group; gr++) {
				i = enc->th[cl][gr - 1] - k;
				p = (double)pm->freq[i] / (pm->cumfreq[pm->size - k]);
				cost += -log(p);
				k += i;
				if (k > MAX_UPARA) break;
			}
		}
		cost /= log(2.0);
		if (cost < min_cost) {
			min_cost = cost;
			min_m = m;
		}
	}
	set_spmodel(pm, MAX_UPARA + 2, min_m);
	p = log(pm->cumfreq[MAX_UPARA + 2]);
	if (fp == NULL) {
		if (flag == 1){
			for (i = 0; i < MAX_UPARA + 2; i++) {
				enc->th_cost[i] = (p - log(pm->freq[i])) / log(2.0);
			}
		}
		bits = (int)min_cost;
	} else {
		rc_encode(fp, enc->rc, min_m, 1, 16);
		for (cl = 0; cl < enc->num_class; cl++) {
			k = 0;
			for (gr = 1; gr < enc->num_group; gr++) {
				i = enc->th[cl][gr - 1] - k;
				rc_encode(fp, enc->rc, pm->cumfreq[i],  pm->freq[i],
					pm->cumfreq[pm->size - k]);
				k += i;
				if (k > MAX_UPARA) break;
			}
		}
		if (enc->num_pmodel > 1) {
			for (gr = 0; gr < enc->num_group; gr++) {
				pm = enc->pmlist[gr];
				rc_encode(fp, enc->rc, pm->id, 1, enc->num_pmodel);
			}
		}
		bits = (int)enc->rc->code;
		enc->rc->code = 0;
	}
	return (bits);
}

#if MULT_PEAK_MODE
//***************************************************
#if 0
int encode_mask(FILE *fp, ENCODER *enc)
{
	int tly, tlx, bry, brx, cl, y, x, yy, xx, bits;

	for(y = 0; y < enc->height; y += WIN_BSIZE){
		for(x = 0; x < enc->width; x += WIN_BSIZE){
			if (enc->mask[y][x] == 0){
				cl = enc->class[y][x];
				tly = y - (NUM_MASK-1);
				tlx = x - (NUM_MASK-1);
				bry = y + WIN_BSIZE + (NUM_MASK-1);
				brx = x + WIN_BSIZE + (NUM_MASK-1);
				if(tly < 0) tly = 0;
				if(tlx < 0) tlx = 0;
				if(bry > enc->height) bry = enc->height;
				if(brx > enc->width) brx = enc->width;
				for(yy = tly; yy < bry; yy++){
					for(xx = tlx; xx < brx; xx++){
						if (enc->class[yy][xx] != cl) goto enc;
					}
				}
//				printf("skip(%d,%d)",y,x);
				continue;
			}
enc:
			rc_encode(fp, enc->rc, enc->mask[y][x], 1, NUM_MASK);
		}
	}
	bits = enc->rc->code;
	enc->rc->code = 0;
	return(bits);
}
#endif

int encode_mask(FILE *fp, ENCODER *enc, int flag)
{
	int i, j, k, numidx, blksize, x, y, bits, *index;
	uint *hist;
	cost_t cost;
	PMODEL *pm;
	double p, c;
	int mtf_code[NUM_MASK];
	cost_t *mask_cost;

#if (!OPT_SIDEINFO)
	if (fp == NULL) return(0);
#endif

	blksize = WIN_BSIZE;
	numidx = ((enc->height + (WIN_BSIZE - 1)) / WIN_BSIZE)
			* ((enc->width + (WIN_BSIZE - 1)) / WIN_BSIZE);

	hist = (uint *)alloc_mem(NUM_MASK * sizeof(uint));
	index = (int *)alloc_mem(numidx * sizeof(int));

	for (i = 0; i < NUM_MASK; i++) {
		hist[i] = 0;
		enc->mtfbuf[i] = i;
	}
	numidx = 0;
	for (y = 0; y < enc->height; y += blksize) {
		for (x = 0; x < enc->width; x += blksize) {
			// if (y >= enc->height || x >= enc->width) continue;
			k = enc->mask[y][x];
			mtf_classlabel(enc->mask, enc->mtfbuf, y, x, blksize, enc->width, NUM_MASK);
			i = enc->mtfbuf[k];
			index[numidx++] = i;
			hist[i]++;
		}
	}

	bits = 0;
	/* Arithmetic */
	mask_cost = (cost_t *)alloc_mem(NUM_MASK * sizeof(cost_t));
	/* context modeling for quad-tree flag */

	/* quantization of log-transformed probability */
	c = 0.0;
	for (i = 0; i < NUM_MASK; i++) {
		c += (double)hist[i];
	}
	for (i = 0; i < NUM_MASK; i++) {
		p = (double)hist[i] / c;
		if (p > 0.0) {
			mtf_code[i] = (int)(-log(p) / log(2.0) * (PMMASK_LEVEL / PMMASK_MAX));
			if (mtf_code[i] >= PMMASK_LEVEL) {
				mtf_code[i] = PMMASK_LEVEL - 1;
			}
		} else {
			mtf_code[i] = PMMASK_LEVEL - 1;
		}
		p = exp(-log(2.0) * ((double)mtf_code[i] + 0.5)
			* PMMASK_MAX / PMMASK_LEVEL);
		mask_cost[i] = -log(p) / log(2.0);
		hist[i] = (uint)(p * (1 << 10));
		if (hist[i] <= 0) hist[i] = 1;
	}
	if (fp == NULL) {
		cost = 0.0;
		for (j = 0; j < numidx; j++) {
			i = index[j];
				cost += mask_cost[i];
		}
//		if (flag == 1) {
//			for (i = 0; i < NUM_MASK; i++)
//				enc->mask_cost[i] = mask_cost[i];
//		}
		bits = (int)cost;
	} else {	/* actually encode */
		PMODEL cpm[1];
		/* additional info. */
		pm = &enc->spm;

		set_spmodel(pm, PMMASK_LEVEL, -1);
		for (i = 0; i < NUM_MASK; i++) {
			j = mtf_code[i];
			rc_encode(fp, enc->rc, pm->cumfreq[j], pm->freq[j],
				pm->cumfreq[pm->size]);
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
		cpm->cumfreq = &(cpm->freq[cpm->size]);
		cpm->cumfreq[0] = 0;
		for (i = 0; i < NUM_MASK; i++) {
			cpm->freq[i] = hist[i];
			cpm->cumfreq[i + 1] = cpm->cumfreq[i] + cpm->freq[i];
		}
		for (j = 0; j < numidx; j++) {
			i = index[j];
				rc_encode(fp, enc->rc, cpm->cumfreq[i], cpm->freq[i],
				cpm->cumfreq[cpm->size]);
		}
		bits += (int)enc->rc->code;
		enc->rc->code = 0;
	}
	free(index);
	free(hist);
	return (bits);
}

cost_t calc_side_info(ENCODER *enc, cost_t cost){
	cost_t sc;

	#if CHECK_DEBUG
		printf("(%d,", (int)cost);
	#endif

	cost += sc = encode_predictor(NULL, enc, 1);
	#if CHECK_DEBUG
		printf("%d(%d),", (int)sc, enc->num_class);
	#endif

	cost += sc = encode_threshold(NULL, enc, 1);
	#if CHECK_DEBUG
		printf("%d,", (int)sc);
	#endif

#if TEMPLATE_MATCHING_ON
	cost += sc = encode_w_gr_threshold(NULL, enc, 1);
	#if CHECK_DEBUG
		printf("%d,", (int)sc);
	#endif
#endif
	cost += sc = encode_class(NULL, enc, 1);
	#if CHECK_DEBUG
		printf("%d,", (int)sc);
	#endif

	#if CHECK_DEBUG
		printf("| %d)", (int)cost);
	#endif
	return(cost);
}


int encode_image(FILE *fp, ENCODER *enc)	//多峰性確率モデル
{
	int x, y, e, u, base, bits, cumbase;
	PMODEL *pm;
	enc->function_number = 8;
	bits = 0;
	/* Arithmetic */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			// u = enc->upara[y][x];
			u = calc_uenc(enc, y, x);
			set_mask_parameter(enc,y,x,u);
			e = enc->encval[y][x];
			if (mask->num_peak == 1){
				base = mask->base[0];
				pm = mask->pm[0];
				cumbase = pm->cumfreq[base];
				enc->rc->y = y;
				enc->rc->x = x;
				rc_encode(fp, enc->rc,
					pm->cumfreq[base + e] - cumbase,
					pm->freq[base + e],
					pm->cumfreq[base + enc->maxval + 1] - cumbase);
			}else{
				pm = &enc->mult_pm;
				set_pmodel_mult(pm,mask,enc->maxval+1);
				#if CHECK_PMODEL
					if(y==check_y && x==check_x) printmodel(pm,enc->maxval+1);
				#endif
				enc->rc->y = y;
				enc->rc->x = x;
				rc_encode(fp, enc->rc
,					pm->cumfreq[e],
					pm->freq[e],
					pm->cumfreq[enc->maxval + 1]);
			}
			// printf("%d,%d,prd,%d(%d),org,%d,conv:%d\n", y, x, enc->prd_class[y][x][enc->class[y][x]] >> (enc->coef_precision-1), enc->class[y][x], enc->org[y][x], enc->err[y][x]);
		}
	}
	rc_finishenc(fp, enc->rc);
	bits += (int)enc->rc->code;
	return (bits);
}

#else
int encode_image(FILE *fp, ENCODER *enc)	//単峰性確率モデル
{
	int x, y, e, prd, base, bits, gr, cumbase;
	PMODEL *pm;

	bits = 0;
	/* Arithmetic */
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			gr = enc->group[y][x];
			prd = enc->prd[y][x];
			prd = CLIP(0, enc->maxprd, prd);
			e = enc->encval[y][x];
			base = enc->bconv[prd];
			pm = enc->pmlist[gr] + enc->fconv[prd];
			cumbase = pm->cumfreq[base];
			rc_encode(fp, enc->rc,
				pm->cumfreq[base + e] - cumbase,
				pm->freq[base + e],
				pm->cumfreq[base + enc->maxval + 1] - cumbase);
		}
	}
	rc_finishenc(fp, enc->rc);
	bits += (int)enc->rc->code;
	return (bits);
}
#endif

#if AUTO_DEL_CL
int opcl_find_class(ENCODER *enc, int k, int tly, int tlx, int bry, int brx, int bufsize)
{
	cost_t cost, min_cost;
	int x, y, cl, min_cl, prd, org;
	char *class_p;
	int *prd_p, *err_p, *org_p;
//	int *prd_p, *prdbuf_p, *err_p, *errbuf_p;
	min_cost = 1E8;
	min_cl = 0;

	for	(cl = 0; cl < enc->num_class; cl++) {
		if (cl == k) continue;
//		bufptr = bufsize * (tly % bufsize) + tlx % bufsize;
		for (y = tly; y < bry; y++) {
			class_p = &enc->class[y][tlx];
			prd_p = &enc->prd[y][tlx];
//			prdbuf_p = &prdbuf[cl][bufptr];
			err_p = &enc->err[y][tlx];
//			errbuf_p = &errbuf[cl][bufptr];
//			bufptr += bufsize;
			org_p = &enc->org[y][tlx];
			for (x = tlx; x < brx; x++) {
				*class_p++ = cl;
				org = *org_p++;
				prd = enc->prd_class[y][x][cl];
				*prd_p++ = prd;
  			prd = CLIP(0, enc->maxprd, prd);
	  		prd >>= (enc->coef_precision - 1);
  			*err_p++ = enc->econv[org][prd];
//				*prd_p++ = *prdbuf_p++;
//				*err_p++ = *errbuf_p++;
			}
		}
		cost = calc_cost2(enc, tly, tlx, bry, brx);
#if OPTIMIZE_MASK_LOOP
			cost += rap_area_cost(enc, tly, tlx, bry, brx);
#endif
			cost += enc->class_cost[enc->mtfbuf[cl]];
		if (cost < min_cost) {
			min_cost = cost;
			min_cl = cl;
		}
	}

//	bufptr = bufsize * (tly % bufsize) + tlx % bufsize;
	for (y = tly; y < bry; y++) {
		class_p = &enc->class[y][tlx];
		prd_p = &enc->prd[y][tlx];
//		prdbuf_p = &prdbuf[min_cl][bufptr];
		err_p = &enc->err[y][tlx];
//		errbuf_p = &errbuf[min_cl][bufptr];
//		bufptr += bufsize;
		org_p = &enc->org[y][tlx];
		for (x = tlx; x < brx; x++) {
			*class_p++ = min_cl;
					org = *org_p++;
				prd = enc->prd_class[y][x][min_cl];
				*prd_p++ = prd;
			prd = CLIP(0, enc->maxprd, prd);
  		prd >>= (enc->coef_precision - 1);
			*err_p++ = enc->econv[org][prd];
//			*prd_p++ = *prdbuf_p++;
//			*err_p++ = *errbuf_p++;
		}
	}
	return (min_cl);
}

void opcl_sub(ENCODER *enc, int k, int *blk, int tly, int tlx, int blksize, int width, int level)
{
//	int cl, min_cl, x, y, bry, brx;
	int min_cl, x, y, bry, brx;
	char **qtmap;
	// cost_t *err_cost, min_cost;

	brx = (tlx + blksize < enc->width) ? (tlx + blksize) : enc->width;
	bry = (tly + blksize < enc->height) ? (tly + blksize) : enc->height;
	if (tly >= bry || tlx >= brx) return;
	if (level > 0) {
		qtmap = enc->qtmap[level - 1];
		y = (tly / MIN_BSIZE) >> level;
		x = (tlx / MIN_BSIZE) >> level;
		if (qtmap[y][x] == 1) {
			blksize >>= 1;
			opcl_sub(enc, k, blk, tly, tlx, blksize, width, level - 1);
			opcl_sub(enc, k, blk, tly, tlx + blksize, blksize, width, level - 1);
			opcl_sub(enc, k, blk, tly + blksize, tlx, blksize, width, level - 1);
			opcl_sub(enc, k, blk, tly + blksize, tlx + blksize, blksize, brx, level - 1);
			return;
		}
	}
	mtf_classlabel(enc->class, enc->mtfbuf, tly, tlx, blksize, width, enc->num_class);
	// err_cost = enc->err_cost[*blk];
	// min_cost = 1E5;
	min_cl = 0;
#if OPTIMIZE_MASK_LOOP
	min_cl = opcl_find_class(enc, k, tly, tlx, bry, brx, BASE_BSIZE);
#else
	for (cl = 0; cl < enc->num_class; cl++) {
		if (cl == k) continue;
		cost = err_cost[cl];
		cost += enc->class_cost[enc->mtfbuf[cl]];
		if (cost < min_cost) {
			min_cost = cost;
			min_cl = cl;
		}
	}
#endif
	for (y = tly; y < bry; y++) {
		for (x = tlx; x < brx; x++) {
			enc->class[y][x] = min_cl;
		}
	}
	(*blk)++;
}

cost_t opcl(ENCODER *enc, int k, int *blk, int restore)
{
	int x, y, i, j, *th_s, *prd_s, **prd_cl_s, blksize, level;
	char *uq_s, **class;
	cost_t cost;

	class = (char **)alloc_2d_array(enc->height, enc->width, sizeof(char));
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			class[y][x] = enc->class[y][x];
		}
	}

	for (i = 0; i < enc->num_class; i++) enc->mtfbuf[i] = i;
	if( enc->optimize_loop > 1 && enc->quadtree_depth >= 0) {
		blksize = MAX_BSIZE;
	}else {
		blksize = BASE_BSIZE;
	}
	level =  enc->quadtree_depth;
	for (y = 0; y < enc->height; y += blksize) {
		for (x = 0; x < enc->width; x += blksize) {
			opcl_sub(enc, k, blk, y, x, blksize, enc->width, level);
		}
	}
	th_s = (int *)alloc_mem(enc->num_group * sizeof(int));
	prd_s = (int *)alloc_mem(enc->max_prd_order * sizeof(int));
	uq_s = (char *)alloc_mem((MAX_UPARA + 1) * sizeof(char));
	prd_cl_s = (int **)alloc_2d_array(enc->height, enc->width,
			sizeof(int));

	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
	 		if (enc->class[y][x] > k) {
				enc->class[y][x]--;
			}
			prd_cl_s[y][x] = enc->prd_class[y][x][k];
			for (i = k; i < enc->num_class - 1; i++) {
				enc->prd_class[y][x][i] = enc->prd_class[y][x][i + 1];
			}
		}
	}
	for (i = 0; i < enc->num_group; i++) {
		th_s[i] = enc->th[k][i];
	}
	for (i = 0; i < enc->max_prd_order; i++) {
		prd_s[i] = enc->predictor[k][i];
	}
	for (i = 0; i < MAX_UPARA + 1; i++) {
		uq_s[i] = enc->uquant[k][i];
	}
	for (i = k; i < enc->num_class - 1; i++) {
		for (j = 0; j < enc->num_group; j++) {
			enc->th[i][j] = enc->th[i + 1][j];
		}
		for (j = 0; j < enc->max_prd_order; j++) {
			enc->predictor[i][j] = enc->predictor[i + 1][j];
		}
		for (j = 0; j < enc->max_prd_order; j++) {
			enc->nzconv[i][j] = enc->nzconv[i + 1][j];
		}
		enc->num_nzcoef[i] = enc->num_nzcoef[i + 1];
		for (j = 0; j < MAX_UPARA + 1; j++) {
			enc->uquant[i][j] = enc->uquant[i + 1][j];
		}
	}
#if AUTO_PRD_ORDER
	set_prd_pels(enc);
#endif
	enc->num_class--;
	predict_region(enc, 0, 0, enc->height, enc->width);
	cost = calc_cost(enc, 0, 0, enc->height, enc->width);
	cost += encode_class(NULL, enc, 0);
	cost += encode_predictor(NULL, enc, 0);
	cost += encode_threshold(NULL, enc, 0);
	if (restore == 1) {
		enc->num_class++;
		for (y = 0; y < enc->height; y ++) {
			for (x = 0; x < enc->width; x ++) {
				// if (enc->class[y][x] >= k) {
					// enc->class[y][x]++;
				// }
				enc->class[y][x] = class[y][x];
				for (i = enc->num_class - 2; i >= k; i--) {
					enc->prd_class[y][x][i + 1] = enc->prd_class[y][x][i];
				}
		   	enc->prd_class[y][x][k] = prd_cl_s[y][x];
			}
		}
		for (i = enc->num_class - 2; i >= k; i--) {
			for (j = 0; j < enc->num_group; j++) {
				enc->th[i + 1][j] = enc->th[i][j];
			}
			for (j = 0; j < enc->max_prd_order; j++) {
				enc->predictor[i + 1][j] = enc->predictor[i][j];
			}
			for (j = 0; j < MAX_UPARA + 1; j++) {
				enc->uquant[i + 1][j] = enc->uquant[i][j];
			}
		}
		for (i = 0; i < enc->num_group; i++) {
			enc->th[k][i] = th_s[i];
		}
		for (i = 0; i < enc->max_prd_order; i++) {
			enc->predictor[k][i] = prd_s[i];
		}
		for (i = 0; i < MAX_UPARA + 1; i++) {
			enc->uquant[k][i] = uq_s[i];
		}
#if AUTO_PRD_ORDER
		set_prd_pels(enc);
#endif
	}
	free(th_s);
	free(prd_s);
	free(uq_s);
  	free(prd_cl_s);
	return(cost);
}

RESTORE_SIDE* init_renew_adc(ENCODER *enc){
	RESTORE_SIDE * side;
	side = (RESTORE_SIDE *)alloc_mem(sizeof(RESTORE_SIDE));
	side->th_s = (int **)alloc_2d_array(enc->num_class, enc->num_group, sizeof(int));
	// init_2d_array(side->th_s, enc->num_class, enc->num_group, 0);
	side->prd_s = (int **)alloc_2d_array(enc->num_class, enc->max_prd_order, sizeof(int));
	// init_2d_array(side->prd_s, enc->num_class, enc->max_prd_order, 0);
	side->uq_s = (char **)alloc_2d_array(enc->num_class, (MAX_UPARA + 1), sizeof(char));
	// init_2d_array(side->uq_s, enc->num_class, MAX_UPARA+1, 0);
	side->prd_cl_s = (int ***)alloc_3d_array(enc->height, enc->width, enc->num_class, sizeof(int));
	// init_3d_array(side->prd_cl_s, enc->height, enc->width, enc->num_class, 0);
	side->class = (char **)alloc_2d_array(enc->height, enc->width, sizeof(char));
	// init_2d_array(side->class, enc->height, enc->width, 0);
	return(side);
}

void save_info(ENCODER *enc, RESTORE_SIDE *r_side, int restore){
	int y, x, i, j;
	if(restore==1){
		enc->num_class = r_side->num_class_s;
		for(y=0; y<enc->height; y++){
			for(x=0; x<enc->width; x++){
				enc->class[y][x] = r_side->class[y][x];
				// r_side->class[y][x] = enc->class[y][x];
				for(i=0; i<enc->num_class; i++){
					enc->prd_class[y][x][i] = r_side->prd_cl_s[y][x][i];
					// r_side->prd_cl_s[y][x][i] = enc->prd_class[y][x][i];
				}
			}
		}

		for(i=0; i<enc->num_class; i++){
			for(j=0; j<enc->num_group; j++){
				enc->th[i][j] = r_side->th_s[i][j];
				// r_side->th_s[i][j] = enc->th[i][j];
			}
			for(j=0; j<enc->max_prd_order; j++){
				enc->predictor[i][j] = r_side->prd_s[i][j];
				// r_side->prd_s[i][j] = enc->predictor[i][j];
			}
			for(j=0; j<MAX_UPARA; j++){
				enc->uquant[i][j] = r_side->uq_s[i][j];
				// r_side->uq_s[i][j] = enc->uquant[i][j];
			}
		}
#if AUTO_PRD_ORDER
		set_prd_pels(enc);
#endif
		optimize_class(enc);
		save_prediction_value(enc);
		predict_region(enc, 0, 0, enc->height, enc->width);
	} else {
		r_side->num_class_s = enc->num_class;
		for(y=0; y<enc->height; y++){
			for(x=0; x<enc->width; x++){
				r_side->class[y][x] = enc->class[y][x];
				for(i=0; i<enc->num_class; i++){
					r_side->prd_cl_s[y][x][i] = enc->prd_class[y][x][i];
				}
			}
		}

		for(i=0; i<enc->num_class; i++){
			for(j=0; j<enc->num_group; j++){
				r_side->th_s[i][j] = enc->th[i][j];
			}
			for(j=0; j<enc->max_prd_order; j++){
				r_side->prd_s[i][j] = enc->predictor[i][j];
			}
			for(j=0; j<MAX_UPARA; j++){
				r_side->uq_s[i][j] = enc->uquant[i][j];
			}
		}
	}
}


cost_t auto_del_class(ENCODER *enc, cost_t pre_cost)
{
	int x, y, k, del_cl, blk;
	cost_t cost, min_cost;
	char **class=0;
	class = (char **)alloc_2d_array(enc->height, enc->width, sizeof(char));

	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			class[y][x] = enc->class[y][x];
		}
	}
	min_cost = 1E10;
	del_cl = 0;
	for (k = 0; k < enc->num_class; k++) {	//各クラスを削除してみて一番コストが小さくなったクラスを見つける
		if(k == enc->temp_cl)	continue;
		blk = 0;
		cost = opcl(enc, k, &blk, 1);
		if (cost < min_cost) {
			min_cost = cost;
			del_cl = k;
		}
	}
	if (pre_cost > min_cost) {
		blk = 0;
		#if CHECK_CLASS
			printf("%d ", del_cl);
		#endif
		cost = opcl(enc, del_cl, &blk, 0);
#if MULT_PEAK_MODE == 0
		for (y = del_cl; y < enc->num_class; y++) {
			for (x = 0; x < blk; x++) {
				enc->err_cost[x][y] = enc->err_cost[x][y + 1];
			}
		}
#endif
		printf("D");
		optimize_class(enc);
	} else {
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				enc->class[y][x] = class[y][x];
			}
		}
	}
	predict_region(enc, 0, 0, enc->height, enc->width);
	cost = calc_cost(enc, 0, 0, enc->height, enc->width);
	cost = calc_side_info(enc, cost);
	return(cost);
}
#endif

int main(int argc, char **argv)
{
	cost_t cost, min_cost, side_cost, sc;
	int i, j, k, l, x, y, xx, yy, cl, gr, bits, **prd_save, **th_save, sw;
	int header_info, class_info, pred_info, th_info, mask_info, err_info, w_gr_info=0;
	int num_class_save;
	char **class_save, **mask_save;
	char **qtmap_save[QUADTREE_DEPTH];
	double rate;
	IMAGE *img;
	ENCODER *enc;
	PMODEL **pmlist_save;
	double elapse = 0.0;
	int f_mmse = 0;
	int f_optpred = 0;
	int quadtree_depth = QUADTREE_DEPTH;
	int num_class = NUM_CLASS;
	int num_group = NUM_GROUP;
	int prd_order = PRD_ORDER;
	int coef_precision = COEF_PRECISION;
	int num_pmodel = NUM_PMODEL;
	int pm_accuracy = PM_ACCURACY;
	int max_iteration = MAX_ITERATION;
	char num_threads = NUM_THREADS;
	char *infile, *outfile;
	FILE *fp;

	cpu_time();
	setbuf(stdout, 0);
	infile = outfile = NULL;
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
		case 'M':
			num_class = atoi(argv[++i]);
			if (num_class <= 0 || num_class > MAX_CLASS) {
				num_class = NUM_CLASS;
			}
			break;
		case 'K':
			prd_order = atoi(argv[++i]);
			if (prd_order <= 0 || prd_order > MAX_PRD_ORDER) {
				prd_order = PRD_ORDER;
			}
			break;
		case 'P':
			coef_precision = atoi(argv[++i]);
			if (coef_precision <= 0 || coef_precision > 16) {
				coef_precision = COEF_PRECISION;
			}
			break;
		case 'V':
			num_pmodel = atoi(argv[++i]);
			if (num_pmodel <= 0 || num_pmodel > 64) {
				num_pmodel = NUM_PMODEL;
			}
			break;
		case 'A':
			pm_accuracy = atoi(argv[++i]);
			if (pm_accuracy < 0 || pm_accuracy > 6) {
				pm_accuracy = PM_ACCURACY;
			}
			break;
		case 'I':
			max_iteration = atoi(argv[++i]);
			if (max_iteration <= 0) {
				max_iteration = MAX_ITERATION;
			}
			break;
		case 'm':
			f_mmse = 1;
			break;
		case 'o':
			f_optpred = 1;
			break;
		case 'f':
			quadtree_depth = -1;
			break;
		case 'n':
			num_threads = atoi(argv[++i]);
			break;
		default:
			fprintf(stderr, "Unknown option: %s!\n", argv[i]);
			exit (1);
			}
		} else {
			if (infile == NULL) {
				infile = argv[i];
			} else {
				outfile = argv[i];
			}
		}
	}

	if (pm_accuracy > coef_precision) pm_accuracy = coef_precision;
	if (infile == NULL || outfile == NULL) {
		printf(BANNER"\n", 0.01 * VERSION);
		printf("usage: encmrp [options] infile outfile\n");
		printf("options:\n");
		printf("    -M num  Number of predictors [%d]\n", num_class);
		printf("    -K num  Prediction order [%d]\n", prd_order);
		printf("    -P num  Precision of prediction coefficients (fractional bits) [%d]\n", coef_precision);
		printf("    -V num  Number of probability models [%d]\n", num_pmodel);
		printf("    -A num  Accuracy of probability models [%d]\n", pm_accuracy);
		printf("    -I num  Maximum number of iterations [%d]\n", max_iteration);
		printf("    -m      Use MMSE predictors\n");
		printf("    -f      Fixed block-size for adaptive prediction\n");
		printf("    -o      Further optimization of predictors (experimental)\n");
		printf("infile:     Input file (must be in a raw PGM format)\n");
		printf("outfile:    Output file\n");
		exit(0);
	}
	omp_set_num_threads(num_threads);

	img = read_pgm(infile);

	fp = fileopen(outfile, "wb");

	k = img->width * img->height;
	if (num_class < 0) {
		num_class = (int)(10.4E-5 * k + 13.8);
		if (num_class > MAX_CLASS) num_class = MAX_CLASS;
	}
	if (prd_order < 0) {
		prd_order = (int)(12.0E-5 * k + 17.2);
		for (i = 1; i < 8; i++) {
			if (prd_order < (i + 1) * (i+1)) {
				prd_order = i * (i+1);
				break;
			}
		}
		if (i >= 8) prd_order = 72;
	}
#if AUTO_DEL_CL
	num_class = MAX_CLASS;
#endif

#if TEMPLATE_MATCHING_ON
	// num_class++;
	int temp_peak_num_save=0, *w_gr_save = 0;
#endif

	printf("%s -> %s (%dx%d)\n", infile, outfile, img->width, img->height);
	printf("------------------------------------------------------------------\n");
#if AUTO_PRD_ORDER
	// printf("M = %d, K = %d, P = %d, V = %d, A = %d, l = %d, m = %d, o = %d, f = %d\n",
	printf("NUM_CLASS\t= %d\nMAX_PRD_ORDER\t= %d\ncoef_precision\t= %d\nnum_pmodel\t= %d\npm_accuracy\t= %d\nmax_iteration\t= %d\nf_mmse\t\t= %d\nf_optpred\t= %d\nquadtree_depth\t= %d\nParallel Threads= %d\nTemplateM\t= %d\n",
		num_class, MAX_PRD_ORDER, coef_precision, num_pmodel, pm_accuracy, max_iteration, f_mmse, f_optpred, quadtree_depth, num_threads, TEMPLATE_MATCHING_ON);

	#if TEMPLATE_MATCHING_ON
		printf("TM_CLASS_NUM\t= %d\n", TEMPLATE_CLASS_NUM);
	#endif
#else
	printf("M = %d, K = %d, P = %d, V = %d, A = %d\n",
		num_class, prd_order, coef_precision, num_pmodel, pm_accuracy);
#endif
	printf("------------------------------------------------------------------\n");
	enc = init_encoder(img, num_class, num_group, prd_order, coef_precision,
		quadtree_depth, num_pmodel, pm_accuracy);
	enc->pmodels = init_pmodels(enc->num_group, enc->num_pmodel,
		enc->pm_accuracy, NULL, enc->sigma,
		enc->maxval + 1);

#if defined(_WIN32) && (LOG_LIST_MODE || LOG_PUT_OUT_ENC)
	if (set_directory()) {
		exit(1);
	}
#endif

#if LOG_LIST_MODE
	init_log_sheet(enc, infile);
#endif

	set_cost_model(enc, f_mmse);
	init_class(enc);
	count_cl(enc);
	num_class_save = num_class;
	prd_save = (int **)alloc_2d_array(enc->num_class, enc->max_prd_order, sizeof(int));
	th_save = (int **)alloc_2d_array(enc->num_class, enc->num_group, sizeof(int));
	class_save = (char **)alloc_2d_array(enc->height, enc->width, sizeof(char));
	mask_save = (char **)alloc_2d_array(enc->height, enc->width, sizeof(char));

	if (enc->quadtree_depth > 0) {
		y = (enc->height + MAX_BSIZE - 1) / MAX_BSIZE;
		x = (enc->width + MAX_BSIZE - 1) / MAX_BSIZE;
		for (i = enc->quadtree_depth - 1; i >= 0; i--) {
			qtmap_save[i] = (char **)alloc_2d_array(y, x, sizeof(char));
			y <<= 1;
			x <<= 1;
		}
	}
	pmlist_save = (PMODEL **)alloc_mem(enc->num_group * sizeof(PMODEL *));

#if TEMPLATE_MATCHING_ON
	tempm_array = (int ***)alloc_3d_array(enc->height, enc->width, MAX_DATA_SAVE_DOUBLE, sizeof(int));
	init_3d_array(tempm_array, enc->height, enc->width, MAX_DATA_SAVE_DOUBLE, 0);
	exam_array = (int ***)alloc_3d_array(enc->height, enc->width, TEMPLATE_CLASS_NUM, sizeof(int));
	init_3d_array(exam_array, enc->height, enc->width, TEMPLATE_CLASS_NUM, 0);
	TemplateM(enc, outfile);
	w_gr_save = (int *)alloc_mem(enc->num_group * sizeof(int));
	for(gr=0; gr<enc->num_group-1; gr++){
		enc->w_gr[gr] = 0;
	}
	enc->w_gr[enc->num_group-1] = MAX_UPARA;
#endif

	/* 1st loop */
	//単峰性確率モデルによる算術符号化
	enc->optimize_loop = 1;
	min_cost = INT_MAX;
	for (i = j = 0; i < max_iteration; i++) {
		printf("[%2d] cost =", i);
		cost = design_predictor(enc, f_mmse);	//予測器の設計(クラス情報を含んで再設計)
		printf(" %d ->", (int)cost);
		cost = optimize_group(enc);	//閾値に対する分散の最適化
		printf(" %d ->", (int)cost);
		cost = optimize_class(enc);	//クラス情報の最適化
		printf(" %d", (int)cost);
		if (cost < min_cost) {
			printf(" *\n");
			min_cost = cost;
			j = i;
			for (y = 0; y < enc->height; y++) {
				for (x = 0; x < enc->width; x++) {
					class_save[y][x] = enc->class[y][x];	//クラス情報の保存
				}
			}
			for (cl = 0; cl < enc->num_class; cl++) {
				for (k= 0; k < enc->max_prd_order; k++) {
					prd_save[cl][k] = enc->predictor[cl][k];	//予測器の保存
				}
				for (k= 0; k < enc->num_group; k++) {
					th_save[cl][k] = enc->th[cl][k];		//閾値の保存
				}
			}
		} else {
			printf("\n");
		}
		if (i - j >= EXTRA_ITERATION) break;
		elapse += cpu_time();
	}	//loop fin


	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			enc->class[y][x] = class_save[y][x];	//コストが最小となるクラスをenc構造体に保存
		}
	}
	for (cl = 0; cl < enc->num_class; cl++) {
		for (k= 0; k < enc->max_prd_order; k++) {
			enc->predictor[cl][k] = prd_save[cl][k];	//クラスに対する予測器をenc構造体に保存
		}
		i = 0;
		for (k= 0; k < enc->num_group; k++) {
			enc->th[cl][k] = th_save[cl][k];		//閾値をenc(ry
			for (; i < enc->th[cl][k]; i++) {
				enc->uquant[cl][i] = k;
			}
		}
	}
	set_cost_rate(enc);	//レート配列を更新
#if AUTO_PRD_ORDER
	set_prd_pels(enc);	//予測器の情報を更新
#endif
	save_prediction_value(enc);	//クラスごとの予測値の保存

	predict_region(enc, 0, 0, enc->height, enc->width);	//クラスに対応した予測値を保存
	cost = calc_cost(enc, 0, 0, enc->height, enc->width);	//現状のコストを算出
	printf("cost = %d\n", (int)cost);

/* init mask */
	init_mask();	//MASK構造体のメモリ確保
	for (i = 0; i < enc->height; i++){
		for (j = 0; j < enc->width; j++){
			enc->mask[i][j] = INIT_MASK;	//マスクの初期化(=0)
		}
	}
	// printf("INIT_MASK = %d\n", INIT_MASK);	//マスクサイズのパラメータ(1,9,25...) / 0であればマスクを用いた多峰性確率モデルではなく，従来の単峰性確率モデル
	set_weight_flag(enc);	//各画素毎にマスクのサイズに応じた隣のブロックにかかる画素の数を算出

	/* 2nd loop */
	//マスクによる多峰性確率モデルの作成および算術符号化
	enc->optimize_loop = 2;
	min_cost = INT_MAX;
	sw = 0;
	l=0;
#if TEMPLATE_MATCHING_ON
	enc->temp_peak_num = TEMPLATE_CLASS_NUM;
	mask->temp_cl = enc->temp_cl;
#endif
#if RENEW_ADC
	int cost_save=0, flg=0, before_cost=0;
	RESTORE_SIDE *min_cost_side, *before_side;
	min_cost_side = init_renew_adc(enc);
	before_side = init_renew_adc(enc);
#endif
	for (i = j = 0; i < max_iteration; i++) {
		printf("(%2d)cost=", i);
		if (f_optpred) {
			cost = optimize_predictor(enc);	//予測器の最適化
			printf(" %d", (int)cost);
		}
		side_cost = sc = encode_predictor(NULL, enc, 1);	//予測器の符号量を見積もる
		printf("[%d]->", (int)sc);
#if TEMPLATE_MATCHING_ON
		cost = optimize_template(enc);
		side_cost += sc = encode_w_gr_threshold(NULL, enc, 1);
		printf(" %d(%2d)[%d]->", (int)cost, enc->temp_peak_num, (int)sc);
#endif
#if OPTIMIZE_MASK_LOOP
		cost = optimize_group_mult(enc);
#else
		cost = optimize_group(enc);	//閾値に対する分散の再決定
#endif
		side_cost += sc = encode_threshold(NULL, enc, 1);	//閾値の符号量を見積もる
		printf(" %d[%d]->", (int)cost, (int)sc);
		cost = optimize_class(enc);		//クラス情報の最適化
		side_cost += sc = encode_class(NULL, enc, 1);	//クラス情報の符号量を見積もる
		printf(" %d[%d](%d)", (int)cost, (int)sc, (int)side_cost);
#if CHECK_CLASS
		for(cl=0; cl<enc->num_class; cl++){
			if(enc->num_nzcoef[cl] == -1){
				printf("[TEMP_CL: %d]" ,cl);
				break;
			}
		}
#endif
#if OPTIMIZE_MASK_LOOP
		// if (sw != 0) {
		cost = optimize_mask(enc);
		printf("-> %d", (int)cost);
		// }else{
			// set_weight_flag(enc);
		// }
#endif
		cost += side_cost;
#if AUTO_DEL_CL
		if (sw != 0) {	//コスト削減に一度でも失敗した場合に入る
			if( enc->num_class > 1 && l < 3) {
			#if PAST_ADC
				sw = enc->num_class;
				cost = auto_del_class(enc, cost);	//使用していないクラスの削除
				cl = 1;
				while (sw != enc->num_class && cl < 5) {	//削除に成功する間
					sw = enc->num_class;
					cost = auto_del_class(enc, cost);
					cl++;
				}
			#elif RENEW_ADC
				// sw = enc->num_class;
				save_info(enc, min_cost_side, 0);
				yy = xx =0;
				cost_save = min_cost;
				before_cost = cost;
				flg = 0;
				while(yy - xx < (EXTRA_ITERATION / 3) && yy < MAX_DEL_CLASS) {
					#if CHECK_DEBUG
						printf("\n");
					#endif
					sw = enc->num_class;
					cost = auto_del_class(enc, cost);
					if(cost < cost_save){
						#if CHECK_DEBUG
							printf(" *");
						#endif
						cost_save = cost;
						xx = yy;
						save_info(enc, min_cost_side, 0);
						flg = 1;
					} else if(cost < before_cost && flg != 1){
						#if CHECK_DEBUG
							printf(" +");
						#endif
						before_cost = cost;
						save_info(enc, before_side, 0);
						flg = 2;
					} else {
						if(flg ==4){
							break;
						} else if(flg == 3){
							flg++;
						} else if(flg == 0){
							flg = 3;
						}
					}
					yy++;
					#if CHECK_DEBUG
						printf("[%d]", flg);
					#endif
					if(sw == enc->num_class)	break;
					if(enc->num_class <= 1)	break;
					// break;
				}
				if(flg == 1){
					save_info(enc, min_cost_side, 1);
					l=0;
				} else if(flg == 2){
					save_info(enc, before_side, 1);
				} else {
					save_info(enc, min_cost_side, 1);
					l++;
				}
				// cost = calc_cost2(enc, 0, 0, enc->height, enc->width);
				cost = calc_side_info(enc, calc_cost(enc, 0, 0, enc->height, enc->width));
			#endif
				printf("->%d[%d]", (int)cost, enc->num_class);
			}
		}
#endif

#if OPTIMIZE_MASK_LOOP
		set_weight_flag(enc);
#endif
		printf("--> %d", (int)cost);
		if (cost < min_cost) {
			printf(" *\n");
			min_cost = cost;
			sw = 0;
			j = i;
			if (f_optpred) {
				num_class_save = enc->num_class;
				for (y = 0; y < enc->height; y++) {
					for (x = 0; x < enc->width; x++) {
						class_save[y][x] = enc->class[y][x];
						mask_save[y][x] = enc->mask[y][x];
					}
				}
				if (enc->quadtree_depth > 0) {
					y = (enc->height + MAX_BSIZE - 1) / MAX_BSIZE;
					x = (enc->width + MAX_BSIZE - 1) / MAX_BSIZE;
					for (k = enc->quadtree_depth - 1; k >= 0; k--) {
						for (yy = 0; yy < y; yy++) {
							for (xx = 0; xx < x; xx++) {
								qtmap_save[k][yy][xx] = enc->qtmap[k][yy][xx];
							}
						}
						y <<= 1;
						x <<= 1;
					}
				}
				for (gr = 0; gr < enc->num_group; gr++) {
					pmlist_save[gr] = enc->pmlist[gr];
				}
				for (cl = 0; cl < enc->num_class; cl++) {
					for (k= 0; k < enc->max_prd_order; k++) {
						prd_save[cl][k] = enc->predictor[cl][k];
					}
					for (k= 0; k < enc->num_group; k++) {
						th_save[cl][k] = enc->th[cl][k];
					}
				}
#if AUTO_DEL_CL
				l--;
#endif
#if TEMPLATE_MATCHING_ON
				temp_peak_num_save = enc->temp_peak_num;
				// w_gr_save = enc->w_gr;
				for(gr=0; gr<enc->num_group; gr++){
					w_gr_save[gr] = enc->w_gr[gr];
				}
#endif
			}
		} else {
			sw = 1;
			printf("\n");
		}
		if (f_optpred) {
			if (i - j >= EXTRA_ITERATION) break;
		} else {
			if (i > j) break;
		}
		elapse += cpu_time();
	}	//loop fin

#if RENEW_ADC
	free(min_cost_side->th_s);
	free(min_cost_side->prd_s);
	free(min_cost_side->uq_s);
	free(min_cost_side->prd_cl_s);
	free(min_cost_side->class);
	free(min_cost_side);
	free(before_side->th_s);
	free(before_side->prd_s);
	free(before_side->uq_s);
	free(before_side->prd_cl_s);
	free(before_side->class);
	free(before_side);
#endif

	if (f_optpred) {
		enc->num_class = num_class_save;
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				enc->class[y][x] = class_save[y][x];
				enc->mask[y][x] = mask_save[y][x];
			}
		}
		if (enc->quadtree_depth > 0) {
			y = (enc->height + MAX_BSIZE - 1) / MAX_BSIZE;
			x = (enc->width + MAX_BSIZE - 1) / MAX_BSIZE;
			for (k = enc->quadtree_depth - 1; k >= 0; k--) {
				for (yy = 0; yy < y; yy++) {
					for (xx = 0; xx < x; xx++) {
						enc->qtmap[k][yy][xx] = qtmap_save[k][yy][xx];
					}
				}
				y <<= 1;
				x <<= 1;
			}
		}
		for (gr = 0; gr < enc->num_group; gr++) {
			enc->pmlist[gr] = pmlist_save[gr];
		}
		for (cl = 0; cl < enc->num_class; cl++) {
			for (k= 0; k < enc->max_prd_order; k++) {
				enc->predictor[cl][k] = prd_save[cl][k];
			}
			i = 0;
			for (k= 0; k < enc->num_group; k++) {
				enc->th[cl][k] = th_save[cl][k];
				for (; i < enc->th[cl][k]; i++) {
					enc->uquant[cl][i] = k;
				}
			}
		}
#if TEMPLATE_MATCHING_ON
		enc->temp_peak_num = temp_peak_num_save;
		// enc->w_gr = w_gr_save;
		for(gr=0; gr<enc->num_group; gr++){
			enc->w_gr[gr] = w_gr_save[gr];
		}
#endif

#if AUTO_PRD_ORDER
		set_prd_pels(enc);
#endif
		save_prediction_value(enc);
		enc->function_number = 100;
		predict_region(enc, 0, 0, enc->height, enc->width);
		calc_cost(enc, 0, 0, enc->height, enc->width);
		printf("5\n");
#if OPTIMIZE_MASK
#if 0
        cost = optimize_mask(enc);
        printf(" optimize_mask: cost = %d (%d)\n", (int)cost, (int)side_cost);

		cost = optimize_group_mult(enc);
	  sc = encode_threshold(NULL, enc, 1);
		printf("optimize_group_mult: cost = %d[%d] \n", (int)cost, (int)sc);

        cost = optimize_mask(enc);
        printf(" optimize_mask: cost = %d (%d)\n", (int)cost, (int)side_cost);

#endif

#endif
	}
	remove_emptyclass(enc);

#if AUTO_PRD_ORDER
	set_prd_pels(enc);
#endif
//Start Encode
	bits = header_info = write_header(enc, fp);
	printf("header info.\t:%10d bits\n", header_info);
	enc->rc = rc_init();
	putbits(fp, 7, 0);	/* byte alignment for the rangecoder */
	bits += class_info = encode_class(fp, enc, 1);
	printf("class info.\t:%10d bits\n", class_info);
	bits += pred_info = encode_predictor(fp, enc, 1);
	printf("predictors\t:%10d bits\n", pred_info);
	bits += th_info = encode_threshold(fp, enc, 1);
	printf("thresholds\t:%10d bits\n", th_info);
#if OPTIMIZE_MASK
	bits += mask_info = encode_mask(fp, enc, 1);
	printf("mask_info\t:%10d bits\n", mask_info);
#endif
#if TEMPLATE_MATCHING_ON
	bits += w_gr_info = encode_w_gr_threshold(fp ,enc, 1);
	printf("w_gr\t\t:%10d bits\n", w_gr_info);
#endif
	bits += err_info = encode_image(fp, enc);
	printf("pred. errors\t:%10d bits\n", err_info);
	printf("------------------------------\n");
	printf("total\t\t:%10d bits\n", bits);
	rate = (double)bits / (enc->height * enc->width);
	printf("coding rate\t:%10.5f b/p\n", rate);
	fclose(fp);
	elapse += cpu_time();
	printf("cpu time: %.2f sec.\n", elapse);

	enc->optimize_loop = -1;	// Coding Done
	calc_ratio_of_model_to_rate(enc);

#if LOG_LIST_MODE
	finish_log_sheet(enc, header_info, class_info, pred_info, th_info, err_info, w_gr_info, bits, rate, elapse);
#endif

#if LOG_PUT_OUT_ENC
	print_predictor(enc->predictor, enc->max_prd_order, enc->num_class, enc->max_coef, outfile);
	#if TEMPLATE_MATCHING_ON
		print_threshold(enc->th, enc->num_group, enc->num_class, enc->pmlist, NULL, enc->w_gr, outfile);
	#else
		print_threshold(enc->th, enc->num_group, enc->num_class, enc->pmlist, NULL, NULL, outfile);
	#endif
	// print_class(enc->class, enc->num_class, enc->height, enc->width, outfile);
	print_class_color(enc->class, enc->num_class, enc->height, enc->width, outfile);
	output_class_map(enc->class, enc->num_class, enc->height, enc->width, outfile);
	print_block_size(enc->org, enc->qtmap, enc->quadtree_depth, enc->height, enc->width, outfile);
	print_class_and_block(enc->class, enc->num_class, enc->qtmap, enc->quadtree_depth, enc->height, enc->width, outfile);
	print_mask(enc->mask, enc->height, enc->width, outfile);
	print_amp_chara(enc->predictor, enc->max_prd_order, enc->num_class, enc->height, enc->width, outfile);
	print_rate_map(enc, outfile);
	calc_var_upara(enc, outfile);
	print_rate_compare_map(enc, outfile);
	print_rate_compare_class_map(enc, outfile);
#endif
	return (0);
}
