#if defined(_WIN32)
# include <windows.h>	// For "make_direcotry"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mrp.h"

extern CPOINT dyx[];
extern MASK *mask;
extern  int win_sample[];
extern int mask_x[], mask_y[];

char pallet[256][3] = {
{255, 128, 128}, {255, 255, 128}, {128, 255, 128}, {0, 255, 128},   {128, 255, 255},
{0, 128, 255},   {255, 128, 192}, {255, 128, 255}, {255, 0, 0},     {255, 255, 0},
{128, 255, 0},   {0, 255, 64},    {0, 255, 255},   {0, 128, 192},  {192, 128, 0}, //{0, 128, 192},
{255, 0, 255},   {128, 64, 64},   {255, 128, 64},  {0, 255, 0},     {0, 128, 128},
{0, 64, 128},    {128, 128, 255}, {128, 0, 64},    {255, 0, 128},   {128, 0, 0},
{255, 128, 0},   {0, 128, 0},     {0, 128, 64},    {0, 0, 255},     {0, 0, 160},
{128, 0, 128},   {128, 0, 255},   {32, 32, 32},    {33, 33, 33},    {34, 34, 34},
{35, 35, 35},    {36, 36, 36},    {37, 37, 37},    {38, 38, 38},    {39, 39, 39},
{40, 40, 40},    {41, 41, 41},    {42, 42, 42},    {43, 43, 43},    {44, 44, 44},
{45, 45, 45},    {46, 46, 46},    {47, 47, 47},    {48, 48, 48},    {49, 49, 49},
{50, 50, 50},    {51, 51, 51},    {52, 52, 52},    {53, 53, 53},    {54, 54, 54},
{55, 55, 55},    {56, 56, 56},    {57, 57, 57},    {58, 58, 58},    {59, 59, 59},
{60, 60, 60},    {61, 61, 61},    {62, 62, 62},    {63, 63, 63},    {64, 64, 64},
{65, 65, 65},    {66, 66, 66},    {67, 67, 67},    {68, 68, 68},    {69, 69, 69},
{70, 70, 70},    {71, 71, 71},    {72, 72, 72},    {73, 73, 73},    {74, 74, 74},
{75, 75, 75},    {76, 76, 76},    {77, 77, 77},    {78, 78, 78},    {79, 79, 79},
{80, 80, 80},    {81, 81, 81},    {82, 82, 82},    {83, 83, 83},    {84, 84, 84},
{85, 85, 85},    {86, 86, 86},    {87, 87, 87},    {88, 88, 88},    {89, 89, 89},
{90, 90, 90},    {91, 91, 91},    {92, 92, 92},    {93, 93, 93},    {94, 94, 94},
{95, 95, 95},    {96, 96, 96},    {97, 97, 97},    {98, 98, 98},    {99, 99, 99},
{100, 100, 100}, {101, 101, 101}, {102, 102, 102}, {103, 103, 103}, {104, 104, 104},
{105, 105, 105}, {106, 106, 106}, {107, 107, 107}, {108, 108, 108}, {109, 109, 109},
{110, 110, 110}, {111, 111, 111}, {112, 112, 112}, {113, 113, 113}, {114, 114, 114},
{115, 115, 115}, {116, 116, 116}, {117, 117, 117}, {118, 118, 118}, {119, 119, 119},
{120, 120, 120}, {121, 121, 121}, {122, 122, 122}, {123, 123, 123}, {124, 124, 124},
{125, 125, 125}, {126, 126, 126}, {127, 127, 127}, {128, 128, 128}, {129, 129, 129},
{130, 130, 130}, {131, 131, 131}, {132, 132, 132}, {133, 133, 133}, {134, 134, 134},
{135, 135, 135}, {136, 136, 136}, {137, 137, 137}, {138, 138, 138}, {139, 139, 139},
{140, 140, 140}, {141, 141, 141}, {142, 142, 142}, {143, 143, 143}, {144, 144, 144},
{145, 145, 145}, {146, 146, 146}, {147, 147, 147}, {148, 148, 148}, {149, 149, 149},
{150, 150, 150}, {151, 151, 151}, {152, 152, 152}, {153, 153, 153}, {154, 154, 154},
{155, 155, 155}, {156, 156, 156}, {157, 157, 157}, {158, 158, 158}, {159, 159, 159},
{160, 160, 160}, {161, 161, 161}, {162, 162, 162}, {163, 163, 163}, {164, 164, 164},
{165, 165, 165}, {166, 166, 166}, {167, 167, 167}, {168, 168, 168}, {169, 169, 169},
{170, 170, 170}, {171, 171, 171}, {172, 172, 172}, {173, 173, 173}, {174, 174, 174},
{175, 175, 175}, {176, 176, 176}, {177, 177, 177}, {178, 178, 178}, {179, 179, 179},
{180, 180, 180}, {181, 181, 181}, {182, 182, 182}, {183, 183, 183}, {184, 184, 184},
{185, 185, 185}, {186, 186, 186}, {187, 187, 187}, {188, 188, 188}, {189, 189, 189},
{190, 190, 190}, {191, 191, 191}, {192, 192, 192}, {193, 193, 193}, {194, 194, 194},
{195, 195, 195}, {196, 196, 196}, {197, 197, 197}, {198, 198, 198}, {199, 199, 199},
{200, 200, 200}, {201, 201, 201}, {202, 202, 202}, {203, 203, 203}, {204, 204, 204},
{205, 205, 205}, {206, 206, 206}, {207, 207, 207}, {208, 208, 208}, {209, 209, 209},
{210, 210, 210}, {211, 211, 211}, {212, 212, 212}, {213, 213, 213}, {214, 214, 214},
{215, 215, 215}, {216, 216, 216}, {217, 217, 217}, {218, 218, 218}, {219, 219, 219},
{220, 220, 220}, {221, 221, 221}, {222, 222, 222}, {223, 223, 223}, {224, 224, 224},
{225, 225, 225}, {226, 226, 226}, {227, 227, 227}, {228, 228, 228}, {229, 229, 229},
{230, 230, 230}, {231, 231, 231}, {232, 232, 232}, {233, 233, 233}, {234, 234, 234},
{235, 235, 235}, {236, 236, 236}, {237, 237, 237}, {238, 238, 238}, {239, 239, 239},
{240, 240, 240}, {241, 241, 241}, {242, 242, 242}, {243, 243, 243}, {244, 244, 244},
{245, 245, 245}, {246, 246, 246}, {247, 247, 247}, {248, 248, 248}, {249, 249, 249},
{250, 250, 250}, {251, 251, 251}, {252, 252, 252}, {253, 253, 253}, {254, 254, 254},
{0, 0, 0}
};

/* Write Predictive Coef */
void print_predictor(int **predictor, int prd_order, int num_class, int max_coef, char *outfile)
{
	int i, j, cl, k, count, count_all, **p_prd;
	int max_mhd = 10;
	char *name;
	char file[256];
	FILE *fp;

	name = strrchr( outfile, BS);
	name++;
	sprintf(file, LOG_PRED_DIR"%s_prd.csv", name);
	fp = fileopen(file, "wb");
	p_prd = (int **)alloc_2d_array(max_mhd + 1, max_mhd << 1, sizeof(int));
	count_all = 0;

	for (cl = 0; cl < num_class; cl++) {
		count = 0;
		fprintf(fp, "\nClass,%d,\n", cl);
		for (i = 0; i < max_mhd + 1; i++) {
			for (j = 0; j < max_mhd << 1; j++) {
				p_prd[i][j] = max_coef + 1;
			}
		}
		for (k = 0; k < prd_order; k++) {
			p_prd[max_mhd + dyx[k].y][max_mhd + dyx[k].x] = predictor[cl][k];
		}
		for (i = 0; i < max_mhd + 1; i++) {
			for (j = 0; j < max_mhd << 1; j++) {
				if (p_prd[i][j] == max_coef + 1) fprintf(fp, ",");
				else{
					 fprintf(fp, "%d,", p_prd[i][j]);
					 if (p_prd[i][j] != 0) count++;
				}
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\nNz_Order,%d,\n", count);
		count_all += count;
	}
	fprintf(fp, "\nAVE,%f,\n", (double)count_all/num_class);
	fclose(fp);
	free(p_prd);
	return;
}

/* Write Threshold */
void print_threshold(int **th, int num_group, int num_class,
					 PMODEL **pmlist, int *pm_idx, char *outfile)
{
	int cl, gr;
	char *name;
	char file[256];
	FILE *fp;

	name = strrchr( outfile, BS);
	name++;
	sprintf(file, LOG_TH_DIR"%s_th.csv", name);
	fp = fileopen(file, "wb");
	for (cl = 0; cl < num_class; cl++){
		fprintf(fp, "[%d],", cl);
	}
	fprintf(fp, "\n");
	for (cl = 0; cl < num_class; cl++) {
		fprintf(fp, "[%d],", cl);
		for (gr = 0; gr < num_group - 1; gr++) {
			fprintf(fp, "%d,", th[cl][gr]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	if (pmlist != NULL) {
		for (gr = 0; gr < num_group; gr++) {
			fprintf(fp, "%d,", pmlist[gr]->id);
		}
		fprintf(fp, "\n");
	}
	if (pm_idx != NULL) {
		for (gr = 0; gr < num_group; gr++) {
			fprintf(fp, "%d,", pm_idx[gr]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	return;
}

/* Draw Class Map */
void print_class(char **class, int num_class, int height, int width, char *outfile)
{
	int i, j, step;
	char *name;
	char file[256];
	FILE *fp;

	name = strrchr( outfile, BS);
	name++;
	sprintf(file, LOG_CL_DIR"%s_class.pgm", name);
	fp = fileopen(file, "wb");
	step = 255 / num_class;
	fprintf(fp, "P5\n%d %d\n255\n", width, height);
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			putc(class[i][j] * step, fp);
		}
	}
	fclose(fp);
	return;
}

void output_class_map(char **class, int num_class, int height, int width, char *outfile)
{
	int i, j;
	int count_class[num_class];
	char *name;
	char file[256];
	FILE *fp;

	memset(&count_class, 0, sizeof(count_class));

	name = strrchr( outfile, BS);
	name++;
	sprintf(file, LOG_CL_DIR"%s_class_map.csv", name);
	fp = fileopen(file, "wb");
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			fprintf(fp , "%d,", class[i][j]);
			count_class[(int)class[i][j]]++;
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n\n\nClass Histogram\n");
	for(i=0; i<num_class; i++){
		fprintf(fp, "%d,%d,%f%%\n", i, count_class[i], (double)count_class[i] * 100 / (height * width));
	}
	fclose(fp);
	return;
}

void print_class_color(char **class, int num_class, int height, int width, char *outfile)
{
	int y, x;
	char *name;
	char file[256];
	FILE *fp;

	name = strrchr( outfile, BS);
	if (name == NULL) {
		name = outfile;
	}else {
		name++;
	}
	sprintf(file, LOG_CL_DIR"%s_class.ppm", name);
	fp = fileopen(file, "wb");
	fprintf(fp, "P6\n%d %d\n255\n", width, height);
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			putc(pallet[(int)(class[y][x])][0],fp);
			putc(pallet[(int)(class[y][x])][1],fp);
			putc(pallet[(int)(class[y][x])][2],fp);
		}
	}
	fclose(fp);
	return;
}

/* Draw Block Size */
void print_block_size(int **org, char ***qtmap, int quadtree_depth,
					  int height, int width, char *outfile)
{
	int y, x, i, j, depth, bs;
	int **level;
	char *name;
	char file[256];
	FILE *fp;

	name = strrchr(outfile, BS);
	if (name == NULL) {
		name = outfile;
	}else {
		name++;
	}
	sprintf(file, LOG_VBS_DIR"%s_blocksize.pgm", name);
	fp = fileopen(file, "wb");
	fprintf(fp, "P5\n%d %d\n255\n", width, height);
	level = (int **)alloc_2d_array(height, width, sizeof(int));
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			level[y][x] = 0;
		}
	}
	if (quadtree_depth > 0) {
		for (depth = 0; depth < quadtree_depth; depth++) {
			for (y = 0; y < height; y++) {
				for (x = 0; x < width; x++) {
					j = y / (MAX_BSIZE >> depth);
					i = x / (MAX_BSIZE >> depth);
					if (qtmap[quadtree_depth - depth - 1][j][i] == 1) {
						level[y][x]++;
					}
				}
			}
		}
	}else {
		depth = MAX_BSIZE / BASE_BSIZE / 2;
		for (y = 0; y < height; y++) {
			for (x = 0; x < width; x++) {
				level[y][x] = depth;
			}
		}
	}
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			depth = level[y][x];
			bs = MAX_BSIZE >> depth;
			if (y % bs == 0 || x %  bs == 0 || y == height - 1 || x == width - 1) {
				putc(255, fp);
			}else {
				putc(org[y][x], fp);
			}
		}
	}
	fclose(fp);
	free(level);
	return;
}

void print_class_and_block(char **class, int num_class, char ***qtmap, int quadtree_depth,
						   int height, int width, char *outfile)
{
	int y, x, i, j, depth, bs;
	int **level;
	char *name;
	char file[256];
	FILE *fp;

	name = strrchr( outfile, BS);
	if (name == NULL) {
		name = outfile;
	}else {
		name++;
	}
	sprintf(file, LOG_CL_DIR"%s_class_block.ppm", name);
	fp = fileopen(file, "wb");
	fprintf(fp, "P6\n%d %d\n255\n", width, height);
	level = (int **)alloc_2d_array(height, width, sizeof(int));
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			level[y][x] = 0;
		}
	}
	if (quadtree_depth > 0) {
		for (depth = 0; depth < quadtree_depth; depth++) {
			for (y = 0; y < height; y++) {
				for (x = 0; x < width; x++) {
					j = y / (MAX_BSIZE >> depth);
					i = x / (MAX_BSIZE >> depth);
					if (qtmap[quadtree_depth - depth - 1][j][i] == 1) {
						level[y][x]++;
					}
				}
			}
		}
	}else {
		depth = MAX_BSIZE / BASE_BSIZE / 2;
		for (y = 0; y < height; y++) {
			for (x = 0; x < width; x++) {
				level[y][x] = depth;
			}
		}
	}
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			depth = level[y][x];
			bs = MAX_BSIZE >> depth;
			if (y % bs == 0 || x %  bs == 0 || y == height - 1 || x == width - 1) {
				putc(255, fp);
				putc(255, fp);
				putc(255, fp);
			}else {
				putc(pallet[(int)(class[y][x])][0],fp);
				putc(pallet[(int)(class[y][x])][1],fp);
				putc(pallet[(int)(class[y][x])][2],fp);
			}
		}
	}
	fclose(fp);
	free(level);
	return;
}

void print_mask(char **mask, int height, int width, char *outfile)
{
	int i, j, step;
	char *name;
	char file[256];
	FILE *fp;

	name = strrchr( outfile, BS);
	name++;
	sprintf(file, LOG_CL_DIR"%s_mask.pgm", name);
	fp = fileopen(file, "wb");
	step = 255 / NUM_MASK;
	fprintf(fp, "P5\n%d %d\n255\n", width, height);
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			putc(mask[i][j] * step, fp);
		}
	}
	fclose(fp);
	return;
}



/* Write Predictor Characteristic (Amplitude) */
void print_amp_chara(int **predictor, int prd_order, int num_class,
					 int height, int width, char *outfile)
{
    int fx, fy, j, cl, delta_y, delta_x;
    int *coef;
    double Re, Im, H;
    char *name;
    char file[256];
    FILE *fp;

    name = strrchr( outfile, BS);
	if (name == NULL) {
		name = outfile;
	}else {
		name++;
	}
    sprintf(file, LOG_AMP_CH_DIR"%s_amp_chara.csv", name);
    fp = fileopen(file, "wb");

    for(cl = 0; cl < num_class; cl++) {
        coef = &predictor[cl][0];
        fprintf(fp, "Predictor,%d\n\n", cl);
        fprintf(fp, " , , ");
        for(fx = 0; fx < (width / 2); fx += 5) {
            fprintf(fp, "%d,", fx);
        }
        fprintf(fp, "\n");
        for(fy = 0; fy < (height / 2); fy += 5) {
            fprintf(fp, " ,%d,", fy);
            for(fx = 0; fx < (width / 2); fx += 5) {
                Re = Im = 0.0;
                for(j = 0; j < prd_order; j++) {
                    delta_y = dyx[j].y;
                    delta_x = dyx[j].x;
                    Re += coef[j]
                      * cos(2 * M_PI * (((double)fx / width) * delta_x
                                        + ((double)fy / height) * delta_y));
                    Im += coef[j]
                      * sin(2 * M_PI * (((double)fx / width) *delta_x
                                        + ((double)fy / height) * delta_y));
                }
                H = sqrt(pow(Re, 2) + pow(Im, 2));
                fprintf(fp, "%f,", H);
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n\n");
    }
    fclose(fp);
    return;
}

void init_mask_parameter(ENCODER *enc,int y, int x, int u)
{
    int ty, tx, cl, i, peak, sample;
		int m_gr,m_prd,m_frac;
    int count_cl[enc->num_class];

    sample = win_sample[(int)enc->mask[y][x]]; //window nai no gaso suu

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
				count_cl[cl]++; //mask nai no class count
    }

   for(cl = peak = 0; cl < enc->num_class; cl++){
   	if (count_cl[cl]!=0){
			mask->class[peak] =cl;
	  	mask->weight[peak] =( (count_cl[cl] << W_SHIFT )/ sample);
	  	m_gr = enc->uquant[cl][u];
	  	m_prd = enc->prd_class[y][x][cl];
	  	m_prd = CLIP(0, enc->maxprd, m_prd);
	  	mask->base[peak] = enc->bconv[m_prd];
	  	m_frac = enc->fconv[m_prd];
	  	mask->pm[peak] = enc->pmlist[m_gr] + m_frac;
	  	peak++;
	  }
   }
   mask->num_peak = peak; //peak no kazu
}

/* Draw Predictive error rate */
#if MULT_PEAK_MODE
void print_rate_map(ENCODER *enc, char *outfile)
{
	int y, x, step;
	int org, u;
	int max_rate = 15;
	double cost;
	char *name;
	char file[256];
	PMODEL *pm;
	FILE *fp;

	name = strrchr( outfile, BS);
	if (name == NULL) {
		name = outfile;
	}else {
		name++;
	}
	sprintf(file, LOG_RATE_DIR"%s_rate.pgm", name);
	fp = fileopen(file, "wb");
	step = 255 / max_rate;
	fprintf(fp, "P5\n%d %d\n255\n", enc->width, enc->height);
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			org = enc->org[y][x];
			u = enc->upara[y][x];
		  init_mask_parameter(enc,y,x,u);
		  pm = &enc->mult_pm;
		  set_pmodel_mult(pm,mask,enc->maxval+1);
			cost = (log(pm->cumfreq[enc->maxval+1])-log(pm->freq[org]))/log(2.0);
			cost *= step;
			if(cost < 0) cost = 0;
			else if(cost > 255) cost = 255;
			putc((int)cost, fp);
		}
	}
	fclose(fp);
	return;
}
#else
void print_rate_map(ENCODER *enc, char *outfile)
{
	int y, x, step;
	int gr, prd, base, frac, org;
	int max_rate = 15;
	double cost;
	char *name;
	char file[256];
	PMODEL *pm;
	FILE *fp;

	name = strrchr( outfile, BS);
	if (name == NULL) {
		name = outfile;
	}else {
		name++;
	}
	sprintf(file, LOG_RATE_DIR"%s_rate.pgm", name);
	fp = fileopen(file, "wb");
	step = 255 / max_rate;
	fprintf(fp, "P5\n%d %d\n255\n", enc->width, enc->height);
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			org = enc->org[y][x];
			gr = enc->group[y][x];
			prd = enc->prd[y][x];
			if(prd < 0) prd = 0;
			else if(prd > enc->maxprd) prd = enc->maxprd;
			base = enc->bconv[prd];
			frac = enc->fconv[prd];
			pm = enc->pmlist[gr] + frac;
			cost = pm->cost[base + org] + pm->subcost[base];
			cost *= step;
			if(cost < 0) cost = 0;
			else if(cost > 255) cost = 255;
			putc((int)cost, fp);
		}
	}
	fclose(fp);
	return;
}
#endif

void print_rate_compare_map(ENCODER *enc, char *outfile)
{
	int y, x;
	int gr, prd, base, frac, org,u;
	double cost,cost_mult;
	char *name;
	char file[256];
	PMODEL *pm;
	FILE *fp;

	name = strrchr( outfile, BS);
	if (name == NULL) {
		name = outfile;
	}else {
		name++;
	}
	sprintf(file, LOG_RATE_DIR"%s_rate_compare.pgm", name);
	fp = fileopen(file, "wb");
	fprintf(fp, "P5\n%d %d\n255\n", enc->width, enc->height);
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			org = enc->org[y][x];
			gr = enc->group[y][x];
			prd = enc->prd[y][x];
			if(prd < 0) prd = 0;
			else if(prd > enc->maxprd) prd = enc->maxprd;
			base = enc->bconv[prd];
			frac = enc->fconv[prd];
			pm = enc->pmlist[gr] + frac;
			cost = pm->cost[base + org] + pm->subcost[base];

			u = enc->upara[y][x];
		  init_mask_parameter(enc,y,x,u);
		  pm = &enc->mult_pm;
		  set_pmodel_mult(pm,mask,enc->maxval+1);
			cost_mult = (log(pm->cumfreq[enc->maxval+1])-log(pm->freq[org]))/log(2.0);
			cost = cost - cost_mult;
			cost = 1000*cost;
			cost = CLIP(0,255,128+cost);
//			printf("(y %d,x %d)deltacost = %d\n",y,x,(int)cost);
//			if((int)cost <= 0) cost = 0; //tahousei lose
//			else if((int)cost > 0) cost = 255; //tahousei Win!
			putc((int)cost, fp);
		}
	}
	fclose(fp);
	return;
}

void print_rate_compare_class_map(ENCODER *enc, char *outfile)
{
	int y, x, tly, tlx, bry, brx, i, j, depth, blksize;
	int gr, prd, base, frac, org,u,putcost;
	int **level;
	double cost,cost_mult,sum_cost;
	double **block_cost;
	char *name;
	char file[256];
	PMODEL *pm;
	FILE *fp;

	name = strrchr( outfile, BS);
	if (name == NULL) {
		name = outfile;
	}else {
		name++;
	}

	sprintf(file, LOG_RATE_DIR"%s_rate_compare_minblk.pgm", name);
	fp = fileopen(file, "wb");
	fprintf(fp, "P5\n%d %d\n255\n", enc->width, enc->height);
	block_cost = (double **)alloc_2d_array(enc->height, enc->width, sizeof(double));
	level = (int **)alloc_2d_array(enc->height, enc->width, sizeof(int));

	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			level[y][x] = 0;
		}
	}
	if (enc->quadtree_depth > 0) {
		for (depth = 0; depth < enc->quadtree_depth; depth++) {
			for (y = 0; y < enc->height; y++) {
				for (x = 0; x < enc->width; x++) {
					j = y / (MAX_BSIZE >> depth);
					i = x / (MAX_BSIZE >> depth);
					if (enc->qtmap[enc->quadtree_depth - depth - 1][j][i] == 1) {
						level[y][x]++;
					}
				}
			}
		}
	}else {
		depth = MAX_BSIZE / BASE_BSIZE / 2;
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				level[y][x] = depth;
			}
		}
	}

	for (tly = 0; tly < enc->height; tly += MIN_BSIZE) {
		bry = (tly + MIN_BSIZE < enc->height) ? (tly + MIN_BSIZE) : enc->height;

		for (tlx = 0; tlx < enc->width; tlx += MIN_BSIZE) {
			brx = (tlx + MIN_BSIZE < enc->width) ? (tlx + MIN_BSIZE) : enc->width;
			sum_cost = 0;
			for (y = tly; y < bry; y++){
				for (x = tlx; x < brx; x++){
					org = enc->org[y][x];
					gr = enc->group[y][x];
					prd = enc->prd[y][x];
					if(prd < 0) prd = 0;
					else if(prd > enc->maxprd) prd = enc->maxprd;
					base = enc->bconv[prd];
					frac = enc->fconv[prd];
					pm = enc->pmlist[gr] + frac;
					cost = pm->cost[base + org] + pm->subcost[base];

					u = enc->upara[y][x];
					init_mask_parameter(enc,y,x,u);
					pm = &enc->mult_pm;
					set_pmodel_mult(pm,mask,enc->maxval+1);
					cost_mult = (log(pm->cumfreq[enc->maxval+1])-log(pm->freq[org]))/log(2.0);
					sum_cost += cost - cost_mult;
				}
			}
			for (y = tly; y < bry; y++){
				for (x = tlx; x < brx; x++){
					block_cost[y][x] = sum_cost;
					//printf("(y %d,x %d)block_cost = %lf\n",y,x,block_cost[y][x]);
				}
			}
		}//tlx
	}//tly

	for (y = 0; y <  enc->height; y++){
		for (x = 0; x <  enc->width; x++){
			putcost = block_cost[y][x];
			if(putcost < 5) putcost = 0;
			else if(putcost >= 5) putcost = 255; //tahousei 5bit Win!
			putc((int)putcost, fp);
		}
	}
	fclose(fp);

#if 1
	sprintf(file, LOG_RATE_DIR"%s_rate_compare_minblk2.pgm", name);
	fp = fileopen(file, "wb");
	fprintf(fp, "P5\n%d %d\n255\n", enc->width, enc->height);

	for (y = 0; y <  enc->height; y++){
		for (x = 0; x <  enc->width; x++){
			putcost = block_cost[y][x];
			putcost = 10*putcost;
			putcost = CLIP(0,255,128+putcost);
			putc((int)putcost, fp);
		}
	}
	fclose(fp);
#endif


	depth = QUADTREE_DEPTH;

	for (blksize = MIN_BSIZE<<1; blksize <= MAX_BSIZE; blksize <<=1 ){
		for (tly = 0; tly < enc->height; tly += blksize) {
			bry = (tly + blksize < enc->height) ? (tly + blksize) : enc->height;
			for (tlx = 0; tlx < enc->width; tlx += blksize) {
				brx = (tlx + blksize < enc->width) ? (tlx + blksize) : enc->width;
				if (level[tly][tlx] < depth){
					sum_cost = block_cost[tly][tlx] + block_cost[tly + (blksize>>1)][tlx] +
						block_cost[tly][tlx +(blksize>>1)] + block_cost[tly + (blksize>>1)][tlx + (blksize>>1)];
					for (y = tly; y < bry; y++){
						for (x = tlx; x < brx; x++){
							block_cost[y][x] = sum_cost;
							//printf("(y %d,x %d)block_cost = %lf\n",y,x,block_cost[y][x]);
						}
					}
				}
			}
		}
		depth--;
	}

	sprintf(file, LOG_RATE_DIR"%s_rate_compare_class.pgm", name);
	fp = fileopen(file, "wb");
	fprintf(fp, "P5\n%d %d\n255\n", enc->width, enc->height);



	for (y = 0; y <  enc->height; y++){
		for (x = 0; x <  enc->width; x++){
			putcost = block_cost[y][x];
			if(putcost < 5) putcost = 0;
			else if(putcost >= 5) putcost = 255; //tahousei 5bit Win!
			putc((int)putcost, fp);
		}
	}
	fclose(fp);

#if 1
	sprintf(file, LOG_RATE_DIR"%s_rate_compare_class2.pgm", name);
	fp = fileopen(file, "wb");
	fprintf(fp, "P5\n%d %d\n255\n", enc->width, enc->height);

	for (y = 0; y <  enc->height; y++){
		for (x = 0; x <  enc->width; x++){
			putcost = block_cost[y][x];
			putcost = 10*putcost;
			putcost = CLIP(0,255,128+putcost);
			putc((int)putcost, fp);
		}
	}
	fclose(fp);
#endif

	free(block_cost);
	free (level);
	return;
}


/* Write upara - variance of predictive error */
void calc_var_upara( ENCODER *enc, char *outfile)
{
    int y, x, u, prd, base, pels, upara;
    int maximum_upara = 0;
    double ave, var;
	char *name;
	char file[256];
    FILE *fp;

	name = strrchr( outfile, BS);
	if (name == NULL) {
		name = outfile;
	}else {
		name++;
	}
    sprintf( file, LOG_VARI_DIR"%s_vari_upara.csv", name);
    fp = fileopen( file, "wb");

    for( y = 0; y < enc->height; y++) {
        for( x = 0; x < enc->width; x++) {
            if( enc->upara[y][x] > maximum_upara)
              maximum_upara = enc->upara[y][x];
        }
    }

	fprintf(fp, "Upara,Num_Pels,Variance,Deviation\n");
    for( u = 0; u <= maximum_upara; u++) {
        ave = 0;
        pels = 0;
        var = 0.0;
        for( y = 0; y < enc->height; y++) {
            for( x = 0; x < enc->width; x++) {
                upara = enc->upara[y][x];

                if( u != upara)
                  continue;

                prd = enc->prd[y][x];
                if( prd < 0) prd = 0;
                else if( prd > enc->maxprd) prd = enc->maxprd;
                base = enc->bconv[prd];

                ave += base + enc->org[y][x];
                pels++;
            }
        }
        if( pels <= 0) {
            ave = 0;
        }else {
            ave /= pels;
        }
        for( y = 0; y < enc->height; y++) {
            for( x = 0; x < enc->width; x++) {
                upara = enc->upara[y][x];

                if( u != enc->upara[y][x])
                  continue;

                prd = enc->prd[y][x];
                if( prd < 0) prd = 0;
                else if( prd > enc->maxprd) prd = enc->maxprd;
                base = enc->bconv[prd];

                var += pow( base + enc->org[y][x] - ave, 2.0);
            }
        }
        if( pels <= 0) {
            var = 0;
            fprintf( fp, "%d,%d, , \n", u, pels);
        }else {
            var /= pels;
            fprintf( fp, "%d,%d,%.4f,%.4f\n", u, pels, var, sqrt(var));
        }
    }
    fprintf( fp,"\n");
    return;
}

void init_log_sheet(ENCODER *enc, char *outfile)
{
	char *name;
	char date[64];
	time_t ct = time(NULL);
	struct tm *lst = localtime(&ct);
	FILE *fp;

	if ((name = strrchr(outfile, BS)) == NULL) {
		name = outfile;
	}else {
		name++;
	}

	sprintf(date, "%4d/%2d/%2d_%2d:%2d",
		1900 + lst->tm_year, 1 + lst->tm_mon, lst->tm_mday, lst->tm_hour, lst->tm_min);

	if ((fp = fopen(LOG_LIST, "rb")) == NULL) {
		fp = fileopen(LOG_LIST, "wb");
		fprintf(fp, "Date,Input,Height,Width,Init_Class,Use_Class,Prd_Order,\
			Header[bits],Class[bits],Predictor[bits],Threshold[bits],\
			Pred. Errors[bits],Total info.[bits],Coding Rate[bits/pel],Time[s], ,");

		fprintf(fp, "Auto_Del_CL,	Auto_Set_Coef,BASE_BSIZE,QUADTREE_DEPTH,\
			MIN_BSIZE,MAX_BSIZE,COEF_PRECISION,PM_ACCURACY,NUM_GROUP,\
			UPEL_DIST,MAX_ITERATION,CTX_WEIGHT,TEMPLATE_MATCHING,COST FUNCTION,MANHATTAN_SORT,\
			Search Window,");
		// fprintf(fp, "\n");
		fclose(fp);
	}

	fp = fileopen(LOG_LIST, "ab");
	fprintf(fp, "\n%s,%s,%d,%d,%d,", date, name, enc->height, enc->width, enc->num_class);
	fclose(fp);

	return;
}

void finish_log_sheet(ENCODER *enc, int header_info, int class_info, int pred_info,
		int th_info, int err_info, int total_info, double rate, double time)
{
	FILE *fp;

	fp = fileopen(LOG_LIST, "ab");
	fprintf(fp, "%d,%d,", enc->num_class, enc->max_prd_order);
	fprintf(fp, "%d,%d,%d,%d,%d,%d,", header_info, class_info, pred_info, th_info, err_info, total_info);
	fprintf(fp, "%f,%f, ,", rate, time);

	if(AUTO_DEL_CL) {
		fprintf(fp, "ON,");
	}else {
		fprintf(fp, "OFF,");
	}

	if(AUTO_PRD_ORDER) {
		fprintf(fp, "ON,");
	}else {
		fprintf(fp, "OFF,");
	}

	fprintf(fp, "%d,", BASE_BSIZE);
	fprintf(fp, "%d,%d,%d,", enc->quadtree_depth, MIN_BSIZE, MAX_BSIZE);
	fprintf(fp, "%d,%d,%d,%d,%d,", enc->coef_precision, enc->pm_accuracy,enc->num_group, UPEL_DIST, MAX_ITERATION);

	if(CTX_WEIGHT) {
		if(MHD_WEIGHT) {
			fprintf(fp, "mhd,");
		}else {
			fprintf(fp, "euclid,");
		}
	}else {
		fprintf(fp, "Non,");
	}

	if(TEMPLATE_MATCHING_ON){
		fprintf(fp, "ON,");
#if TEMPLATE_MATCHING_ON
		if(ZNCC){
			fprintf(fp, "ZNCC,");
		} else if(ZSAD){
			fprintf(fp, "ZSAD,");
		} else {
			fprintf(fp, "SSD,");
		}

		if(MANHATTAN_SORT){
			fprintf(fp, "ON,");
		} else {
			fprintf(fp, "OFF,");
		}

		fprintf(fp, "%d*%d,", X_SIZE*2+1, Y_SIZE);
#endif
	} else {
		fprintf(fp, "OFF,");
	}

	// fprintf(fp, "\n");
	fclose(fp);

	return;
}

#if defined(_WIN32)
int check_pathname(LPTSTR lpPathName)
{
	if (lpPathName[1] != ':' || lpPathName[2] != '\\') return 0;
	if (strpbrk(lpPathName + 2, TEXT("/,:;*\"<>|"))) return 0;
	if (strstr(lpPathName, TEXT("\\\\"))) return 0;
	if (strstr(lpPathName, TEXT(" \\"))) return 0;

	return 1;
}

int make_directory_sub(LPTSTR lpPathName)
{
	char lpPrePath[MAX_PATH], lpCorrPath[MAX_PATH];
	char *tmp;
	int flag;

	strcpy(lpCorrPath, lpPathName);
	tmp = strrchr(lpCorrPath, '\\');
	if (!tmp) return 0;
	tmp[0] = '\0';

	flag = CreateDirectory(lpCorrPath, NULL);
	if (!flag) {
		strcpy(lpPrePath, lpCorrPath);
		make_directory_sub(lpPrePath);
		flag = CreateDirectory(lpCorrPath, NULL);
	}
	return (flag);
}

int make_directory(LPTSTR lpPathName)
{
	int flag = 0;

	if (check_pathname(lpPathName)) {
		CreateDirectory(lpPathName, NULL);
		if (ERROR_ALREADY_EXISTS == GetLastError()){
			flag = 1;
		}else {
			flag = make_directory_sub(lpPathName);
		}
	}
	return (flag);
}

int set_directory(void)
{
	int flag = 1;

	flag *= make_directory(LOG_RATE_DIR);
	flag *= make_directory(LOG_VARI_DIR);
	flag *= make_directory(LOG_CL_DIR);
	flag *= make_directory(LOG_VBS_DIR);
	flag *= make_directory(LOG_TH_DIR);
	flag *= make_directory(LOG_PRED_DIR);
	flag *= make_directory(LOG_AMP_CH_DIR);

	if (flag == 0) {
		fprintf(stderr, "Can't make directory!\n");
		return (1);
	}
	return (0);
}
#endif

#if TEMPLATE_MATCHING_ON
void TemplateM_Log_Output(ENCODER *enc, char *outfile, int ***tempm_array, int ***exam_array){
	int y, x, k;
	FILE *fp;
	char *name, file[256];

	name = strrchr(outfile, BS);
	name++;
	sprintf(file, LOG_TEMP_DIR"%s_temp.csv", name);
	if ((fp = fopen(file, "rb")) != NULL) {
		if(remove(file) != 0){
			printf("Remove Error!!![%s]\n", file);
		}
	}
	fp = fileopen(file, "wb");
	// fprintf(fp, "Y_SIZE,X_SIZE,AREA,MAX_DATA_SAVE,MAX_DATA_SAVE_DOUBLE,MAX_MULTIMODAL\n");
	fprintf(fp, "%d,%d,%d,%d,%d\n",
		Y_SIZE, X_SIZE, AREA, MAX_DATA_SAVE, MAX_DATA_SAVE_DOUBLE);
	// fprintf(fp, "NAS_ACCURACY,TEMPLATE_CLASS_NUM\n");
	fprintf(fp, "%d,%d,%d,%d\n", NAS_ACCURACY, TEMPLATE_CLASS_NUM, MANHATTAN_SORT, ZNCC);

	fprintf(fp, "\n");
	for(y=0; y<enc->height; y++){
		for(x=0; x<enc->width; x++){
			for(k=0; k<MAX_DATA_SAVE_DOUBLE; k++){
				fprintf(fp, "%d,", tempm_array[y][x][k]);
			}
			for(k=0; k<MAX_DATA_SAVE; k++){
				fprintf(fp, "%d,", enc->array[y][x][k]);
			}
			for(k=0; k<TEMPLATE_CLASS_NUM; k++){
				fprintf(fp, "%d,", exam_array[y][x][k]);
			}
			fprintf(fp,"\n");
		}
	}
	fclose(fp);
	return;
}

void TemplateM_Log_Input(ENCODER *enc, char *outfile, int ***tempm_array, int ***exam_array){
	int y, x, k, y_size_check, x_size_check, area_check, max_data_save_check,
		max_data_save_double_check, nas_accuracy_check, manhattan_sort_check, zncc_check,
		template_class_num_check, flg=0;
	FILE *fp;
	char *name, file[256];

	name = strrchr(outfile, BS);
	name++;
	sprintf(file, LOG_TEMP_DIR"%s_temp.csv", name);
	if (( fp = fileopen(file, "rb")) == NULL){
		printf("[%s] is NOT exist\n", file);
		exit(1);
	}
	fscanf(fp, "%d,%d,%d,%d,%d\n", &y_size_check, &x_size_check, &area_check, &max_data_save_check,
		&max_data_save_double_check);
	if(y_size_check != Y_SIZE){
		printf("Y_SIZE is NOT coincide!!!\n");
		flg++;
	}
	if(x_size_check != X_SIZE){
		printf("X_SIZE is NOT coincide!!!\n");
		flg++;
	}
	if(area_check != AREA){
		printf("AREA is NOT coincide!!!\n");
		flg++;
	}
	if(max_data_save_check != MAX_DATA_SAVE){
		printf("MAX_DATA_SAVE is NOT coincide!!!\n");
		flg++;
	}
	if(max_data_save_double_check != MAX_DATA_SAVE_DOUBLE){
		printf("MAX_DATA_SAVE_DOUBLE is NOT coincide!!!\n");
		flg++;
	}

	fscanf(fp, "%d,%d,%d,%d\n", &nas_accuracy_check, &template_class_num_check,
		&manhattan_sort_check, &zncc_check);
	if(nas_accuracy_check != NAS_ACCURACY){
		printf("NAS_ACCURACY is NOT coincide!!!\n");
		flg++;
	}
	if(template_class_num_check < TEMPLATE_CLASS_NUM){
		printf("TEMPLATE_CLASS_NUM is NOT coincide!!!\n");
		flg++;
	}
	if(manhattan_sort_check != MANHATTAN_SORT){
		printf("MANHATTAN_SORT is NOT coincide!!!\n");
		flg++;
	}
	if(zncc_check != ZNCC){
		printf("ZNCC is NOT coincide!!!\n");
		flg++;
	}
	if(flg!=0)exit(1);

	fscanf(fp, "\n");

	for(y=0; y<enc->height; y++){
		for(x=0; x<enc->width; x++){
			for(k=0; k<MAX_DATA_SAVE_DOUBLE; k++){
				fscanf(fp, "%d,", &tempm_array[y][x][k]);
			}

			for(k=0; k<MAX_DATA_SAVE; k++){
				fscanf(fp, "%d,", &enc->array[y][x][k]);
			}

			for(k=0; k<TEMPLATE_CLASS_NUM; k++){
				fscanf(fp, "%d,", &exam_array[y][x][k]);
			}
			fscanf(fp,"\n");
		}
	}
	fclose(fp);
	return;
}

#endif
