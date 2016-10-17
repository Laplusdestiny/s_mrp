/****** LOG ********************************/
#define LOG_LIST_MODE		1
#define LOG_PUT_OUT_ENC		1
#define LOG_PUT_OUT_DEC		1

#if defined(_WIN32)		// Windows OS ?
# define BOUNDARY		"\\"
# define BS				'\\'
# define DIR			"C:\\test\\"	// Directory for output (Set this value according to your environment.)
				// [ATTENTION] This program does not support relative path!
#else
# define BOUNDARY		"/"
# define BS				'/'
# define DIR			"./Info/" // Directory for output (Must not change!)
#endif

// For Encoder only
#define LOG_LIST		DIR"Log_List.csv"
#define LOG_RATE_DIR		DIR"Rate_Map"BOUNDARY
#define LOG_VARI_DIR		DIR"Var_Upara"BOUNDARY

// For Encoder and Decoder
#define LOG_CL_DIR		DIR"Class_Map"BOUNDARY
#define LOG_VBS_DIR		DIR"Block_Size"BOUNDARY
#define LOG_TH_DIR		DIR"Threshold"BOUNDARY
#define LOG_PRED_DIR		DIR"Predictor"BOUNDARY
#define LOG_AMP_CH_DIR	DIR"Amp_Chara"BOUNDARY
#define LOG_TEMP_DIR		DIR"Temp_Map"BOUNDARY

/****** MRP-VERSION ************************/
#define MAGIC_NUMBER	('M' << 8) + 'R'
#define BANNER		"ENCMRP/DECMRP version %.2f ( Jan. 2016 )"
#define VERSION		611

/****** OPTIMIZE ***************************/
#define OPT_SIDEINFO		1 // 1 : side-info into consideration (standard), 0 : neglect side-info
#define MAX_ITERATION 	100	//100
#define EXTRA_ITERATION	10
#define AUTO_DEL_CL		1
#define AUTO_PRD_ORDER	1
#if AUTO_DEL_CL
	#define RENEW_ADC	1
	#if RENEW_ADC
		#define	PAST_ADC	0
		#define MAX_DEL_CLASS	5
	#else
		#define PAST_ADC	1
	#endif
#endif
/****** MULT PEAK **************************/
#define MULT_PEAK_MODE		1
#define OPTIMIZE_MASK 		1
#define OPTIMIZE_MASK_LOOP		1
#define WIN_BSIZE			32
#define NUM_MASK			5
#define W_SHIFT			7
#define INIT_MASK			0   // 0 1 2 3 4
#define MAX_PEAK_NUM		25

/***** BLOCK_SIZE **************************/
#define BASE_BSIZE		8
#define QUADTREE_DEPTH	4	//可変ブロックサイズにするときは -f フラグを取る
#define MAX_BSIZE		32
#define MIN_BSIZE		(MAX_BSIZE >> QUADTREE_DEPTH)

/***** CLASS *******************************/
#define NUM_CLASS 		-1	//負の値にすると画像サイズに依存したクラス数
#define MAX_CLASS		63

/***** PREDICTOR ***************************/
#define COEF_PRECISION  6
#define PRD_ORDER       -1
#define BASE_PRD_ORDER  20
#define MAX_PRD_ORDER   110

/***** GROUP *******************************/
#define NUM_GROUP       16

/***** UPARA *******************************/
#define MAX_UPARA       512
#define UPEL_DIST       3
#define NUM_UPELS       (UPEL_DIST * (UPEL_DIST + 1))
#define CTX_WEIGHT		1	// 1 : weight on (standard) , 0 : weight off
#if CTX_WEIGHT
	#define MHD_WEIGHT	1	// 1 : Manhattan Distance (standard), 0 : Euclid Distance
#endif

/***** PMODEL ******************************/
#define PM_ACCURACY		3
#define NUM_PMODEL 		16
#define MIN_FREQ		1
#define PMCLASS_MAX		16
#define PMCLASS_LEVEL	32
#define PMMASK_MAX		5
#define PMMASK_LEVEL	10
#define NUM_ZMODEL		49
#define TOT_ZEROFR		(1 << 10)
#define MAX_SYMBOL		1024	// must be >> MAX_UPARA

/***** RangeCoder **************************/
#define HAVE_64BIT_INTEGER	1
#if HAVE_64BIT_INTEGER
	#define RANGE_SIZE 		64
	#if defined(__INTEL_COMPILER) || defined(_MSC_VER) || defined(__BORLANDC__)
		#define range_t unsigned __int64
	#else
		#define range_t unsigned long long int
	#endif
	#define MAX_TOTFREQ (1 << 20)	/* must be < RANGE_BOT */
#else
	#define RANGE_SIZE 32
	#define range_t unsigned int
	#define MAX_TOTFREQ (1 << 14)	/* must be < RANGE_BOT */
#endif
#define RANGE_TOP  ((range_t)1 << (RANGE_SIZE - 8))
#define RANGE_BOT  ((range_t)1 << (RANGE_SIZE - 16))

/***** TYPE DEFINE *************************/
#define uint            unsigned int
#define img_t           unsigned char
#define cost_t          double

/***** PI **********************************/
#ifndef M_PI
#	define M_PI 3.14159265358979323846
#endif

/***** MACRO DEFINE ************************/
#define CLIP(min, max, i)       ((i < min) ? min : ((i > max) ? max : i))

/***** TIME ********************************/
#define HAVE_CLOCK

/*****TEMPLATE MATCHING ********************/
#define TEMPLATE_MATCHING_ON 	1
#if TEMPLATE_MATCHING_ON

#define ZNCC			0

#if ZNCC
#define	ZSAD			0
#else
#define ZSAD			1
#endif

#define MANHATTAN_SORT	0	//市街地距離で近い順に事例を更に並び替える
#define AREA			6
#define Y_SIZE			20
#define X_SIZE			20	//調査結果より80*20がいいかも？
#define NAS_ACCURACY	1000
#define MAX_DATA_SAVE	50
#define MAX_DATA_SAVE_DOUBLE 	MAX_DATA_SAVE*4
#define W_GR 			7
#define WEIGHT_CN		2	//ラプラス関数
#define TEMPLATE_CLASS_NUM	30
#define TEMPLATEM_LOG_OUTPUT	1	//テンプレートマッチングの結果を書き出す

#define PREDICT_FROM_EXAMPLES	1

#else
#define TEMPLATE_CLASS_NUM	0

#endif

#define TEMPLATE_FLAG	2 << COEF_PRECISION

/*********DEBUG******************************/
#define CHECK_TM 		0
#define CHECK_TM_DETAIL	0
#define CHECK_DEBUG 	0
#define CHECK_PMODEL	0
#define CHECK_CLASS		0
#define CHECK_PREDICTOR	0
#define check_y			26
#define check_x			202
#define F_NUM			8

/****************** PARARELL ***************/
#define  OPENMP_ON		0
#define NUM_THREADS		8	//depend on server threads num

/***** STRUCTURE ***************************/

typedef struct {
	int height;
	int width;
	int maxval;
	img_t **val;
} IMAGE;

typedef struct {
	int size;
	int id;
	uint *freq;
	uint *cumfreq;
	float *cost;
	float *subcost;
	double norm;
} PMODEL;

typedef struct {
	char num_peak;
	char *class;
	int *base;
	PMODEL **pm;
	int *weight;
	int freq;
	int cumfreq;
#if TEMPLATE_MATCHING_ON
	char temp_cl;
#endif
} MASK;

typedef struct {
	range_t low;
	range_t code;
	range_t range;
	int y, x;
} RANGECODER;

typedef struct {
	int y, x;
} CPOINT;

typedef struct{
	int id;
	int by;
	int bx;
	int sum;
	int ave_o;

#if MANHATTAN_SORT
	int mhd;
#endif

#if ZNCC
	double s_devian;
#endif

} TM_Member;

typedef struct{
	int num_class_s;
	int **th_s;
	int **prd_s;
	char **uq_s;
	int ***prd_cl_s;
	char **class;
} RESTORE_SIDE;

typedef struct {
	int height;
	int width;
	int maxval;
	int num_class;
	int num_group;
	int prd_order;
	int max_prd_order;
	int coef_precision;
	int num_pmodel;
	int pm_accuracy;
	int maxprd;
	int max_coef;
	int quadtree_depth;
	int optimize_loop;
	int **predictor;
	int **nzconv;
	int *num_nzcoef;
	int **th;
	int **upara;
	int **prd;
	int **encval;
	int **err;
	int **org;
	int *ctx_weight;
	int ***roff;
	int ***prd_class;
	int ***weight;	//マスクを用いた確率モデルの高さの重み
	char **mask;
	int qtctx[QUADTREE_DEPTH << 3];
	char **qtmap[QUADTREE_DEPTH];
	char **class;
	char **group;
	char **uquant;
	int **econv;
	img_t *bconv;
	img_t *fconv;
	PMODEL ***pmodels;
	PMODEL **pmlist;
	PMODEL spm;
	PMODEL mult_pm;
	RANGECODER *rc;
	double *sigma;
	int *mtfbuf;
	int *cl_hist;
	uint **cl_pos;
	int *coef_m;
#if AUTO_PRD_ORDER
	int prd_mhd;
	int *ord2mhd;
	int *num_search;
	int *zero_m;
	int *zero_fr;
	cost_t ***coef_cost;
#else
	cost_t **coef_cost;
#endif
	cost_t *th_cost;
	cost_t *class_cost;
	cost_t qtflag_cost[QUADTREE_DEPTH << 3];
#if AUTO_DEL_CL
	cost_t **err_cost;
#endif
#if TEMPLATE_MATCHING_ON
	int **temp_num;
	int ***array;
	char temp_peak_num;
	int *w_gr;
	int ***temp_roff;
#endif
	char temp_cl;
	char function_number;
} ENCODER;

typedef struct {
	int version;
	int height;
	int width;
	int maxval;
	int num_class;
	int num_group;
	int max_prd_order;
	int num_pmodel;
	int pm_accuracy;
	int maxprd;
	int max_coef;
	int coef_precision;
	int quadtree_depth;
	int **predictor;
	int **nzconv;
	int *num_nzcoef;
	int **th;
	char **mask;
	int **err;
	int ***roff;
	int **econv;
	int *ctx_weight;
	char **qtmap[QUADTREE_DEPTH];
	char **class;
	int *pm_idx;
	PMODEL ***pmodels;
	PMODEL spm;
				PMODEL mult_pm;
	RANGECODER *rc;
	double *sigma;
	int *mtfbuf;
#if AUTO_PRD_ORDER
	int prd_mhd;
	int *zero_fr;
	int *ord2mhd;
#endif
#if TEMPLATE_MATCHING_ON
	int **temp_num;
	int *array;
	char temp_peak_num;
	int *w_gr;
	int ***temp_roff;
#endif
	char temp_cl;
	int **org;
} DECODER;

/***** FUNC - common.c ****************************/
FILE *fileopen(char *, char *);
void *alloc_mem(size_t);
void **alloc_2d_array(int, int, int);
void ***alloc_3d_array(int, int, int, int);
IMAGE *alloc_image(int, int, int);
PMODEL ***init_pmodels(int, int, int, int *, double *, int);
void printmodel(PMODEL *, int);
void set_pmodel_mult(PMODEL *, MASK *, int);
void set_spmodel(PMODEL *, int, int);
int *init_ctx_weight(void);
void mtf_classlabel(char **, int *, int, int, int, int, int);
double cpu_time(void);
void init_array(int *, int , int);
void init_2d_array(int **, int , int, int);
void init_3d_array(int ***, int , int, int, int);
int cmp(const void*, const void*);
int round_int(double);

/***** FUNC - rc.c ********************************/
RANGECODER *rc_init(void);
void rc_encode(FILE *, RANGECODER *, uint, uint, uint);
void rc_finishenc(FILE *, RANGECODER *);
int rc_decode(FILE *, RANGECODER *, PMODEL *, int, int);
void rc_startdec(FILE *, RANGECODER *);

/***** FUNC - log.c *******************************/
void print_predictor(int **, int, int, int, char *);
void print_threshold(int **, int, int, PMODEL **, int *, int*,  char *);
void print_class(char **, int, int, int, char *);
void output_class_map(char **, int, int, int, char *);
void print_class_color(char **, int, int, int, char *);
void print_class_and_block(char **, int, char ***, int, int, int, char *);
void print_mask(char **, int, int, char *);
void print_amp_chara(int **, int, int, int, int, char *);
//Lower funcs are can use in encoder only.
void print_rate_map(ENCODER *, char *);
void print_rate_compare_map(ENCODER *, char *);
void print_rate_compare_class_map(ENCODER *, char *);
void print_block_size(int **, char ***, int, int, int, char *);
void calc_var_upara( ENCODER *, char *);
void init_log_sheet(ENCODER *, char *);
void finish_log_sheet(ENCODER *, int, int, int, int, int, int, double, double);
void TemplateM_Log_Output(ENCODER *, char *, int ***, int ***);
void TemplateM_Log_Input(ENCODER *, char *, int ***, int ***);

#if defined(_WIN32)
	int set_directory(void);
#endif

