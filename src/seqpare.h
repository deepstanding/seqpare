//=====================================================================================
//Common structs, parameters, functions
//Based on AIList (John Feng, 2018)
//by Selena Feng  12/05/2019--1/25/2020
//-------------------------------------------------------------------------------------
#ifndef __SEQPARE_H__
#define __SEQPARE_H__
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <glob.h>
#include <zlib.h>
#include <assert.h>
//-------------------------------------------------------------------------------------
#include "khash.h"
#include "kseq.h"
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAXC 10							//max number of components
//-------------------------------------------------------------------------------------
typedef struct {
	int32_t N1;
	int32_t N2;
	int64_t tc;	//total counts
	double teo;	//total effective overlaps
} seqpare_t;

typedef struct {						//as database
    uint32_t start;      				//region start: 0-based
    uint32_t end;    					//region end: not inclusive
    int32_t q_idx;						//index of the query interval
    int32_t cnts;						//intersection counts
    float s_max;						//similarity index
} gdata_t;

typedef struct {						//as query
    uint32_t start;      				//region start: 0-based
    uint32_t end;    					//region end: not inclusive
} gdata0_t;

typedef struct{
	char *name;    						//name of the contig
	int64_t nr, mr;						//number of regions
	gdata_t *glist;						//regions data
	int nc, lenC[MAXC], idxC[MAXC];		//components
	uint32_t *maxE;						//augmentation	
} ctg_t;

typedef struct{
	char *name;    						//name of the contig
	int64_t nr, mr;						//number of regions
	gdata0_t *glist;						//regions data
	int nc, lenC[MAXC], idxC[MAXC];		//components
	uint32_t *maxE;						//augmentation	
} ctg0_t;

typedef struct {	
	ctg_t *ctg;        					// list of contigs (of size _n_ctg_)
	int32_t nctg, mctg; 				// number and max number of contigs
	void *hc;             				// dict for converting contig names to int    
} ailist_t;

typedef struct {	
	ctg0_t *ctg;        					// list of contigs (of size _n_ctg_)
	int32_t nctg, mctg; 				// number and max number of contigs
	void *hc;             				// dict for converting contig names to int    
} ailist0_t;

//-------------------------------------------------------------------------------------
//Parse a line of BED file
char *parse_bed(char *s, int32_t *st_, int32_t *en_);

//Initialize ailist_t
ailist_t *ailist_init(void);
ailist0_t *ailist0_init(void);

//read .BED file
ailist_t* readBED(const char* fn);
ailist0_t* readBED0(const char* fn);

//Add a gdata_t interval
void ailist_add(ailist_t *ail, const char *chr, uint32_t s, uint32_t e, int32_t v);
void ailist0_add(ailist0_t *ail, const char *chr, uint32_t s, uint32_t e, int32_t v);

//Construct ailist: decomposition and augmentation
void ailist_construct(ailist_t *ail, int cLen);
void ailist_construct0(ailist_t *ail, int cLen);
void ailist0_construct(ailist0_t *ail, int cLen);

//Get chr index
int32_t get_ctg(const ailist_t *ail, const char *chr);

//Binary search
uint32_t bSearch(gdata_t* As, uint32_t idxS, uint32_t idxE, uint32_t qe);

//Query ailist intervals
uint32_t ailist_query(ailist_t *ail, char *chr, uint32_t qs, uint32_t qe, uint32_t *mr, uint32_t **ir);

//Query ailist intervals: 1 to 1
void seq_compare_11(char *file1, char *file2, seqpare_t *sm);

//Query ailist intervals: 1 to n (folder of n .bed fiels)
void seq_compare_1n(char *iPath, char *qfile, int32_t *nf, seqpare_t **sps);
void seq_compare_mn(char *iPath1, char *iPath2, int32_t *nf, int32_t *mf, seqpare_t ***sps);
void seq_compare_nn(char *iPath, int32_t *nf, seqpare_t ***sps);

//compare two lists
void seq_compare(ailist_t *ail1, ailist0_t *ail2, seqpare_t *sp);
void seqpare0(ailist_t *ail1, ailist0_t *ail2, int gid1, int gid2, int iq);

//Free ailist data
void ailist_destroy(ailist_t *ail);
void ailist0_destroy(ailist0_t *ail);
//-------------------------------------------------------------------------------------
//The following section taken from Dr Heng Li's cgranges
// (https://github.com/lh3/cgranges)

KSTREAM_INIT(gzFile, gzread, 0x10000)
/**************
 * Radix sort *
 **************/
#define RS_MIN_SIZE 64
#define RS_MAX_BITS 8

#define KRADIX_SORT_INIT(name, rstype_t, rskey, sizeof_key) \
	typedef struct { \
		rstype_t *b, *e; \
	} rsbucket_##name##_t; \
	void rs_insertsort_##name(rstype_t *beg, rstype_t *end) \
	{ \
		rstype_t *i; \
		for (i = beg + 1; i < end; ++i) \
			if (rskey(*i) < rskey(*(i - 1))) { \
				rstype_t *j, tmp = *i; \
				for (j = i; j > beg && rskey(tmp) < rskey(*(j-1)); --j) \
					*j = *(j - 1); \
				*j = tmp; \
			} \
	} \
	void rs_sort_##name(rstype_t *beg, rstype_t *end, int n_bits, int s) \
	{ \
		rstype_t *i; \
		int size = 1<<n_bits, m = size - 1; \
		rsbucket_##name##_t *k, b[1<<RS_MAX_BITS], *be = b + size; \
		assert(n_bits <= RS_MAX_BITS); \
		for (k = b; k != be; ++k) k->b = k->e = beg; \
		for (i = beg; i != end; ++i) ++b[rskey(*i)>>s&m].e; \
		for (k = b + 1; k != be; ++k) \
			k->e += (k-1)->e - beg, k->b = (k-1)->e; \
		for (k = b; k != be;) { \
			if (k->b != k->e) { \
				rsbucket_##name##_t *l; \
				if ((l = b + (rskey(*k->b)>>s&m)) != k) { \
					rstype_t tmp = *k->b, swap; \
					do { \
						swap = tmp; tmp = *l->b; *l->b++ = swap; \
						l = b + (rskey(tmp)>>s&m); \
					} while (l != k); \
					*k->b++ = tmp; \
				} else ++k->b; \
			} else ++k; \
		} \
		for (b->b = beg, k = b + 1; k != be; ++k) k->b = (k-1)->e; \
		if (s) { \
			s = s > n_bits? s - n_bits : 0; \
			for (k = b; k != be; ++k) \
				if (k->e - k->b > RS_MIN_SIZE) rs_sort_##name(k->b, k->e, n_bits, s); \
				else if (k->e - k->b > 1) rs_insertsort_##name(k->b, k->e); \
		} \
	} \
	void radix_sort_##name(rstype_t *beg, rstype_t *end) \
	{ \
		if (end - beg <= RS_MIN_SIZE) rs_insertsort_##name(beg, end); \
		else rs_sort_##name(beg, end, RS_MAX_BITS, (sizeof_key - 1) * RS_MAX_BITS); \
	}

/*********************
 * Convenient macros *
 *********************/

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#define CALLOC(type, len) ((type*)calloc((len), sizeof(type)))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

#define EXPAND(a, m) do { \
		(m) = (m)? (m) + ((m)>>1) : 16; \
		REALLOC((a), (m)); \
	}while (0) 

#endif
