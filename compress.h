#ifndef _COMPRESS_H
#define _COMPRESS_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <stddef.h>
#include <stdlib.h>
#include <sys/stat.h>
#define block_size 2000000//200M
//#define REFERENCE_SIZE 1*1024*1024*1024  //1G

#define OFFSET 25

#define LINESIZE       1000000
#define N_LINK         2   //set 2 can be smaller

#define BITSET(a,b) ( a|=(1<<b) )  //set 1
#define BITCLK(a,b) ( (a)&=~(1<<(b)) )  //set 0
#define BITTEST(a,b) ( (a)&(1<<(b)) )  //test 1 or 0

extern size_t read_len;
/* The two types:
0 : read only "rb"
1 : write only "wb"
 */


FILE *m_fopen(char *file_name,int types);
void* m_malloc(size_t n);
void init_compression(char *fastq_name,char *ref_name);
size_t get_file_size(char *file);
char* get_meta(char *meta_temp,char *metadta);
void get_meta_qscore(char *file);
void quit_compression_dispose();
char m_toupper(char ch);
void split_map_result(char *ref_name);
char *get_pos(char *pos_temp,char *map_result);
char *get_flag(char *flag,char *map_result);
char *get_seq(char *seq_temp,char *map_result);
char *get_cigar(char *samcigar,char *map_result);
char *get_mismatch(char *match,const char *seq,const char *base,char *pos);
char *new_get_base(char *base,char *add_,char *samcigar,char *ref_gen,char *pos,char *base_);
char *new_split_cigar(char *add_temp,char *seq,char *samcigar,char *flag,char *clip_seq);
char *m_strcat(char *dst, const char* src);
void LWFQZip_compression(char *input_name,char *ref_name,int assemble_flag,unsigned long long totalread,int BlockNum,int highestCom_flag);
int rll_quality_score(char *buf, char *str);
#endif // _COMPRESS_H
