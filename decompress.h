#ifndef _DECOMPRESS_H
#define _DECOMPRESS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <stddef.h>
#include <stdlib.h>

#define LINESIZE      100000000//4096
#define N_LINK        3
#define BLOCK         2000000

#include "compress.h"

extern size_t read_len;

char *get_fastq_name(char *fastq_name,char *str);
void init_decompression(char *fastq_name,char *ref_name_);
void FQZip_decompression(char *fastq_name, char *compressedfile,int assemble_flag,int highestCom_flag);
void recover_base(char *ref_name_);
void recover_meta(char *fastq_name);
char *variation_recover(char *base_temp,char *_cor);
char *palindrome_conversion(char *des,char *str);
void recover_fastq(char *fastq_name);
char *get_meta_variation_len(char *len,char *meta);
char *get_meta_variation(char *variation,char *meta);
void quit_decompression_dispose(int assemble_flag);
void LWFQZip_decompression(char *input_name,char *ref_name,int assemble_flag,int highestCom_flag);
char* decode_quality_score(char *Instream, FILE *outfp);
void recover_qs();
#endif // _DECOMPRESS_H
