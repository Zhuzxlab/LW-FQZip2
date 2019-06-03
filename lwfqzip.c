#include "lwfqzip.h"
#include <pthread.h>

#define ReadLen 300000
#define MaxBlockNum 50
#define _FILE_OFFSET_BITS 64

int BlockNum = 10; //Default file block && number of thread
unsigned long long totleread=0,mappedread=0,exactmatch=0,doublemap=0,totalbase=0, mappedbase=0;
char input_name[5000]={'\0'},ref_name[5000]={'\0'};

struct M
{
	char Minput[500];
	char MOutFile[MaxBlockNum][10000];
}M1;

struct G
{
	char input[500];
	char orders[10000];
}G1[MaxBlockNum];

struct P
{
	char filename[500];
	char file[500];
	long int readpos;
	long int readle;
	FILE *subfastq;
	FILE *fastq;

}P1[MaxBlockNum];

void *map_thread(void * args)
{
	struct G *G2 = args;
	system(G2->orders);
	return NULL;
}

void *block_thread(void * args)
{
	struct P *P2 = args;
	int subline = 0;
	char tempStr1[ReadLen] = { '\0' }; //temporary row number
	P2->subfastq = fopen(P2->filename, "wb");
	if ((P2->fastq = fopen(P2->file, "rb")) == NULL) { printf("can not find this file:%s\n", P2->file); exit(0); }
	fseek(P2->fastq,P2->readpos,SEEK_SET);
	while (fgets(tempStr1, ReadLen,P2->fastq) != NULL){
		if (subline  == ((P2->readle))*4)
			break;
		subline++;
		fputs(tempStr1, P2->subfastq);
	}
	fclose(P2->subfastq);
	fclose(P2->fastq);
	return NULL;
}

void prekeycontent(char  *refname)
{
	char prekeyorder[10000] = {0};
	char prekeyorderAT[10000] = {0},prekeyorderCG[10000] = {0},prekeyorderGT[10000] = {0},prekeyorderAC[10000] = {0};
	char prekeyorderTAG[10000] = {0},prekeyorderGCA[10000] = {0},prekeyorderTAT[10000] = {0},prekeyorderACT[10000] = {0};
	char ref_name[500] = {'\0'};

	strcpy(ref_name, refname);

	sprintf(prekeyorderAT,"%s%s%s%s%s%s","grep -o ","'","AT","' ",ref_name," | wc -l ");
	sprintf(prekeyorderCG,"%s%s%s%s%s%s","grep -o ","'","CG","' ",ref_name," | wc -l ");
	sprintf(prekeyorderGT,"%s%s%s%s%s%s","grep -o ","'","GT","' ",ref_name," | wc -l ");
	sprintf(prekeyorderAC,"%s%s%s%s%s%s","grep -o ","'","AC","' ",ref_name," | wc -l ");
	sprintf(prekeyorderTAG,"%s%s%s%s%s%s","grep -o ","'","TAG","' ",ref_name," | wc -l ");
	sprintf(prekeyorderGCA,"%s%s%s%s%s%s","grep -o ","'","GCA","' ",ref_name," | wc -l ");
	sprintf(prekeyorderTAT,"%s%s%s%s%s%s","grep -o ","'","TAT","' ",ref_name," | wc -l ");
	sprintf(prekeyorderACT,"%s%s%s%s%s%s","grep -o ","'","ACT","' ",ref_name," | wc -l ");

	printf("The counts of prekey AT is:\n");
	system(prekeyorderAT);
	printf("The counts of prekey CG is:\n");
	system(prekeyorderCG);
	printf("The counts of prekey GT is:\n");
	system(prekeyorderGT);
	printf("The counts of prekey AC is:\n");
	system(prekeyorderAC);
	printf("The counts of prekey TAG is:\n");
	system(prekeyorderTAG);
	printf("The counts of prekey GCA is:\n");
	system(prekeyorderGCA);
	printf("The counts of prekey TAT is:\n");
	system(prekeyorderTAT);
	printf("The counts of prekey ACT is:\n");
	system(prekeyorderACT);
}

void Output_print(int block){
	
	int i,j,b;
	b=block;
	FILE *readout=NULL;
	unsigned long long digit[b][6];
	if ((readout = fopen("Output.txt","rb")) == NULL) {  //whether the file exists
		printf("can not find this file:%s\n!", "Output.txt");
		exit(0);
	}
	i=0;
	while(i<b){ fscanf(readout, "%lld,%lld,%lld,%lld,%lld,%lld",&digit[i][0],&digit[i][1],&digit[i][2],&digit[i][3],&digit[i][4],&digit[i][5]);   i++; }
	for(j=0;j<b;j++){ totleread+=digit[j][0]; mappedread+=digit[j][1]; exactmatch+=digit[j][2];doublemap+=digit[j][3]; totalbase+=digit[j][4]; mappedbase+=digit[j][5]; }
	fclose(readout);

	printf("Total readCounts: %lld\n",totleread);
	printf("Mapped reads= %d\n",mappedread);
	printf("Unmapped read= %d\n",totleread-mappedread);
	printf("Mappedread rate: %f\n",(float)mappedread/totleread);
	//printf("Exact Match= %d\n",exactmatch);
	//printf("Inexact Match= %d\n",mappedread-exactmatch);
	//printf("N*= %d\n",doublemap);
	//printf("Total base:%ld\nMapped base:%ld\nMapped base rate:%f\n", totalbase, mappedbase,(float)mappedbase / totalbase);
	remove("Output.txt");
}
/*
void Merge_file(char  *inputname,char  *refname,char outfile[][10000]){
	
	int k;
	char pos_name[500]={'\0'},cigar_name[500]={'\0'},cor_name[500]={'\0'},add_name[500]={'\0'},pos_order[5000]={'\0'},cigar_order[5000]={'\0'},cor_order[5000]={'\0'},add_order[5000]={'\0'};
	char pos_outfile[MaxBlockNum][10000] = {0},cigar_outfile[MaxBlockNum][10000] = {0},cor_outfile[MaxBlockNum][10000] = {0},add_outfile[MaxBlockNum][10000] = {0};
	FILE *pos_w=NULL,*cigar_w=NULL,*add_w=NULL,*cor_w=NULL;
	
	strcpy(pos_name,inputname); strcat(pos_name,".pos.txt");
	strcpy(cigar_name,inputname); strcat(cigar_name,".cigar.txt");
	strcpy(cor_name,inputname); strcat(cor_name,".cor.txt");
	strcpy(add_name,inputname); strcat(add_name,".add.txt");
	strcpy(pos_order, "cat");
	strcpy(cigar_order, "cat");
	strcpy(cor_order, "cat");
	strcpy(add_order, "cat");
	
	pos_w=fopen(pos_name,"wb");
	cor_w=fopen(cor_name,"wb");
	cigar_w=fopen(cigar_name,"wb");
	add_w=fopen(add_name,"wb");

	for (k = 0; k < BlockNum; k++) {
		char s[100];
		sprintf(s, "%d", k);
		strcpy(pos_outfile[k], outfile[k]); strcat(pos_outfile[k],".pos.txt");
		strcpy(cigar_outfile[k], outfile[k]); strcat(cigar_outfile[k], ".cigar.txt");
		strcpy(cor_outfile[k], outfile[k]); strcat(cor_outfile[k], ".cor.txt");
		strcpy(add_outfile[k], outfile[k]); strcat(add_outfile[k], ".add.txt");
		
		strcat(pos_order," "); strcat(pos_order,pos_outfile[k]);
		strcat(cigar_order," "); strcat(cigar_order,cigar_outfile[k]);
		strcat(cor_order," "); strcat(cor_order,cor_outfile[k]);
		strcat(add_order," "); strcat(add_order,add_outfile[k]);
	}

	strcat(pos_order," > "); strcat(pos_order,pos_name);
	strcat(cigar_order," > "); strcat(cigar_order,cigar_name);
	strcat(cor_order," > "); strcat(cor_order,cor_name);		
	strcat(add_order," > "); strcat(add_order,add_name);
	
	system(pos_order);
	system(cigar_order);
	system(cor_order);
	system(add_order);
	for (k = 0; k < BlockNum; k++) {
		remove(pos_outfile[k]);
		remove(cigar_outfile[k]);
		remove(cor_outfile[k]);
		remove(add_outfile[k]);
	}
	fclose(pos_w);
	fclose(cigar_w);
	fclose(add_w);
	fclose(cor_w);
}
*/
void *Merge_pos_file(void * args){
	
	struct M *Mp = args;
	int k;
	char pos_name[500]={'\0'},pos_order[5000]={'\0'};
	char pos_outfile[MaxBlockNum][10000] = {0};
	
	strcpy(pos_name,Mp->Minput); strcat(pos_name,".pos.txt");
	strcpy(pos_order, "cat");
	
	for (k = 0; k < BlockNum; k++) {
		char s[100];
		sprintf(s, "%d", k);
		strcpy(pos_outfile[k], Mp->MOutFile[k]); strcat(pos_outfile[k],".pos.txt");
		strcat(pos_order," "); strcat(pos_order,pos_outfile[k]);
	}

	strcat(pos_order," > "); strcat(pos_order,pos_name);
	system(pos_order);
	for (k = 0; k < BlockNum; k++) {
		remove(pos_outfile[k]);
	}
}

void *Merge_cigar_file(void * args){
	
	struct M *Mc = args;
	int k;
	char cigar_name[500]={'\0'},cigar_order[5000]={'\0'};
	char cigar_outfile[MaxBlockNum][10000] = {0};
	
	strcpy(cigar_name,Mc->Minput); strcat(cigar_name,".cigar.txt");
	strcpy(cigar_order, "cat");

	for (k = 0; k < BlockNum; k++) {
		char s[100];
		sprintf(s, "%d", k);	
		strcpy(cigar_outfile[k], Mc->MOutFile[k]); strcat(cigar_outfile[k], ".cigar.txt");		
		strcat(cigar_order," "); strcat(cigar_order,cigar_outfile[k]);
	}
	strcat(cigar_order," > "); strcat(cigar_order,cigar_name);
	system(cigar_order);

	for (k = 0; k < BlockNum; k++) {
		remove(cigar_outfile[k]);
	}
}

void *Merge_cor_file(void * args){
	
	struct M *Mr = args;
	int k;
	char cor_name[500]={'\0'},cor_order[5000]={'\0'};
	char cor_outfile[MaxBlockNum][10000] = {0};
		
	strcpy(cor_name,Mr->Minput); strcat(cor_name,".cor.txt");
	strcpy(cor_order, "cat");
	
	for (k = 0; k < BlockNum; k++) {
		char s[100];
		sprintf(s, "%d", k);
		strcpy(cor_outfile[k], Mr->MOutFile[k]); strcat(cor_outfile[k], ".cor.txt");		
		strcat(cor_order," "); strcat(cor_order,cor_outfile[k]);
	}

	strcat(cor_order," > "); strcat(cor_order,cor_name);		
	
	system(cor_order);
	for (k = 0; k < BlockNum; k++) {
		remove(cor_outfile[k]);
	}
}


void *Merge_add_file(void * args){
	
	struct M *Ma = args;
	int k;
	char add_name[500]={'\0'},add_order[5000]={'\0'};
	char add_outfile[MaxBlockNum][10000] = {0};

	strcpy(add_name,Ma->Minput); strcat(add_name,".add.txt");
	strcpy(add_order, "cat");
	

	for (k = 0; k < BlockNum; k++) {
		char s[100];
		sprintf(s, "%d", k);
		strcpy(add_outfile[k], Ma->MOutFile[k]); strcat(add_outfile[k], ".add.txt");
		strcat(add_order," "); strcat(add_order,add_outfile[k]);
	}
		
	strcat(add_order," > "); strcat(add_order,add_name);
	
	system(add_order);
	for (k = 0; k < BlockNum; k++) {
		remove(add_outfile[k]);
	}
}

void *Merge_qs_file(void * args){
	
	struct M *Mq = args;
	int k;
	char qs_name[500]={'\0'},qs_order[5000]={'\0'};
	char qs_outfile[MaxBlockNum][10000] = {0};
		
	strcpy(qs_name,Mq->Minput); strcat(qs_name,".qs.txt");
	strcpy(qs_order, "cat");

	for (k = 0; k < BlockNum; k++) {
		char s[100];
		sprintf(s, "%d", k);
		strcpy(qs_outfile[k], Mq->MOutFile[k]); strcat(qs_outfile[k], ".qs.txt");		
		strcat(qs_order," "); strcat(qs_order,qs_outfile[k]);
	}

	strcat(qs_order," > "); strcat(qs_order,qs_name);		
	
	system(qs_order);
	for (k = 0; k < BlockNum; k++) {
		remove(qs_outfile[k]);
	}
}
/*
void *Merge_meta_file(void * args){
	
	struct M *Mm = args;
	int k,i,l;
	char meta_name[500]={'\0'},meta_order[5000]={'\0'},meta_info[ReadLen]={'\0'};
	char meta_outfile[MaxBlockNum][10000] = {0};
	FILE *meta_w=NULL,*meta_r=NULL;
	size_t buffer_sizes=get_file_size(Mm->Minput),buf_counts=0;
	if(buffer_sizes<block_size)
		buf_counts = buffer_sizes;
	else buf_counts = block_size;
	char *_buffer=m_malloc(buf_counts+ReadLen);
	char *meta_temp=m_malloc(ReadLen);
	memset(_buffer,'\0',sizeof(char)*(buf_counts+ReadLen));
	memset(meta_temp,'\0',sizeof(char)*ReadLen);
		
	strcpy(meta_name,Mm->Minput); strcat(meta_name,".meta.txt");
	strcpy(meta_order, "cat");
	
	meta_r=fopen(Mm->Minput,"rb");
	meta_w=fopen(meta_name,"wb");
	
	fgets(meta_info,ReadLen,meta_r);
	strcpy(_buffer,meta_info);
	//get_meta(meta_mid,_buffer);
	//sscanf(_buffer,"%*s%s",meta_mid);
	for(i=0,l=0;i<strlen(_buffer);i++,l++)
	{
		if(_buffer[i]==' '&&_buffer[i+1]=='l'&&_buffer[i+2]=='e')//remove "length="
			break;
		else *(meta_temp+l)=*(_buffer+i);
	}
	fputs(meta_temp,meta_w);
	fputc('\n',meta_w);
	fclose(meta_w);
	fclose(meta_r);
	
	for (k = 0; k < BlockNum; k++) {
		char s[100];
		sprintf(s, "%d", k);
		strcpy(meta_outfile[k], Mm->MOutFile[k]); strcat(meta_outfile[k], ".meta.txt");		
		strcat(meta_order," "); strcat(meta_order,meta_outfile[k]);
	}

	strcat(meta_order," >> "); strcat(meta_order,meta_name);		
	
	system(meta_order);
	for (k = 0; k < BlockNum; k++) {
		remove(meta_outfile[k]);
	}
}
*/

long int file_size(char *filename){
	
	struct stat buf;
	stat (filename, &buf);
	long int filesize;
	filesize = buf.st_size;
	return filesize;
}

long long ifnormal(char *inputname,long int pos){//修改....chx
	
	FILE *fastqr = NULL;
	//int countline =0;
	//char a;
	char tempStr[ReadLen] = { '\0' },readtmp[500] = { '\0' },readline[500] = { '\0' };
	long long linecounts = 0;
	if ((fastqr = fopen(inputname, "rb")) == NULL) {  printf("can not find this file:%s\n!", inputname); exit(0); }
	while(fgets(tempStr,ReadLen,fastqr)!=NULL){
		linecounts++;
	}
	linecounts = linecounts>>2;
	return linecounts;
}

void multireadfile(long linecounts,int BlockNum,char *FASTQName,long int filesize){	//......修改过....chx
				
	long readcounts, subreadcounts,subreadcounts1, lastblockcounts,filepos,adjustpos,adjustline = 0;
	char adjusttempStr[ReadLen] = { '\0' },adjusttempStr1[ReadLen] = { '\0' },adjusttempStr2[ReadLen] = { '\0' };
	char u,ua,a ;
	int poscount,i;
	FILE *fastqr = NULL;
	
	if ((fastqr = fopen(FASTQName, "rb")) == NULL) { printf("can not find this file:%s\n!", FASTQName); exit(0); }
	fgets(adjusttempStr1,ReadLen,fastqr);
	fclose(fastqr);
	
	readcounts = linecounts<<2;
	subreadcounts = readcounts / BlockNum;		//每块有多少行
	if(readcounts%BlockNum==0){ lastblockcounts = subreadcounts; }		//刚好平分到每个块
	else{ lastblockcounts = subreadcounts+(readcounts - subreadcounts*BlockNum); }	//没有整除的话最后多出一个块来
	if ((fastqr = fopen(FASTQName, "rb")) == NULL) { printf("can not find this file:%s\n!", FASTQName); exit(0); }

	P1[0].readpos = 0;			//readpos是这块开始的位置
	
	for(poscount=1;poscount<BlockNum;poscount++){
		while(fgets(adjusttempStr,ReadLen,fastqr)!=NULL){
			adjustline++;
			if(adjustline == subreadcounts){
				while(1){
					if((adjusttempStr[0]==adjusttempStr1[0])&&(adjusttempStr[1]==adjusttempStr1[1])&&(adjusttempStr[2]==adjusttempStr1[2])&&(adjusttempStr[3]==adjusttempStr1[3])&&adjusttempStr[4]==adjusttempStr1[4]){
						//printf("adjustline是:%d\n",adjustline);
						fgets(adjusttempStr,ReadLen,fastqr);
						fgets(adjusttempStr,ReadLen,fastqr);
						fgets(adjusttempStr,ReadLen,fastqr);
						adjustline = adjustline+3;
						P1[poscount].readpos = ftell(fastqr);
						P1[poscount-1].readle = adjustline>>2;
					//	printf("P1[%d].readle是：%d",poscount-1,P1[poscount-1].readle);
						adjustline = 0;
						break;
					}
					else{
						fgets(adjusttempStr,ReadLen,fastqr);
						adjustline++;
					}
				}
				break;
			}	
		}		
		//printf("P1【%d】的readpos传完\n",poscount);
		//printf("ftell是%d\n",P1[poscount].readpos);
		/*filepos=(filesize/BlockNum) * poscount;		//filesize是整个文件的大小；filepos是每块的起始位置
		adjustpos=filepos;				//adjustpos是作为调整后的位置
		while(adjustpos){
			adjustpos++;
			fseek(fastqr,adjustpos,SEEK_SET);
			if(a= fgetc(fastqr)=='\n'){
				int posline=0;
				fgets(adjusttempStr, ReadLen, fastqr);
				char* prefix_meta_1;			//找出元数据的头几个字符作为评判标准
				for(int i=0;i<5;i++){
					prefix_meta_1[i] = adjusttempStr[i];
				}
				if(strcmp(prefix_meta,prefix_meta_1) == 0){	//进入循环	strcmp判断两个char*是否相等，相等返回0
					posline = 1;
				}
				if(posline){	//posline没有实际作用，仅作为是否进入if的条件
					long int subline;
					long int tmpline;
					sscanf(adjusttempStr,"%[^ ]",adjusttempStr1);
					sscanf(adjusttempStr1,"%*[^.].%s",adjusttempStr2);
					adjustline= atoi(adjusttempStr2);
					while(adjustline != (subreadcounts * poscount+1)){
						if(adjustline > (subreadcounts * poscount+1)){
							tmpline = adjustline - (subreadcounts * poscount+1);
							subline = tmpline<<2;
							while(subline){
								adjustpos--;
								fseek(fastqr,adjustpos,SEEK_SET);
								if(u= fgetc(fastqr)=='\n') subline--;
							}
							fgets(adjusttempStr, ReadLen, fastqr);
							sscanf(adjusttempStr,"%[^ ]",adjusttempStr1);
							sscanf(adjusttempStr1,"%*[^.].%s",adjusttempStr2);
							adjustline= atoi(adjusttempStr2);
						}
						else{
							tmpline = (subreadcounts * poscount+1) - adjustline;
							subline = tmpline<<2;
							while(subline){
								adjustpos++;
								fseek(fastqr,adjustpos,SEEK_SET);
								if(u= fgetc(fastqr)=='\n') subline--;
							}
							fgets(adjusttempStr, ReadLen, fastqr);
							sscanf(adjusttempStr,"%[^ ]",adjusttempStr1);
							sscanf(adjusttempStr1,"%*[^.].%s",adjusttempStr2);
							adjustline= atoi(adjusttempStr2);
						}
					}
					P1[poscount].readpos = adjustpos+1;
					break;
				}
			}
		}*/
	}
	fclose(fastqr);
	//printf("multireadfile该函数结束!\n");
}

void readfile(long linecounts,int BlockNum,char *FASTQName,long int filesize,char outfile[][10000]){
	
	pthread_t thread[BlockNum];
	long readcounts, subreadcounts, lastblockcounts, subline;
	FILE *fastq = NULL;
	FILE *fastq2 = NULL;
	FILE *subfastq[BlockNum];
	char tempStr[ReadLen] = { '\0' };
	int i;
	char c;

	if ((fastq = fopen(FASTQName, "rb")) == NULL) {  printf("can not find this file:%s\n!", FASTQName); exit(0);}
	while ((c = fgetc(fastq)) != EOF) { if (c == '\n') linecounts++; }
	fclose(fastq);
	readcounts = linecounts>>2;
	subreadcounts = readcounts / BlockNum;
	if(readcounts%BlockNum==0){ lastblockcounts = subreadcounts; }
	else{ lastblockcounts = subreadcounts+(readcounts - subreadcounts*BlockNum); }
	if ((fastq2 = fopen(FASTQName, "rb")) == NULL) { printf("can not find this file:%s\n!", FASTQName); exit(0); }
	for (i = 0; i < BlockNum - 1; i++){
		subline = 0;
		subfastq[i] = fopen(outfile[i], "wb");
		while (fgets(tempStr, ReadLen, fastq2) != NULL){
			subline++;
			fputs(tempStr, subfastq[i]);
			if (subline  == subreadcounts<<2) break;
		}
		fclose(subfastq[i]);
		pthread_create(&thread[i], 0, &map_thread, &G1[i]);
	}
	subline = 0;
	subfastq[BlockNum - 1] = fopen(outfile[BlockNum - 1], "wb");
	while (fgets(tempStr, ReadLen, fastq2) != NULL){
		if (subline == lastblockcounts<<2) break;
		subline++;
		fputs(tempStr, subfastq[BlockNum - 1]);
	}
	fclose(subfastq[BlockNum - 1]);
	pthread_create(&thread[BlockNum - 1], 0, &map_thread, &G1[BlockNum - 1]);
	fclose(fastq2);
	for (i = 0; i < BlockNum; i++) { pthread_join(thread[i], 0); }
}

void assemble_get_fasta(char *FASTQName,long int filesize,float FASTA_rate){
	
	FILE *OriFASTQ=NULL,*FASTA=NULL;
	char OriFASTQName[5000]={'\0'},tempStr[300000]={'\0'};
	float FASTA_Maxsize=filesize*FASTA_rate,FASTA_size=0;	
	long countline=0;				
	int i=0,flag_CG=0;		
	strcpy(OriFASTQName,FASTQName);
	strcpy(ref_name,FASTQName);
	strcat(ref_name,".fasta");
		
	if((OriFASTQ=fopen(OriFASTQName,"rb"))==NULL){
		printf("can not open file\n");
		exit(0);	
	}
	if((FASTA=fopen(ref_name,"wb"))==NULL){
		printf("can not open file\n");
		exit(0);
	}
			
	strcpy(tempStr,">gi|Assemble_ref|\n");
	fputs(tempStr,FASTA);

	while(fgets(tempStr,300000,OriFASTQ))
	{
		countline++;			
		if(FASTA_size+(strlen(tempStr)*sizeof(char))>FASTA_Maxsize)
		break;
		if(countline%4==2)
		{
			flag_CG=0;
			for(i=0;i<strlen(tempStr);i++)					
				if(tempStr[i]=='C'&&tempStr[i+1]=='G'){flag_CG=1;break;}
			if(flag_CG==1)
			{
				for(i=0;i<strlen(tempStr);i++)
				{
					if(tempStr[i]=='\n')
						continue;
					fputc(tempStr[i],FASTA);
				}
				FASTA_size+=(strlen(tempStr)*sizeof(char));
			}
		}
	}
	fclose(OriFASTQ);
	fclose(FASTA);
}

int main(int argc, char *argv[])
{
	bool compress_flag=false;
	bool decompress_flag=false;
	bool statistics_flag=false;
	int assemble_flag=0;
	int highestCom_flag=0;
	float FASTA_rate=0.003;//assemble-based
	char order[MaxBlockNum][10000] = {0};

	opterr = 0;
	int opt;
	while ((opt = getopt(argc, argv, "hvcdsi:gar:m:b:p:k:l:e:t:")) != -1)
	{
		switch (opt)
		{
			
			case 'I':
			case 'i':
				strcpy(input_name, argv[optind - 1]);
				if(!strstr(input_name,".fastq")&&!strstr(input_name,".fq"))
				{
				fprintf(stderr,"Please input correct FASTQ file\n");
				return 0;
				}
				break;
			case 'G':
			case 'g':
			highestCom_flag=1;//highest compression model
			printf("Start the best compression mode!\n");
				break;
			case 'R':
			case 'r':
				strcpy(ref_name, argv[optind-1]);
				break;
			case 'A':
			case 'a':
		    		assemble_flag=1;//assemble-based model
				printf("Start the assemble-based mode (default:0.003)\n");
			if(!decompress_flag)
			{
			if(argv[optind]!=NULL&&atof(argv[optind])!=0)
			{
				FASTA_rate=atof(argv[optind]);
				printf("FASTA extraction rate: %.2f%%\n",atof(argv[optind])*10000);
			}	
			}
				break;
			case 'C':
			case 'c'://compression
				compress_flag = true;
				break;
			case 'D':
			case 'd'://decompression
				decompress_flag = true;
				break;
			case 'S':
			case 's':
				statistics_flag = true;				
				break;
			case 'M':
			case 'm'://max_read_len
				read_len = strtoul(optarg, NULL, 100);
				if (read_len>300000 || read_len<30000)
				{
					fprintf(stderr, "Please input the correct read length.\n");
					fprintf(stderr, "Read length must be 30000--300000.\n");
					exit(0);
				}
			case 'B':
			case 'b':
				BlockNum = atoi(argv[optind - 1]);
				if (BlockNum < 6||BlockNum>64 )
				{
					fprintf(stderr, "the number of thread must be greater than 6 and less than 64\n"); exit(0);
				}
				break;
			case 'H':
			case 'h':
				print_help();
				exit(0);
				break;
			case 'V':
			case 'v':
				print_version();
				break;
			case '?':
				break;
			
			default:
				break;
		}
	}

	if(compress_flag&&decompress_flag)
	{
		fprintf(stderr,"Please input correct <Mode> parameter\n");
		fprintf(stderr,"Please choose only one <Mode> parameter\n");
		return 0;
	}
	if(compress_flag){

		int i,f;
		//temporary row number
		char tempStr[ReadLen] = { '\0' },OutFile[MaxBlockNum][10000] = { 0 },FASTQName[500] = { '\0' },prekeyname[5] = { '\0' };
		long long linecounts = 0; //line counts of file
		long readcounts, subreadcounts, lastblockcounts,filesize,gensize;

		//Sub-block
		strcpy(FASTQName, argv[3]);  //file name
		strcpy(M1.Minput,FASTQName);
		strcpy(prekeyname,"CG");
		FILE *fastqr = NULL;
		FILE *fastqw = NULL;
		long int pos = 1;
		
		//file size
		filesize=file_size(FASTQName);
		gensize=file_size(ref_name);
		printf ("File size (%s): %lld MB.\n",FASTQName, filesize>>20);
		if(!assemble_flag)
			printf ("File size (%s): %ld MB.\n",ref_name, gensize>>20);
		
		if(assemble_flag==1){
			assemble_get_fasta(FASTQName,filesize,FASTA_rate);			
		}
		
		if(statistics_flag){
			prekeycontent(ref_name);
			exit(0);
		}

		for (i = 0; i < BlockNum; i++){
			char s[100];
			sprintf(s, "%d", i);
			strcat(OutFile[i], FASTQName);
			strcat(M1.MOutFile[i], FASTQName);
			strcpy(P1[i].file, FASTQName);
			strcat(OutFile[i], s);
			strcat(M1.MOutFile[i], s);
			strcat(OutFile[i], ".fastq");
			strcat(M1.MOutFile[i], ".fastq");
			strcpy(G1[i].input, OutFile[i]);
			strcpy(P1[i].filename, OutFile[i]);
			strcat(order[i], "./LWMapping ");
			if(assemble_flag==1){	//assemble model
				strcat(order[i], FASTQName);
				strcat(order[i], ".fasta");
			}
			else					//reference-base model
				strcat(order[i], argv[5]);
			strcat(order[i], " ");
			strcat(order[i], G1[i].input);
			strcat(order[i], " ");
			strcpy(G1[i].orders, order[i]);
		}
		
		
		int pi;
		pi = 6;
		while (argc>pi){
			if (*argv[pi] == '-') { //option: -P -k -h
				int o;
				switch (*(argv[pi] + 1)) {
					case 't':
						pi++;
						break;
					case 'B':
					case 'b':
						pi++;
						break;
					case 'P':
					case 'p':
						pi++;
						for (o = 0; o < BlockNum; o++) {
							strcat(order[o], "-p ");
							strcat(order[o], argv[pi]);
							strcat(order[o], " ");
							strcpy(G1[o].orders, order[o]);
						}
						if (strlen(argv[pi]) < 2|| strlen(argv[pi]) > 3 )
						{
							fprintf(stderr, "the length of kmer prefixes must be greater than 1 and less than 4\n"); exit(0);
						}
						//strcpy(prekeyname,argv[pi]);
						break;
					case 'K':
					case 'k':
						pi++;
						for (o = 0; o < BlockNum; o++) {
							strcat(order[o], "-k ");
							strcat(order[o], argv[pi]);
							strcat(order[o], " ");
							strcpy(G1[o].orders, order[o]);
						}
						if (argv[pi] < 8 || argv[pi] > 12 )
						{
							fprintf(stderr, "the length of kmer must be greater than 7 and less than 13\n"); exit(0);
						}
						break;
					case 'E':
					case 'e':
						pi++;
						for (o = 0; o < BlockNum; o++) {
							strcat(order[o], "-e ");
							strcat(order[o], argv[pi]);
							strcat(order[o], " ");
							strcpy(G1[o].orders, order[o]);
						}
						if (atof(argv[pi]) < 0|| atof(argv[pi]) > 1 )
						{
							fprintf(stderr, "the tolerance rate of mismatches must be greater than 0 and less than 1\n"); exit(0);
						}
						break;
					case 'L':
					case 'l':
						pi++;
						for (o = 0; o < BlockNum; o++) {
							strcat(order[o], "-l ");
							strcat(order[o], argv[pi]);
							strcat(order[o], " ");
							strcpy(G1[o].orders, order[o]);
						}
						break;

					case 'O':
					case 'o':
						pi++;
						for (o = 0; o < BlockNum; o++) {
							strcat(order[o], "-o ");
							strcat(order[o], argv[pi]);
							strcat(order[o], " ");
							strcpy(G1[o].orders, order[o]);
						}
						if(!atoi(argv[pi])){ printf("close the complementart palindrome mode. \n"); }
						break;
				}
				pi++;
			}
			else break;
		}
		//if fileszie>5.6G
		if(filesize > 6000000000){ 	linecounts = ifnormal(FASTQName,pos);	}
        	//fileszie>5.6G,multithread read file
		if(linecounts !=0){
			printf("File size greater than 6G,wait...");
			multireadfile(linecounts,BlockNum,FASTQName,filesize);
			
			//Create thread read
			pthread_t thread1[BlockNum];
			for (i = 0; i < BlockNum - 1; i++) { pthread_create(&thread1[i], 0, &block_thread, &P1[i]); }
			for (i = 0; i < BlockNum - 1; i++) { pthread_join(thread1[i], 0); }

			fastqw = fopen(OutFile[BlockNum - 1], "wb");
			if ((fastqr = fopen(FASTQName, "rb")) == NULL) { printf("can not find this file:%s\n!",FASTQName); exit(0); }
			fseek(fastqr,P1[BlockNum - 1].readpos,SEEK_SET);
			while (fgets(tempStr,ReadLen,fastqr) != NULL){
				fputs(tempStr,fastqw);
			}
			
			fclose(fastqw);
			fclose(fastqr);

			//Create thread map
			pthread_t thread[BlockNum];
			for (i = 0; i < BlockNum; i++) { pthread_create(&thread[i], 0, &map_thread, &G1[i]); }
			for (i = 0; i < BlockNum; i++) { pthread_join(thread[i], 0); }
		}

		//fileszie<5.6G or not a normal file
		if(linecounts == 0 ){  readfile(linecounts,BlockNum,FASTQName,filesize,OutFile);  }
		
		//Create thread
		
		//delete block file
		for(f = 0; f < BlockNum; f++) { remove(OutFile[f]); }

		//Merge
		//for (i = 0; i < BlockNum; i++) { strcat(OutFile[i], ".map.txt"); }
		//Merge_file(input_name,ref_name,OutFile);
		pthread_t mt1,mt2,mt3,mt4,mt6;
		pthread_create(&mt1, 0, &Merge_pos_file, &M1);
		pthread_create(&mt2, 0, &Merge_cigar_file, &M1);
		pthread_create(&mt3, 0, &Merge_add_file, &M1);
		pthread_create(&mt4, 0, &Merge_cor_file, &M1);
		//pthread_create(&mt5, 0, &Merge_meta_file, &M1);
		pthread_create(&mt6, 0, &Merge_qs_file, &M1);
		pthread_join(mt1, 0);
		pthread_join(mt2, 0);
		pthread_join(mt3, 0);
		pthread_join(mt4, 0);
		//pthread_join(mt5, 0);
		pthread_join(mt6, 0);
		//Merge_pos_file(input_name,ref_name,OutFile);
		//Merge_cigar_file(input_name,ref_name,OutFile);
		//Merge_add_file(input_name,ref_name,OutFile);
		//Merge_cor_file(input_name,ref_name,OutFile);
		
		//delete block map file
		for (i = 0; i < BlockNum; i++) { strcat(OutFile[i], ".map.txt"); }
		for(f = 0; f < BlockNum; f++) { remove(OutFile[f]); }

		//Output
		Output_print(BlockNum);

		//compress
		LWFQZip_compression(input_name,ref_name,assemble_flag,totleread,BlockNum,highestCom_flag);

	}
	if(decompress_flag){
		if(assemble_flag==1){	
			get_fastq_name(ref_name,input_name);
			strcat(ref_name,".fasta");
		}
		LWFQZip_decompression(input_name,ref_name,assemble_flag,highestCom_flag);
	}

	return 0;
}

void print_help()
{
	fprintf(stderr,
			"LW-FQZip 2 -- Reference-based compression of long-read FASTQ files\n"
			"Usage: LWFQZip2 <mode>...[options] ...\n"
			"Mode:\n"
			"  -c  --compression\n"
			"  -d  --decompression\n\n"

			"Compression/Decompression Options:\n"
			"  -i, --input FASTQ file or compressed file.\n"
			"  -r, --input Reference file.\n"
			"  -h, --help  print this message\n"
			"  -v, --version display program version\n"
			"  -m, --maximal read length,range from 30000 to 300000 (Default: 300000)\n"
			"  -s, --Calculate the counts of the prefixes.\n"
			"  -g, --best compression model but slowest\n for example: LWFQZip2 -c -i input -r reference -g.\n		LWFQZip2 -d -i input.lz -r reference -g.\n"
			"  -a, --assemble-based model, An optional amount (Default: 0.3 percent of the original file size) of reads, which contains the predefined prefix (Default: 'CG', could be combined to be an artificial reference. At the end of the package, this artificial reference is included.\n for example: LWFQZip2 -c -i input -a 0.003(Default: '-a 0.003').\n		LWFQZip2 -d -i input.lz -a.\n"
			
			"Mapping Options Options:\n"
			"  -b, --the number of mapping thread(Default: 100, mininum:  6 )\n"
			"  -p, --specify the kmer prefixes, e.g.,'CG', 'AT', and 'TAG' (Default: '-p CG'). 'AA' is not recommended as a prefix.\n"
			"  -k, --length of a kmer used in locate local alignment. (Default: '-k 8')\n"
			"  -e, --the tolerance rate of mismatches.(Default: '-e 0.05')\n"
			"  -L, --the mini length of a legal alignment.(Default: '-l 12')\n"
			"  -o, --open the complementart palindrome mode.(Default: '-o 1' means open the complementart palindrome mode.)\n"
	       );
}

void print_version()
{
	char *VERSION="2.0";
	printf("LWFQZip %s\n", VERSION);
}
