#include "lwfqzip.h"
#include <pthread.h>

#define ReadLen 300000
#define MaxBlockNum 50
#define _FILE_OFFSET_BITS 64

int BlockNum = 10; //Default file block && number of thread

struct G
{
	char input[50];
	char orders[100];
}G1[MaxBlockNum];

struct P
{
	char filename[50];
	char file[50];
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
	//printf("subfastqname:%s\n",P2->filename);
	//printf("fastqname:%s\n",P2->file);

	P2->subfastq = fopen(P2->filename, "wb");
	if ((P2->fastq = fopen(P2->file, "rb")) == NULL) { printf("can not find this file:%s\n", P2->file); exit(0); }

	//printf("subblockline:%ld\n",P2->readle);
	//printf("blockpos:%ld\n",P2->readpos);
	fseek(P2->fastq,P2->readpos,SEEK_SET);
	while (fgets(tempStr1, ReadLen,P2->fastq) != NULL){
		if (subline  == ((P2->readle)) * 4)
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
	char prekeyorder[100] = {0};
	char prekeyorderAT[100] = {0},prekeyorderCG[100] = {0},prekeyorderGT[100] = {0},prekeyorderAC[100] = {0};
	char prekeyorderTAG[100] = {0},prekeyorderGCA[100] = {0},prekeyorderTAT[100] = {0},prekeyorderACT[100] = {0};
	//char prekeyname[5]  =  {'\0'};
	char ref_name[50] = {'\0'};

	//strcpy(prekeyname, prename);
	strcpy(ref_name, refname);

	sprintf(prekeyorderAT,"%s%s%s%s%s%s","grep -o ","'","AT","' ",ref_name," | wc -l ");
	sprintf(prekeyorderCG,"%s%s%s%s%s%s","grep -o ","'","CG","' ",ref_name," | wc -l ");
	sprintf(prekeyorderGT,"%s%s%s%s%s%s","grep -o ","'","GT","' ",ref_name," | wc -l ");
	sprintf(prekeyorderAC,"%s%s%s%s%s%s","grep -o ","'","AC","' ",ref_name," | wc -l ");
	sprintf(prekeyorderTAG,"%s%s%s%s%s%s","grep -o ","'","TAG","' ",ref_name," | wc -l ");
	sprintf(prekeyorderGCA,"%s%s%s%s%s%s","grep -o ","'","GCA","' ",ref_name," | wc -l ");
	sprintf(prekeyorderTAT,"%s%s%s%s%s%s","grep -o ","'","TAT","' ",ref_name," | wc -l ");
	sprintf(prekeyorderACT,"%s%s%s%s%s%s","grep -o ","'","ACT","' ",ref_name," | wc -l ");


	//printf("prekeyorder is:%s\n",prekeyorder);
	//printf("prekeyorderAT is:%s\nprekeyorderCG is:%s\nprekeyorderGT is:%s\nprekeyorderAC is:%s\n",prekeyorderAT,prekeyorderCG,prekeyorderGT,prekeyorderAC);

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

int main(int argc, char *argv[])
{
	bool compress_flag=false;
	bool decompress_flag=false;
	bool statistics_flag=false;
	

	unsigned long long totleread=0,mappedread=0,exactmatch=0,doublemap=0;
	unsigned long long totalbase=0, mappedbase=0;

	int assemble_flag=0;
	int highestCom_flag=0;
	float FASTA_rate=0.003;//assemble-based
	char input_name[50]={'\0'};
	char ref_name[50]={'\0'};
	char order[MaxBlockNum][100] = {0};

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
				printf("FASTA extraction rate: %.2f%%\n",atof(argv[optind])*100);
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
				read_len = strtoul(optarg, NULL, 10);
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

		FILE *fastq = NULL;
		FILE *fastq2 = NULL;
		FILE *MergeMapfile = NULL;
		FILE *Mergefastq0 = NULL;
		FILE *readout=NULL;
		FILE *subfastq[BlockNum];
		FILE *Mergefastq[BlockNum];

		unsigned long long digit[BlockNum][6];
		int i;

		char OutFile[MaxBlockNum][100] = { 0 };
		char tempStr[ReadLen] = { '\0' }; //temporary row number
		char tempMergerow[ReadLen] = { '\0' };
		char FASTQName[50] = { '\0' };
		char prekeyname[5] = { '\0' };
		char MergeMapName[50] = { '\0' };
		long long linecounts = 0; //line counts of file
		long readcounts, subreadcounts, lastblockcounts;
		long c, subline;

		//Sub-block
		strcpy(FASTQName, argv[3]);  //file name

		FILE *fastqr = NULL;
		FILE *fastqw = NULL;
		char readtmp[50] = { '\0' };
		char readline[50] = { '\0' };
		char a;
		long int pos = 1;
		int countline =0;

		struct stat buf;
		stat (FASTQName, &buf);
		long int filesize;
		filesize = buf.st_size;
		printf ("File size (%s): %lld MB.\n",FASTQName, filesize/1024/1024);
		
		struct stat buf1;
		stat (ref_name, &buf1);
		long int gensize;
		gensize= buf1.st_size;
		if(!assemble_flag)
		printf ("File size (%s): %ld MB.\n",ref_name, gensize/1024/1024);
		
		strcpy(prekeyname,"CG");

		//long int prekeycounts;
		//fscanf(readout, "%ld",&digit[i][0]);

		if(assemble_flag==1)
		{
			FILE *OriFASTQ=NULL,*FASTA=NULL;
			char OriFASTQName[500]={'\0'};
			char tempStr[300000]={'\0'};
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
			float FASTA_Maxsize=filesize*FASTA_rate;
			float FASTA_size=0;
			long countline=0;
			int i=0;
			int flag_CG=0;
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
					if(tempStr[i]=='C'&&tempStr[i+1]=='G')
					{flag_CG=1;break;}
					if(flag_CG==1)
					{
						for(i=0;i<strlen(tempStr);i++)
						{if(tempStr[i]=='\n')
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
		
		if(statistics_flag){
			prekeycontent(ref_name);
			exit(0);
		}

		for (i = 0; i < BlockNum; i++){
			char s[10];
			sprintf(s, "%d", i);
			strcat(OutFile[i], FASTQName);
			strcpy(P1[i].file, FASTQName);
			strcat(OutFile[i], s);
			strcat(OutFile[i], ".fastq");
			strcpy(G1[i].input, OutFile[i]);
			strcpy(P1[i].filename, OutFile[i]);
			strcat(order[i], "./LWMapping ");
			if(assemble_flag==1)	//assemble model
			{
			strcat(order[i], FASTQName);
			strcat(order[i], ".fasta");
			}
			else					//reference-base model
				strcat(order[i], argv[5]);
			strcat(order[i], " ");
			strcat(order[i], G1[i].input);
			strcat(order[i], " ");
			strcpy(G1[i].orders, order[i]);
			//printf("Order%s\n",G1[i].orders);
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

		if(filesize > 6000000000){
			if ((fastqr = fopen(FASTQName, "rb")) == NULL) {  printf("can not find this file:%s\n!", FASTQName); exit(0); }
			fseek(fastqr,-pos,SEEK_END);
			while(pos){
				if(a= fgetc(fastqr)=='\n'){
					countline++;
				}
				pos++;
				fseek(fastqr,-pos,SEEK_END);
				if(countline==5) break;
			}
			pos-=2;
			fseek(fastqr,-pos,SEEK_END);
			fgets(tempStr, ReadLen, fastqr);
			sscanf(tempStr,"%[^ ]",readtmp);
			sscanf(readtmp,"%*[^.].%s",readline);
			linecounts = atoi(readline);
			printf("linecounts:%d\n",linecounts);
			fclose(fastqr);
		}
        	//linecounts = 0;
		if(linecounts !=0){
			readcounts = linecounts;
			//printf("readcounts:%ld\n",readcounts);
			subreadcounts = readcounts / BlockNum;
			//printf("subreadcounts:%ld\n",subreadcounts);
			if(readcounts%BlockNum==0){ lastblockcounts = subreadcounts; }
			else{ lastblockcounts = subreadcounts+(readcounts - subreadcounts*BlockNum); }
			//printf("lastblockcounts:%ld\n",lastblockcounts);

			for (i = 0; i < BlockNum - 1; i++){
				P1[i].readle = subreadcounts;
			//	printf("The %d has %d reads.\n",i,P1[i].readle);
			}
			if ((fastqr = fopen(FASTQName, "rb")) == NULL) { printf("can not find this file:%s\n!", FASTQName); exit(0); }


			long int filepos;
			long int adjustpos;
			long int adjustline;
			char adjusttempStr[ReadLen] = { '\0' };
			char adjusttempStr1[ReadLen] = { '\0' };
			char adjusttempStr2[ReadLen] = { '\0' };
			char u,ua;
			P1[0].readpos = 0;
			int poscount;
			for(poscount=1;poscount<BlockNum;poscount++){
				filepos=(filesize/BlockNum) * poscount;
				adjustpos=filepos;
				while(adjustpos){
					adjustpos++;
					fseek(fastqr,adjustpos,SEEK_SET);
					if(a= fgetc(fastqr)=='\n'){
						int posline=0;
						fgets(adjusttempStr, ReadLen, fastqr);
						ua = adjusttempStr[0];
						for(i=0;i<strlen(adjusttempStr);i++){
							u = adjusttempStr[i];
							if (isblank(u) && ua == '@'){
								posline=1;
							}
						}
						if(posline){
							long int subline;
							long int tmpline;
							sscanf(adjusttempStr,"%[^ ]",adjusttempStr1);
							sscanf(adjusttempStr1,"%*[^.].%s",adjusttempStr2);
							adjustline= atoi(adjusttempStr2);
							//printf("Approximate position:%ld\n",adjustline);
							while(adjustline != (subreadcounts * poscount+1)){
								if(adjustline > (subreadcounts * poscount+1)){
									tmpline = adjustline - (subreadcounts * poscount+1);
									subline = tmpline*4;
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
									subline = tmpline*4;
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
							//printf("Ture position:  %ld\n",adjustline);
							//printf("Ture Pointer position:  %ld\n",adjustpos);
							break;
						}
					}
				}
			}
			fclose(fastqr);

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
			//Create thread
			pthread_t thread[BlockNum];
			for (i = 0; i < BlockNum; i++) { pthread_create(&thread[i], 0, &map_thread, &G1[i]); }
			for (i = 0; i < BlockNum; i++) { pthread_join(thread[i], 0); }
		}
		if(linecounts == 0 )
		{
			pthread_t thread[BlockNum];
			if ((fastq = fopen(FASTQName, "rb")) == NULL) {  printf("can not find this file:%s\n!", FASTQName); exit(0);}
			while ((c = fgetc(fastq)) != EOF) { if (c == '\n') linecounts++; }
			fclose(fastq);
			readcounts = linecounts / 4;
			subreadcounts = readcounts / BlockNum;
			if(readcounts%BlockNum==0){ lastblockcounts = subreadcounts; }
			else{ lastblockcounts = subreadcounts+(readcounts - subreadcounts*BlockNum); }
			if ((fastq2 = fopen(FASTQName, "rb")) == NULL) { printf("can not find this file:%s\n!", FASTQName); exit(0); }
			for (i = 0; i < BlockNum - 1; i++){
				subline = 0;
				subfastq[i] = fopen(OutFile[i], "wb");
				while (fgets(tempStr, ReadLen, fastq2) != NULL){
					subline++;
					fputs(tempStr, subfastq[i]);
					if (subline  == subreadcounts * 4)
						break;
				}
				fclose(subfastq[i]);
				pthread_create(&thread[i], 0, &map_thread, &G1[i]);
			}
			subline = 0;
			subfastq[BlockNum - 1] = fopen(OutFile[BlockNum - 1], "wb");
			while (fgets(tempStr, ReadLen, fastq2) != NULL){
				if (subline / 4 == lastblockcounts)
					break;
				subline++;
				fputs(tempStr, subfastq[BlockNum - 1]);
			}
			fclose(subfastq[BlockNum - 1]);
			pthread_create(&thread[BlockNum - 1], 0, &map_thread, &G1[BlockNum - 1]);
			fclose(fastq2);
			for (i = 0; i < BlockNum; i++) { pthread_join(thread[i], 0); }

		}
		//Create thread
		//pthread_t thread[BlockNum];
		//for (i = 0; i < BlockNum; i++) { pthread_create(&thread[i], 0, &map_thread, &G1[i]); }
		//for (i = 0; i < BlockNum; i++) { pthread_join(thread[i], 0); }			 				
		
		//delete block file
		int f;
		for(f = 0; f < BlockNum; f++) { remove(OutFile[f]); }
		//Merge
		strcpy(MergeMapName, input_name);
		strcat(MergeMapName, ".map.txt");
		int k;
		for (k = 0; k < BlockNum; k++) { strcat(OutFile[k], ".map.txt"); }

		MergeMapfile = fopen(MergeMapName, "wb");
		Mergefastq0 = fopen(OutFile[0], "rb");
		fgets(tempMergerow, ReadLen, Mergefastq0);
		fclose(Mergefastq0);
		fputs(tempMergerow, MergeMapfile);
		fprintf(MergeMapfile, "@PG\tID:map\tPN:map\tVN:0.1\tCL:%s %s %s\n", "./LWMapping ", ref_name, input_name);

		for (k = 0; k < BlockNum; k++) {
			Mergefastq[k] = fopen(OutFile[k], "rb");
			fgets(tempMergerow, ReadLen, Mergefastq[k]); //first line
			fgets(tempMergerow, ReadLen, Mergefastq[k]); //second line
			while (fgets(tempMergerow, ReadLen, Mergefastq[k])){
				fputs(tempMergerow, MergeMapfile);
			}
			fclose(Mergefastq[k]);
		}
		fclose(MergeMapfile);

		//delete block map file
		for(f = 0; f < BlockNum; f++) { remove(OutFile[f]); }

		//Output
		if ((readout = fopen("Output.txt","rb")) == NULL) {  //whether the file exists
			printf("can not find this file:%s\n!", "Output.txt");
			exit(0);
		}
		i=0;
		while(i<BlockNum){
			fscanf(readout, "%lld,%lld,%lld,%lld,%lld,%lld",&digit[i][0],&digit[i][1],&digit[i][2],&digit[i][3],&digit[i][4],&digit[i][5]);
			i++;
		}
		int j;
		for(j=0;j<BlockNum;j++){
			totleread+=digit[j][0]; mappedread+=digit[j][1]; exactmatch+=digit[j][2];doublemap+=digit[j][3]; totalbase+=digit[j][4]; mappedbase+=digit[j][5];
			//printf("digit=%lld\n",digit[j][4]);
		}

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

		//compress
		LWFQZip_compression(input_name,ref_name,assemble_flag,totleread,BlockNum,highestCom_flag);

	}
	if(decompress_flag){
		if(assemble_flag==1)
		{	get_fastq_name(ref_name,input_name);
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
			"  -b, --the number of mapping thread(Default: 10, mininum:  6 )\n"
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
