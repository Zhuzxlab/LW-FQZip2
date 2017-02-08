#include "compress.h"
#include <pthread.h>

FILE *fastq=NULL,*meta=NULL,*qs=NULL,*ref_g=NULL,*combine=NULL;
FILE *pos=NULL,*cigar=NULL,*add=NULL,*cor=NULL,*sam=NULL;
char meta_name[50]={'\0'},qs_name[50]={'\0'};
char pos_name[50]={'\0'},cigar_name[50]={'\0'};
char cor_name[50]={'\0'},add_name[50]={'\0'};
char sam_name[50]={'\0'},combine_name[50]={'\0'};
char map_name[50]={'\0'};
char *input_n, *ref_n;
pthread_t t1, t2;

void* m_malloc(size_t n)
{
	void* p = malloc(n);
	if (p == NULL) {
		fprintf(stderr, "Can not allocate %zu bytes.\n", n);
		exit(EXIT_FAILURE);
	}
	return p;
}

FILE *m_fopen(char *file_name,int types)
{
	/* The two types:
0 : read only "rb"
1 : write only "wb"
*/
	FILE *fp= NULL;
	if((types!= 0)&&(types!= 1)){
		fprintf(stderr, "Error types.");
		exit(EXIT_FAILURE);
	}
	if(types == 0)
		fp = fopen(file_name, "rb");
	else fp = fopen(file_name, "wb");
	if (fp == NULL) {
		fprintf(stderr,"Error opening file.%s",file_name);
		exit(EXIT_FAILURE);
	}
	return fp;
}

void init_compression(char *fastq_name,char *ref_name)
{
	strcpy(meta_name,fastq_name);
	m_strcat(meta_name,".meta.txt");
	strcpy(qs_name,fastq_name);
	m_strcat(qs_name,".qs.txt");
	strcpy(pos_name,fastq_name);
	m_strcat(pos_name,".pos.txt");
	strcpy(cigar_name,fastq_name);
	m_strcat(cigar_name,".cigar.txt");
	strcpy(cor_name,fastq_name);
	m_strcat(cor_name,".cor.txt");
	strcpy(add_name,fastq_name);
	m_strcat(add_name,".add.txt");
	strcpy(sam_name,fastq_name);
	m_strcat(sam_name,".map.txt");

	fastq=m_fopen(fastq_name,0);
	meta=m_fopen(meta_name,1);
	qs=m_fopen(qs_name,1);
	ref_g=m_fopen(ref_name,0);
	pos=m_fopen(pos_name,1);
	cor=m_fopen(cor_name,1);
	cigar=m_fopen(cigar_name,1);
	add=m_fopen(add_name,1);
	sam=m_fopen(sam_name,0);
}
/*
   void file_combine(char *fastq_name)
   {
   strcpy(combine_name,fastq_name);
   m_strcat(combine_name,".combine.txt");
   combine=m_fopen(combine_name,1);


   }
   */


size_t get_file_size(char *file)
{
	size_t file_size=-1;
	struct stat statbuff;
	if(stat(file, &statbuff) < 0){
		return file_size;
	}else{
		file_size = statbuff.st_size;
	}
	return file_size;
}

char* get_meta(char *meta_temp,char *metadta)
{
	sscanf(metadta,"%*s%s",meta_temp);
	return meta_temp;
}

char* get_qs(char *qs_temp,char *qsdta)
{
	sscanf(qsdta,"%s",qs_temp);
	return qs_temp;
}

char *get_pos(char *pos_temp,char *map_result)
{
	sscanf(map_result,"%*s%*s%s",pos_temp);
	return pos_temp;
}



char *get_flag(char *flag,char *map_result)
{
	sscanf(map_result,"%*s%s",flag);
	return flag;
}

char *get_seq(char *seq_temp,char *map_result)
{
	sscanf(map_result,"%*s%*s%*s%*s%s",seq_temp);
	return seq_temp;
}

char *get_cigar(char *samcigar,char *map_result)
{
	sscanf(map_result,"%*s%*s%*s%s",samcigar);
	return samcigar;
}

void get_meta_qscore(char *file)
{
	size_t buffer_size=get_file_size(file),buf_count=0;
	if(buffer_size<block_size)
		buf_count = buffer_size;
	else buf_count = block_size;

	char *_buffer=m_malloc(buf_count+read_len);
	char *read_temp=m_malloc(read_len);
	char *meta_temp=m_malloc(read_len);
	char *qs_temp=m_malloc(read_len);
	char *meta_mid=m_malloc(read_len);
	char *meta_part=m_malloc(read_len);

	memset(_buffer,'\0',sizeof(char)*(buf_count+read_len));
	memset(read_temp,'\0',sizeof(char)*read_len);
	memset(meta_temp,'\0',sizeof(char)*read_len);
	memset(qs_temp,'\0',sizeof(char)*read_len);
	memset(meta_mid,'\0',sizeof(char)*read_len);
	memset(meta_part,'\0',sizeof(char)*read_len);

	size_t fread_num=0,read_num=1,i=0,j=0,k=0,l=0,codelen=0;
	char num[20]={'\0'};
	//Copy the first metadata
	fgets(_buffer,read_len,fastq);
	get_meta(meta_mid,_buffer);

	char *ptr=NULL;
	ptr=strstr(meta_mid,"length=");
	for(i=0,l=0;i<strlen(_buffer);i++,l++)
	{
		if(_buffer[i]==' '&&_buffer[i+1]=='l'&&_buffer[i+2]=='e')//remove "length="
			break;
		else *(meta_temp+l)=*(_buffer+i);
	}
	fputs(meta_temp,meta);
	fputc('\n',meta);
	memset(meta_temp,'\0',strlen(meta_temp));
	//memset(meta_mid,'\0',strlen(meta_mid));

	while((fread_num=fread(_buffer,sizeof(char),buf_count,fastq)))
	{
		_buffer[fread_num]='\0';
		for(i=0;i<fread_num;i++)
		{
			*(read_temp+j)=*(_buffer+i);
			j++;
			if(*(_buffer+i)=='\n')
			{
				read_num++;
				if(read_num%4==1){ //Incremental coding of metadata
					get_meta(meta_temp,read_temp);
					if(ptr==NULL){
						for(k=0;k<strlen(meta_temp);k++){
							if(*(meta_temp+k)!=*(meta_mid+k)){
								sprintf(num,"%zu",k);
								for(l=0;k<strlen(meta_temp);k++,l++)
									*(meta_part+l)=*(meta_temp+k);
							}
						}
						if(k==strlen(meta_temp)){
							memset(num,'\0',strlen(num));
							sprintf(num,"%zu",k);
						}
						memset(meta_mid,'\0',strlen(meta_mid));
						strcpy(meta_mid,meta_temp);
						fputs(num,meta);
						if(k!=strlen(meta_temp)){
							fputc(' ',meta);
							fputs(meta_part,meta);
						}
						fputc('\n',meta);
						memset(num,'\0',strlen(num));
						memset(meta_part,'\0',strlen(meta_part));
					}
					else{
						fputs(meta_temp,meta);
						fputc('\n',meta);
					}
				}
				if(read_num%4==0){
					//codelen is the coded length of qualityScore
					//codelen=rll_quality_score(qs_temp, read_temp);
					//codelen=fwrite(qs_temp,1,codelen,qs);
					get_qs(qs_temp,read_temp);
					fwrite(qs_temp,1,strlen(qs_temp),qs);
					fputc('\n',qs);
				}
				memset(read_temp,'\0',strlen(read_temp));
				j=0;
			}
		}
	}
	//fprintf(stderr,"read_num=%d\n",read_num/4);
	free(_buffer);  _buffer=NULL;
	free(meta_temp); meta_temp=NULL;
	free(meta_mid);  meta_mid=NULL;
	free(meta_part);  meta_part=NULL;
	free(qs_temp);  qs_temp=NULL;
	free(read_temp);  read_temp=NULL;
}

void quit_compression_dispose()
{
	fclose(fastq);
	fclose(meta);
	fclose(qs);
	fclose(ref_g);
	fclose(pos);
	fclose(cigar);
	fclose(add);
	fclose(cor);
	fclose(sam);
}

char m_toupper(char ch)
{
	if (ch>96)	return ch-32;
	else return ch;
}

void split_map_result(char *ref_name)
{
	size_t fread_num=0,i=0,j=0,buffer_size=0,buf_count=0;
	//Load reference genome into memory.
	buffer_size=get_file_size(ref_name);
	buf_count = buffer_size+10000;

	char *ref_temp=m_malloc(buf_count);
	char *ref_gen=m_malloc(buf_count),ref_ch;
	memset(ref_temp,'\0',sizeof(char)*buf_count);
	memset(ref_gen,'\0',sizeof(char)*buf_count);
	//Ignore the title line
	fgets(ref_temp,2000,ref_g);
	j=0;
	while((fread_num=fread(ref_temp,sizeof(char),buf_count, ref_g))){
		for(i=0;i<fread_num;i++){
			if( ( ref_ch=m_toupper(*(ref_temp+i)) )>64){
				*(ref_gen+j)=ref_ch;
				j++;
			}
		}
	}
	free(ref_temp);ref_temp=NULL;
	//Load reference end.

	buffer_size=get_file_size(sam_name),buf_count=0;
	if(buffer_size<block_size)
		buf_count = buffer_size;
	else buf_count = block_size;

	char *map_result=m_malloc(read_len*2);
	char *seq=m_malloc(read_len);
	char *base_temp=m_malloc(read_len);
	char *match=m_malloc(read_len);
	char *add_temp=m_malloc(read_len);
	char *clip_seq=m_malloc(read_len);
	char *base_=m_malloc(read_len);

	memset(base_,'\0',read_len*sizeof(char));
	memset(clip_seq,'\0',sizeof(char)*read_len);
	memset(map_result,'\0',sizeof(char)*read_len*2);
	memset(seq,'\0',sizeof(char)*read_len);
	memset(base_temp,'\0',sizeof(char)*read_len);
	memset(match,'\0',sizeof(char)*read_len);
	memset(add_temp,'\0',sizeof(char)*read_len);

	char flag[10]={'\0'},pos_temp[30]={'\0'},samcigar[2000]={'\0'};

	/*Ignore the two lines in 'read.fastq.map.txt' file */
	fgets(map_result,(read_len*2),sam);
	fgets(map_result,(read_len*2),sam);
	char *buffer=m_malloc(buf_count+read_len);
	memset(buffer,'\0',sizeof(char)*(buf_count+read_len));

	j=0;
	while((fread_num=fread(buffer,sizeof(char),buf_count, sam)))
	{
		buffer[fread_num]='\0';
		for(i=0;i<fread_num;i++)
		{
			*(map_result+j)=*(buffer+i);
			j++;
			if(*(buffer+i)=='\n')
			{
				get_pos(pos_temp,map_result);
				get_seq(seq,map_result);
				get_cigar(samcigar,map_result);
				get_flag(flag,map_result);
				new_split_cigar(add_temp,seq,samcigar,flag,clip_seq);
				new_get_base(base_temp,add_temp,samcigar,ref_gen,pos_temp,base_);
				get_mismatch(match,seq,base_temp,pos_temp);

				int flag_num=atoi(flag);
				//                char exact_map[30]={'\0'};
				//                int seq_len=strlen(seq);
				//                sprintf(exact_map,"%d",seq_len);
				//                exact_map[strlen(exact_map)]='M';
				if(flag_num==16)//Palindrome
					fputc('0',cigar);
				fputs(samcigar,cigar);
				fputc('\n',cigar);

				fputs(add_temp,add);
				fputc('\n',add);
				fputs(pos_temp,pos);
				fputc('\n',pos);
				fputs(match,cor);

				memset(add_temp,'\0',strlen(add_temp));
				memset(base_temp,'\0',strlen(base_temp));
				memset(match,'\0',strlen(match));
				memset(map_result,'\0',strlen(map_result));
				memset(clip_seq,'\0',strlen(clip_seq));
				memset(base_,'\0',strlen(base_));
				j=0;
			}
		}
	}
	quit_compression_dispose();
	free(ref_gen); ref_gen=NULL;
	free(map_result); map_result=NULL;
	free(seq); seq=NULL;
	free(base_temp); base_temp=NULL;
	free(match); match=NULL;
	free(buffer); buffer=NULL;
	free(add_temp); add_temp=NULL;
	free(clip_seq); clip_seq=NULL;
	free(base_);  base_=NULL;
}

char *get_mismatch(char *match,const char *seq,const char *base,char *pos)
{
	if(pos[0]==0)
	{
		strcpy(match,"\n");
		return match;
	}
	size_t i=0,temp=0,first=0;
	char num[30]={'\0'};
	for(i=0;i<strlen(seq);i++)
	{
		if(*(seq+i)!=*(base+i))
		{
			temp=i-first;
			sprintf(num,"%zu",temp);
			first=i;
			m_strcat(match,num);
			match[strlen(match)]=*(seq+i);
			memset(num,'\0',strlen(num));
		}
	}
	match[strlen(match)]='\n';
	return match;
}

char *new_get_base(char *base,char *add_,char *samcigar,char *ref_gen,char *pos,char *base_)
{
	if(add_[strlen(add_)-1]=='\n')
		add_[strlen(add_)-1]='\0';
	char map_info='\0',cigar_num[30]={'\0'};
	size_t map_pos=0,i=0,j=0,map_len=0,k=0,a=0;
	bool flag=false;

	map_pos=atol(pos);
	if(map_pos==0)
	{
		strcpy(base,add_);
		return base;
	}
	else
	{
		for(i=0;i<strlen(samcigar);i++)
		{
			if(*(samcigar+i)=='\n')
				break;
			if(*(samcigar+i)>='0'&&(*(samcigar+i))<='9')
			{
				*(cigar_num+j)=*(samcigar+i);
				j++;
				if(*(samcigar+i+1)>='A'&&(*(samcigar+i+1))<='Z')
					map_info=*(samcigar+i+1);
			}
			else
			{
				map_len=atol(cigar_num);
				memset(cigar_num,'\0',30);
				if(map_info=='M')
					flag=true;
				switch(map_info)
				{
					case 'M':{
							 for(j=map_pos;j<map_pos+map_len;j++){
								 *(base_+k)=ref_gen[j-1];
								 k++;
							 }
							 m_strcat(base,base_);
							 memset(base_,'\0',strlen(base_));
							 k=0;
							 j=0;
							 map_pos=map_pos+map_len;
							 break;
						 }
					case 'S':{
							 for(j=a;j<a+map_len;j++){
								 *(base_+k)=*(add_+j);
								 k++;
							 }
							 m_strcat(base,base_);
							 memset(base_,'\0',strlen(base_));
							 k=0;
							 j=0;
							 a=a+map_len;
							 if(flag)
								 map_pos=map_pos+map_len;
							 break;
						 }
					case 'I':{
							 for(j=a;j<a+map_len;j++){
								 *(base_+k)=*(add_+j);
								 k++;
							 }
							 m_strcat(base,base_);
							 memset(base_,'\0',strlen(base_));
							 k=0;
							 j=0;
							 a=a+map_len;
							 break;
						 }
					case 'D':{
							 j=0;
							 if(flag)
								 map_pos=map_pos+map_len;
							 break;
						 }
					case 'N':{
							 j=0;
							 if(flag)
								 map_pos=map_pos+map_len;
							 break;
						 }
					case 'H':break;
					case 'P':break;
					default:break;
				}
			}
		}
	}
	return base;
}


char *new_split_cigar(char *add_temp,char *seq,char *samcigar,char *flag,char *clip_seq)
{
	size_t i=0,j=0,k=0;
	char NUM[30]={'\0'},unmap[]={"*"},char_cig[100]={'\0'};
	int num_cig[100],num_cig_sum[100];

	if(strcmp(samcigar,unmap)==0){
		strcpy(add_temp,seq);
		return add_temp;
	}
	else
	{
		for(i=0;i<strlen(samcigar);i++){
			num_cig[i]=0;
			num_cig_sum[i]=0;
		}
		for(i=0;i<strlen(samcigar);i++)
		{
			if(*(samcigar+i)>='0'&&(*(samcigar+i))<='9'){
				NUM[j++]=*(samcigar+i);
			}else{
				char_cig[k]=*(samcigar+i);
				num_cig[k]=atoi(NUM);
				if(k>=1){
					num_cig_sum[k]=num_cig[k-1]+num_cig_sum[k-1];
					if(char_cig[k-1]=='D')
						num_cig_sum[k]=num_cig_sum[k]-num_cig[k-1];
				}
				else
					num_cig_sum[k]=0;
				k++;
				memset(NUM,'\0',strlen(NUM));
				j=0;
			}
		}
		for(i=0;i<strlen(char_cig);i++)
		{
			k=0;
			switch(char_cig[i])
			{
				case 'S':{
						 for(j=num_cig_sum[i];j<(num_cig[i]+num_cig_sum[i]);j++){
							 *(clip_seq+k)=*(seq+j);
							 k++;
						 }
						 m_strcat(add_temp,clip_seq);
						 memset(clip_seq,'\0',strlen(clip_seq));
						 break;
					 }
				case 'I':{
						 for(j=num_cig_sum[i];j<(num_cig[i]+num_cig_sum[i]);j++){
							 *(clip_seq+k)=*(seq+j);
							 k++;
						 }
						 m_strcat(add_temp,clip_seq);
						 memset(clip_seq,'\0',strlen(clip_seq));
						 break;
					 }
				case 'M': break;
				case 'D': break;
				case 'N': break;
				case 'H': break;/*Hard clipping.*/
				case 'P': break;
				default:  break;
			}
		}
		return add_temp;
	}
}

char *m_strcat(char *dst, const char* src)
{
	while(*dst)   dst++;
	while((*dst ++ = *src++));
	return --dst ;
}

void  FQZip_compression(char *fastq_name, unsigned long long totalread,int BlockNum, char *ref_name, int assemble_flag,int highestCom_flag)
{


	char meta_lz[100]={'\0'},qs_lz[100]={'\0'};
	char pos_lz[100]={'\0'},cigar_lz[100]={'\0'};
	char cor_lz[100]={'\0'},add_lz[100]={'\0'};
	char fasta_lz[100]={'\0'};
	char order[1000]={'\0'};
	strcpy(meta_lz,meta_name);
	m_strcat(meta_lz,".lz");
	strcpy(qs_lz,qs_name);
	m_strcat(qs_lz,".lz");
	strcpy(pos_lz,pos_name);
	m_strcat(pos_lz,".lz");
	strcpy(cigar_lz,cigar_name);
	m_strcat(cigar_lz,".lz");
	strcpy(cor_lz,cor_name);
	m_strcat(cor_lz,".lz");
	strcpy(add_lz,add_name);
	m_strcat(add_lz,".lz");
	if(assemble_flag==0)	//reference-based model
	sprintf(order,"%s %s %s %s %s %s %s %s %s %s %s %s %s %lld %d %d","./FQZip c  ",pos_name,cigar_name,cor_name,meta_name,add_name,qs_name,pos_lz,cigar_lz,cor_lz,meta_lz,add_lz,qs_lz, totalread,BlockNum,highestCom_flag);
	else		//assemble-based model
	{sprintf(order,"%s %s %s %s %s %s %s %s %s %s %s %s %s %lld %d %s %d","./FQZip c  ",pos_name,cigar_name,cor_name,meta_name,add_name,qs_name,pos_lz,cigar_lz,cor_lz,meta_lz,add_lz,qs_lz, totalread,BlockNum,ref_name,highestCom_flag);	
	strcpy(fasta_lz,ref_name);
	m_strcat(fasta_lz,".lz");
	}
	if (-1 == (system(order)))printf("Please make sure you have installed the FQZip correctly");
	char lwfqzip_name[50]={'\0'};
	char tar_name[500]={'\0'};
	char tar1_name[500]={'\0'};
	m_strcat(tar_name,meta_lz);
	m_strcat(tar_name," ");
	//m_strcat(tar_name,qs_lz);
	m_strcat(tar_name,qs_name);
	m_strcat(tar_name,"*.lz");
	m_strcat(tar_name," ");
	m_strcat(tar_name,pos_lz);
	m_strcat(tar_name," ");
	m_strcat(tar_name,cigar_lz);
	m_strcat(tar_name," ");
	m_strcat(tar_name,cor_lz);
	m_strcat(tar_name," ");
	m_strcat(tar_name,add_lz);
	m_strcat(tar_name," ");
	char dir[500] = { '\0' };
	char file_name[256]={'\0'};
	char *p;
	int i;
	bool t = 0;
	strcpy(dir, fastq_name);
	for (i = 0; i < 256; i++)
	{
		if (dir[i] == '/')
		{
			t = 1;
			break;
		}
	}
	if (t == 1)   //not LWFQZip2 dir
	{
		p = rindex(dir, '/');p++;
		for (i = 0; i < 100; i++)
		{file_name[i]=*p;   *p = '\0';p++;}
		char local_dir[256];
		m_strcat(tar1_name,file_name);
		m_strcat(tar1_name,".meta.txt.lz ");
		m_strcat(tar1_name,file_name);
		m_strcat(tar1_name,".qs.txt*.lz ");
		m_strcat(tar1_name,file_name);
		m_strcat(tar1_name,".cigar.txt.lz ");
		m_strcat(tar1_name,file_name);
		m_strcat(tar1_name,".cor.txt.lz ");
		m_strcat(tar1_name,file_name);
		m_strcat(tar1_name,".pos.txt.lz ");
		m_strcat(tar1_name,file_name);
		m_strcat(tar1_name,".add.txt.lz ");
		if(assemble_flag==1)	//reference-based model
		{
		m_strcat(tar1_name,file_name);
		m_strcat(tar1_name,".fasta.lz ");
		}
		getcwd(local_dir,sizeof(local_dir));
		chdir(dir);
		strcpy(lwfqzip_name,file_name);
		m_strcat(lwfqzip_name,".lz");
		sprintf(order,"%s %s %s ","tar -c -f ",lwfqzip_name,tar1_name);
		if (-1 == (system(order)))printf("Please make sure you have installed the tar tool correctly");
		chdir(local_dir);

	}
	if (t==0)       //LWFQZip2 dir
	{
		if(assemble_flag==1)	//reference-based model
		{
		m_strcat(tar_name,fasta_lz);
		}
		strcpy(lwfqzip_name,fastq_name);
		m_strcat(lwfqzip_name,".lz");
		sprintf(order,"%s %s %s ","tar -c -f ",lwfqzip_name,tar_name);
		if (-1 == (system(order)))printf("Please make sure you have installed the tar tool correctly");
	}
	sprintf(order,"rm %s*",qs_name);
	if (-1 == (system(order)))printf("Can't delete the qs blcoks\n");
	remove(meta_lz);
	remove(qs_lz);
	remove(cigar_lz);
	remove(cor_lz);
	remove(pos_lz);
	remove(add_lz);
	remove(meta_name);
	remove(qs_name);
	remove(cigar_name);
	remove(cor_name);
	remove(pos_name);
	remove(add_name);
	if(assemble_flag==1)
	{
	remove(fasta_lz);
	remove(ref_name);
	}

}


void *thread_1(void *args)
{
	get_meta_qscore(input_n);
}
void *thread_2(void *args)
{
	split_map_result(ref_n);
}

void LWFQZip_compression(char *input_name,char *ref_name,int assemble_flag, unsigned long long totalread,int BlockNum,int highestCom_flag)
{
	char *ptr=NULL;
	ptr=strstr(input_name,".fastq");
	strcpy(map_name,input_name);
	input_n = input_name;
	ref_n = ref_name;
	m_strcat(map_name,".map.txt");
	if(!strstr(input_name,".fastq")&&!strstr(input_name,".fq")){
		fprintf(stderr,"Please input correct FASTQ file\n");
		exit(EXIT_FAILURE);
	}
	init_compression(input_name,ref_name);
	get_meta_qscore(input_name);
	split_map_result(ref_name);
	FQZip_compression(input_name, totalread,BlockNum,ref_name,assemble_flag,highestCom_flag);
	remove(map_name);
}


/*
   align the qulity score
input: read.fq--there should be a qulity score
output: RL code qulity score ---- 8 records or 16 records
// output:  align result or statistic result,such as starray[][]-----should include the file "statisticChar.c"
score range: 33--126
*/
int rll_quality_score(char *buf, char *str) //33 num rule---n is a code parameter
{
	char prech, ch, pprech=0;
	int  i, count=1, len=0, flag=0;  //len means the length of buf

	prech=*str;
	//len+=sprintf(buf+len,"%c",prech); //start character of the line
	*(buf+len++)=prech;

	for(i=1; prech>32; i++ )  //read Quality line
	{
		ch=*(str+i); //no tail
		if(prech==ch)
			count++;
		else
		{

			//if(count>=128) printf("count: %d %c  %d\n%s\n",count,prech,i,str+i);
			if(flag)  //may be double number
			{
				if(count<32)
					//len+=sprintf(buf+len,"%c",count); //C5
					*(buf+len++)=count;
				else if(count<128)
				{
					*(buf+len++)=256-count;
					//len+=sprintf(buf+len,"%c",-count); //-count
				}
				else
				{
					flag=0;
					*(buf+len++)=prech;
					//printf("flag count:%d %c\n",count, prech); //test
					//*(buf+len++)=pprech;
					i--;
					continue;
				}
				if(ch!=pprech) { *(buf+len++)=ch; flag=0; }
			}
			else
			{
				if( count<N_LINK ) //N_LINK == 2
				{
					while( --count >0 ) //change to : *(buf+len++)=prech; when N_LINK==2
						*(buf+len++)=prech;
					*(buf+len++)=ch;
					flag=0;
				}
				else //count > N : RLcode
				{
					if (count<32) //compare to n 33 128 ----char 128<0 255=-1 EOF...
					{
						//len+=sprintf(buf+len,"%c",count); //C5
						*(buf+len++)=count;
					}
					else if( count<128 )
						//len+=sprintf(buf+len,"%c",-count); //C-155
						*(buf+len++)=256-count;

					else
					{
						pprech=0; //no link Code
						//printf("flag=0 count:%d\n",count); //test
						do
						{
							len+=sprintf(buf+len,"%c%c",128,prech); //char 255=-1 so use -n instead
							count-=128;
						}while(count>127);
						if(count>31)
							//len+=sprintf(buf+len,"%c",-count);
							*(buf+len++)=256-count;
						else if(count>=N_LINK)
							//len+=sprintf(buf+len,"%c",count); //C5
							*(buf+len++)=count;
						else
						{
							if(count>0)
								while( --count >0 )
									*(buf+len++)=prech;
							else len--; //count==0, count%128==0 prech is too much
						}
					}
					if(ch==pprech)
						flag=1;
					else
					{
						flag=0;
						*(buf+len++)=ch;
					}

				}
			}
			pprech=prech; //record last code char--if count>128 maybe a bug
			prech=ch;
			count=1;
		}
	}

	if(ch<33) //endline
		len--;
	*(buf+len++)=32; //code end symbol
	return len;
}
