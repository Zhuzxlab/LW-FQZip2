#include "decompress.h"
#include <pthread.h>

FILE *fastq_=NULL,*meta_=NULL,*qs_=NULL,*ref_g_=NULL;
FILE *pos_=NULL,*cigar_=NULL,*add_=NULL,*cor_=NULL;
FILE *base_=NULL,*meta_data=NULL,*qscore=NULL;
char de_fastq_name[50]={'\0'},qscore_name[50]={'\0'};
char meta_name_[50]={'\0'},qs_name_[50]={'\0'};
char pos_name_[50]={'\0'},cigar_name_[50]={'\0'};
char cor_name_[50]={'\0'},add_name_[50]={'\0'};
char base_name_[50]={'\0'},metadata[50]={'\0'};
char fastq_name[50] = { '\0' },fasta_nameC[50]={'\0'};
char *ref_;
pthread_t t1, t2, t3;

char *get_fastq_name(char *fastqname,char *str)
{
	int i=0,j=0;
	for(i=0;i<60;i++)
	{
		if(str[i]=='.'&&str[i+1]=='l'&&(str[i+2]=='w'||str[i+2]=='z'))
			break;
		else fastqname[j++]=str[i];
	}
	return fastqname;
}

void init_decompression(char *fastq_name,char *ref_name_)
{
	strcpy(meta_name_,fastq_name);
	m_strcat(meta_name_,".meta.txt");
	strcpy(qs_name_,fastq_name);
	m_strcat(qs_name_,".qs.txt");
	strcpy(pos_name_,fastq_name);
	m_strcat(pos_name_,".pos.txt");
	strcpy(cigar_name_,fastq_name);
	m_strcat(cigar_name_,".cigar.txt");
	strcpy(cor_name_,fastq_name);
	m_strcat(cor_name_,".cor.txt");
	strcpy(add_name_,fastq_name);
	m_strcat(add_name_,".add.txt");

	strcpy(base_name_,fastq_name);
	m_strcat(base_name_,".base");
	strcpy(metadata,fastq_name);
	m_strcat(metadata,".meta");
	strcpy(qscore_name,fastq_name);
	//
	m_strcat(qscore_name,".qs.txt");

	base_=m_fopen(base_name_,1);
	meta_data=m_fopen(metadata,1);
	//qscore=m_fopen(qscore_name,1);

	ref_g_=m_fopen(ref_name_,0);
	pos_=m_fopen(pos_name_,0);
	meta_=m_fopen(meta_name_,0);
	//qs_=m_fopen(qs_name_,0);
	cor_=m_fopen(cor_name_,0);
	cigar_=m_fopen(cigar_name_,0);
	add_=m_fopen(add_name_,0);
}
void FQZip_decompression(char *fastq_name, char *compressedfile, int assemble_flag,int highestCom_flag)
{


	char meta_nameC[100]={'\0'},qs_nameC[100]={'\0'};
	char pos_nameC[100]={'\0'},cigar_nameC[100]={'\0'};
	char cor_nameC[100]={'\0'},add_nameC[100]={'\0'};
	strcpy(meta_nameC,fastq_name);
	m_strcat(meta_nameC,".meta.txt");
	strcpy(qs_nameC,fastq_name);
	m_strcat(qs_nameC,".qs.txt");
	strcpy(pos_nameC,fastq_name);
	m_strcat(pos_nameC,".pos.txt");
	strcpy(cigar_nameC,fastq_name);
	m_strcat(cigar_nameC,".cigar.txt");
	strcpy(cor_nameC,fastq_name);
	m_strcat(cor_nameC,".cor.txt");
	strcpy(add_nameC,fastq_name);
	m_strcat(add_nameC,".add.txt");
	strcpy(fasta_nameC,fastq_name);
	m_strcat(fasta_nameC,".fasta");

	char meta_lz[100]={'\0'},qs_lz[100]={'\0'};
	char pos_lz[100]={'\0'},cigar_lz[100]={'\0'};
	char cor_lz[100]={'\0'},add_lz[100]={'\0'};
	char fasta_lz[100]={'\0'};
	char order[1000]={'\0'};
	char dir[500] = { '\0' };
	char *p;
	int i;
	bool t = 0;
	strcpy(dir, fastq_name);
	for (i = 0; i < 256; i++)
	{
		if (dir[i] == '/')
		{
			t = 1; break;
		}
	}
	if (t == 1)   //not LWFQZip dir
	{
		p = rindex(dir, '/');
		for (i = 0; i < 30; i++)
			*p = '\0'; p++;
	}
	strcpy(meta_lz,fastq_name);
	m_strcat(meta_lz,".meta.txt.lz");
	strcpy(qs_lz,fastq_name);
	m_strcat(qs_lz,".qs.txt.lz");
	strcpy(pos_lz,fastq_name);
	m_strcat(pos_lz,".pos.txt.lz");
	strcpy(cigar_lz,fastq_name);
	m_strcat(cigar_lz,".cigar.txt.lz");
	strcpy(cor_lz,fastq_name);
	m_strcat(cor_lz,".cor.txt.lz");
	strcpy(add_lz,fastq_name);
	m_strcat(add_lz,".add.txt.lz");
	strcpy(fasta_lz,fastq_name);
	m_strcat(fasta_lz,".fasta.lz");
	if (t == 1)
	{
		sprintf(order, "%s %s %s %s", "tar -x -f ", compressedfile, " -C ", dir);
	}
	if (t == 0)
	{
		sprintf(order, "%s %s ", "tar -x -f ", compressedfile);
	}
	if (-1 == (system(order)))printf("Please make sure you have installed the tar tool correctly");
	if(assemble_flag)//assemble-based model
	sprintf(order,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %d ","./FQZip d  ",pos_lz,cigar_lz,cor_lz,meta_lz,add_lz,qs_lz,pos_nameC,cigar_nameC,cor_nameC,meta_nameC,add_nameC,qs_nameC,fasta_nameC,highestCom_flag);
	else
	sprintf(order,"%s %s %s %s %s %s %s %s %s %s %s %s %s %d ","./FQZip d  ",pos_lz,cigar_lz,cor_lz,meta_lz,add_lz,qs_lz,pos_nameC,cigar_nameC,cor_nameC,meta_nameC,add_nameC,qs_nameC,highestCom_flag);
	if (-1 == (system(order)))printf("Please make sure you have installed the FQZip tool correctly");
	remove(meta_lz);
	remove(qs_lz);
	remove(pos_lz);
	remove(cigar_lz);
	remove(cor_lz);
	remove(add_lz);
	if(assemble_flag)
	{remove(fasta_lz);
	}
	

}


void recover_base(char *ref_name_)
{
	size_t buffer_size=get_file_size(ref_name_),buf_count=0,fread_num=0,i=0,j=0;;
	buf_count = buffer_size+200000;
	char *ref_temp=m_malloc(buf_count);
	char *ref_gen=m_malloc(buf_count),ref_ch;
	memset(ref_temp,'\0',sizeof(char)*buf_count);
	memset(ref_gen,'\0',sizeof(char)*buf_count);

	fgets(ref_temp,200000,ref_g_);
	j=0;
	while((fread_num=fread(ref_temp,sizeof(char),buf_count, ref_g_))){
		for(i=0;i<fread_num;i++){
			if( ( ref_ch=m_toupper(*(ref_temp+i)) )>64){
				*(ref_gen+j)=ref_ch;
				j++;
			}
		}
	}
	free(ref_temp);ref_temp=NULL;
	//printf("load reference end.\n");

	char _pos[30]={'\0'},_cigar[500]={'\0'},_cor[5000]={'\0'};
	char *_add=m_malloc(read_len);
	char *base_temp=m_malloc(read_len);
	char *temp=m_malloc(read_len);
	memset(_add,'\0',sizeof(char)*read_len);
	memset(base_temp,'\0',sizeof(char)*read_len);
	memset(temp,'\0',sizeof(char)*read_len);

	while(fgets(_pos,30,pos_)&&fgets(_cigar,500,cigar_)&&fgets(_add,read_len,add_)&&fgets(_cor,5000,cor_))
	{
		new_get_base(base_temp,_add,_cigar,ref_gen,_pos,temp);
		variation_recover(base_temp,_cor);
		if(_cigar[0]=='0'){
			palindrome_conversion(temp,base_temp);
			fputs(temp,base_);
			fputc('\n',base_);
			memset(temp,'\0',strlen(temp));
		}
		else{
			fputs(base_temp,base_);
			fputc('\n',base_);
		}
		memset(base_temp,'\0',strlen(base_temp));
	}
	fclose(base_);
	free(ref_gen); ref_gen=NULL;
	free(base_temp); base_temp=NULL;
	free(_add); _add=NULL;
	free(temp); temp=NULL;
}

char *variation_recover(char *base_temp,char *_cor)
{
	size_t i=0,j=0,len=0,next=0;
	char NUM[30]={'\0'},_base;

	if(_cor[0]=='\n')
	{
		return base_temp;
	}
	else
	{
		for(i=0;i<strlen(_cor);i++)
		{
			if(_cor[i]>='0'&&_cor[i]<='9')
				NUM[j++]=_cor[i];
			if(_cor[i]>='A'&&_cor[i]<='Z')
			{
				_base=_cor[i];
				len=atoi(NUM);
				next=len+next;
				base_temp[next]=_base;
				j=0;
				memset(NUM,'\0',30);
			}
		}
		return base_temp;
	}
}

char *palindrome_conversion(char *des,char *str)
{
	size_t i=0,j=0;
	char temp;
	j=strlen(str)-1;
	while(j>i)
	{
		temp=str[i];
		str[i]=str[j];
		str[j]=temp;
		j--;
		i++;
	}
	for(i=0;i<strlen(str);i++)
	{
		switch(str[i])
		{
			case 'C':
				{
					str[i]='G';
					break;
				}
			case 'G':
				{
					str[i]='C';
					break;
				}
			case 'T':
				{
					str[i]='A';
					break;
				}
			case 'A':
				{
					str[i]='T';
					break;
				}
			default:break;
		}
	}
	strcpy(des,str);

	return des;
}

char *get_meta_variation_len(char *len,char *meta)
{
	sscanf(meta,"%s",len);
	return len;
}

char *get_meta_variation(char *variation,char *meta)
{
	sscanf(meta,"%*s%s",variation);
	return variation;
}

char *get_sra_id(char *id,char *str)
{
	int i=0,j=0;
	int isLocal=1;
	char *p;
	for(i=0;i<50;i++)
	{
	if (str[i]=='/')
	{isLocal=0;break;}
	}
	if(isLocal)
	{
		for(i=0;i<50;i++)
		if(str[i]=='.')	// local
			break;
		else id[j++]=str[i];
	}
	else
	{
		p = rindex(str, '/');p++;
		for (j = 0; j < 100; j++)
		{
			if(*p=='.')
				break;
			id[j]=*p;
			p++;
		}	
	}
	//printf("sra+id=%s\n",id);
	return id;
}

void recover_meta(char *fastq_name)
{
	char sra_id[50]={'\0'};
	get_sra_id(sra_id,fastq_name);

	char meta_temp[1000]={'\0'},pre_meta[1000]={'\0'},mid_meta[1000]={'\0'};

	size_t buffer_size=get_file_size(meta_name_);
	size_t buf_count=0,fread_num=0,read_num=1;

	if(buffer_size<block_size)
		buf_count = buffer_size;
	else buf_count = block_size;

	char *buffer=m_malloc(buf_count+read_len);
	memset(buffer,'\0',sizeof(char)*(buf_count+read_len));

	//Dispose the first metadata
	char variation_len[5]={'\0'},variation[1000]={'\0'},num[10]={'\0'};
	int len=0,m=0,n=0,i=0,j=0,snum=-1;

	fgets(meta_temp,1000,meta_);
	sscanf(meta_temp,"%*s%s",pre_meta);
	meta_temp[strlen(meta_temp)-1]='\0';
	fputs(meta_temp,meta_data);
	fputs(" length=\n",meta_data);

	memset(meta_temp,'\0',strlen(meta_temp));
	while((fread_num=fread(buffer,sizeof(char),buf_count, meta_)))
	{
		buffer[fread_num]='\0';
		for(i=0;i<fread_num;i++)
		{
			*(meta_temp+j)=*(buffer+i);
			j++;
			if(*(buffer+i)=='\n')
			{
				read_num++;
				sscanf(meta_temp,"%s",variation_len);
				//get_meta_variation_len(variation_len,meta_temp);
				snum=sscanf(meta_temp,"%*s%s",variation);
				len=atoi(variation_len);

				fputc('@',meta_data);
				fputs(sra_id,meta_data);
				fputc('.',meta_data);
				sprintf(num,"%zu",read_num);
				fputs(num,meta_data);
				fputc(' ',meta_data);
				memset(num,'\0',strlen(num));

				if(len==0){
					strcpy(mid_meta,variation);
					strcpy(pre_meta,variation);
				}
				else{
					for(m=0,n=0;m<len;m++,n++){
						mid_meta[m]=pre_meta[n];
					}
					if(snum>=0)
						m_strcat(mid_meta,variation);
					memset(pre_meta,'\0',strlen(pre_meta));
					strcpy(pre_meta,mid_meta);
				}
				fputs(mid_meta,meta_data);
				char *ptr=NULL;
				ptr=strstr(mid_meta,"length");
				if(ptr==NULL)
					fputs(" length=\n",meta_data);
				else fputs("length=\n",meta_data);
				j=0;
				memset(mid_meta,'\0',strlen(mid_meta));
				memset(meta_temp,'\0',strlen(meta_temp));
				//                memset(variation,'\0',strlen(variation));
				//                memset(variation_len,'\0',strlen(variation_len));
			}
		}
	}
	free(buffer);buffer=NULL;
	fclose(meta_data);
}

void quit_decompression_dispose(int assemble_flag)
{
	fclose(meta_);
	//fclose(qs_);
	fclose(ref_g_);
	fclose(pos_);
	fclose(cigar_);
	fclose(add_);
	fclose(cor_);
	remove(meta_name_);
	remove(qs_name_);
	remove(pos_name_);
	remove(cigar_name_);
	remove(cor_name_);
	remove(add_name_);
	remove(base_name_);
	remove(metadata);
	if(assemble_flag)
	remove(fasta_nameC);
	//remove(qscore_name);
}



void *thread_base(void *args)
{
	recover_base(ref_);
}
void *thread_meta(void *args)
{
	recover_meta(fastq_name);
}
void *thread_qs(void *args)
{
	recover_qs();
}



void LWFQZip_decompression(char *input_name,char *ref_name,int assemble_flag,int highestCom_flag)
{
	char *ptr=NULL;
	ptr=input_name;
	if(ptr==NULL){
		fprintf(stderr,"Please input correct LWFQZip compressed file end with 'lw or lz '\n");
		exit(EXIT_FAILURE);
	}
	ref_ = ref_name;
	get_fastq_name(fastq_name,input_name);
	FQZip_decompression(fastq_name,input_name,assemble_flag,highestCom_flag);
	init_decompression(fastq_name,ref_name);
	pthread_create(&t1, 0, &thread_base, NULL);
	pthread_create(&t2, 0, &thread_meta,NULL);
	pthread_join(t1, 0);
	pthread_join(t2, 0);
	//recover_qs();
	recover_fastq(fastq_name);
	quit_decompression_dispose(assemble_flag);
}




//L is the length of each read--L<0 means var length
char* decode_quality_score(char *Instream, FILE *outfp)
{
	int i,Outsize=0, flag=1;
	char nextch, ch, tempch, prech=0, Outstr[LINESIZE]={0};
	//LineCount++;
	//prech=*Instream; //init
	//i=1;
	while( (ch=*Instream++ )>32 ) //EOF == 0 '\0' ? condition: end of line is 32, end of file is 0
	{
		//Outsize+=sprintf(Outstr+Outsize,"%c",ch);
		//Outstr[Outsize++]=ch;
		if((nextch=*Instream++)>32)
		{
			Outstr[Outsize++]=ch;
			prech=ch;		flag=1;
		}
		else if(nextch<32)//next is a number
		{
			do
			{
				if( nextch>0 ) //1-31
				{
					for(i=nextch; i>0; i--)
						Outstr[Outsize++]=ch;
				}
				else if(nextch<0) //32-127 negative
				{
					for(i=-nextch; i>0; i--)
						Outstr[Outsize++]=ch;
				}
				else // ==0 last character
				{
					if(flag) //pre is a symbol
						Outstr[Outsize++]=ch;  //last character:1 last is number 2 last is a symbol
					break;
				}
				tempch=ch;
				ch=prech;
				prech=tempch; flag=0;
			}while( (nextch=*Instream++)<32 );
			if(nextch==32) //output
			{
				Outstr[Outsize++]='\n';
				fwrite(Outstr,1,Outsize,outfp); //Outsize=L+1
				break;
			}
		}
		else //32 is end  output
		{
			Outstr[Outsize++]=ch;
			Outstr[Outsize++]='\n';
			//fprintf(outfp,"%d ",LineCount);
			fwrite(Outstr,1,Outsize,outfp); //Outsize=L+1
			break;
		}
		Instream--;
	}

	//if(LineCount==66629)printf("66629: outsize:%d\n%s",Outsize,Outstr); //test
	return Instream; //return the changed Instream
}

void recover_qs()
{
	char *codebuf  =(char*)calloc(BLOCK*sizeof(char), sizeof(char)); //enough size
	char *ptr=codebuf;
	int i=0;

	size_t filesize=0;
	while( filesize=fread(codebuf,1,BLOCK,qs_) ) //canbe several blocks
	{
		//endptr=codebuf+cfilesize-1; //last character
		i=0;
		while(*(codebuf+filesize-1-i) !=32) i++;
		*(codebuf+filesize-i)='\0';
		fseek(qs_, -i, SEEK_CUR);
		//printf("%d\n",*(codebuf+filesize-1-i));
		while(*codebuf>32)
		{
			codebuf=decode_quality_score(codebuf,qscore);
		}
		codebuf=ptr;
	}
	free(ptr);ptr=NULL;
	fclose(qscore);
}

void recover_fastq(char *fastq_name)
{
	strcpy(de_fastq_name,fastq_name);
	m_strcat(de_fastq_name,".d");
	fastq_=m_fopen(de_fastq_name,1);

	base_=m_fopen(base_name_,0);
	meta_data=m_fopen(metadata,0);
	qscore=m_fopen(qscore_name,0);

	char *meta_temp=m_malloc(1000);
	char *base_temp=m_malloc(read_len);
	char *qs_temp=m_malloc(read_len);
	memset(meta_temp,'\0',sizeof(char)*1000);
	memset(base_temp,'\0',sizeof(char)*read_len);
	memset(qs_temp,'\0',sizeof(char)*read_len);
	size_t readlen=0;
	char num[10]={'\0'};

	while(fgets(meta_temp,1000,meta_data)&&fgets(base_temp,read_len,base_)&&fgets(qs_temp,read_len,qscore))
	{
		readlen=strlen(qs_temp);
		meta_temp[strlen(meta_temp)-1]='\0';
		sprintf(num,"%zu",readlen-1);
		m_strcat(num,"\n");
		m_strcat(meta_temp,num);
		fputs(meta_temp,fastq_);
		memset(num,'\0',strlen(num));

		if(readlen==strlen(base_temp)){
			fputs(base_temp,fastq_);
		}
		else{
			base_temp[strlen(base_temp)-1]='\0';
			fputs(base_temp,fastq_);
			fgets(base_temp,read_len,base_);
			fputs(base_temp,fastq_);
		}
		meta_temp[0]='+';
		fputs(meta_temp,fastq_);

		fputs(qs_temp,fastq_);
	}
	fclose(fastq_);
	fclose(base_);
	fclose(meta_data);
	fclose(qscore);
	free(meta_temp);meta_temp=NULL;
	free(base_temp);base_temp=NULL;
	free(qs_temp);qs_temp=NULL;
}
