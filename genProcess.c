//load the genome file to memery and get the index of PreKey!

#include "head.h"
//#define INFINITY 0x7fffffff
char  *file_name=NULL,  *gen_buf=NULL;  
int   **Indexarray=NULL, *IndexarrayNum=NULL, *InitarraySize=NULL;  //store the IndexPosition

int isfasta(char * p)
{
	int len, i;
	len=strlen(p);
	for(i=len-1; *(p+i)!='.'; i--);
	if( !stringcmp(p+i, ".fasta") || !stringcmp(p+i, ".fa") )
		return 1;
	else return 0;
}
void AddIndexWord(unsigned int Key, int Posi) //Link or array
{
	//Link
	/*struct IndexLink *p;
	  p=(struct IndexLink *) malloc( sizeof(struct IndexLink) );
	  p->Posi=posi;
	  p->Next=*(IndexData+Key)->Next;
	 *(IndexData+Key)->Next=p;
	 return *(IndexData+Key);*/
	//array

}
char * getIndex(int NumberP, int *IndexP, char *buf) //get the Index via Hashtable
{
	int i, j, KeyP, PreKeyLength, IndexHashLength;
	unsigned int  HashWord=0;
	PreKeyLength=strlen(IndexPreKey);
	IndexHashLength=_pow(4,HashStrlen);
	InitarraySize=(int*)calloc( IndexHashLength,sizeof(int) );
	IndexarrayNum=(int*)calloc( IndexHashLength,sizeof(int) );
	Indexarray=(int **) calloc( IndexHashLength,sizeof(int*));

	for(i=0; i<IndexHashLength; i++)
	{
		*(InitarraySize+i)=100;
		*(Indexarray+i)=(int *) calloc( *(InitarraySize+i), sizeof(int) );
	}

	for(i=0; i<NumberP; i++) 
	{
		KeyP=*(IndexP+i);
		for(j=0; j<HashStrlen; j++)   //hashStrlen=8
			HashWord=base_hash( *(buf+KeyP+j),HashWord );
		*(*(Indexarray+HashWord)+ *(IndexarrayNum+HashWord) )= KeyP;  //posi is point to the next char of "PreKey"

		*(IndexarrayNum+HashWord)+=1;
		if( *(IndexarrayNum+HashWord)==*(InitarraySize+HashWord) )  //request more memory
		{
			*(InitarraySize+HashWord) += 100;
			*(Indexarray+HashWord)=(int *)realloc( *(Indexarray+HashWord), sizeof(int)*(*(InitarraySize+HashWord)) );
		}
		HashWord=0;
	}
	return buf;
}
size_t InitsrcFile(FILE *fp)
{
	char firstch, str[2048];
	size_t genome_size=0;
	int  size;		
	fseek(fp,0,SEEK_END);
	genome_size=ftell(fp);  //get file size
	fseek(fp,0,SEEK_SET);

	if( (firstch=fgetc(fp))=='>' )       //if there is a title
	{	
		title=(char *)calloc( 4096,sizeof(char) );
		fscanf(fp,"%s",title);
		fgets(str,2048,fp); 		
	}
	//printf("%s\n",*title );  //test title
	size=sizeof(char);
	fseek(fp, -size, SEEK_CUR); //back a char	

	return genome_size;
}

//load genome to genBaseBuf & getIndex
char * genProcess(char *genFileName, char *BaseBuf, char *PreKey) 
{	
	FILE  *srcFile,  *outFile;
	size_t genSize=0;
	long int   i=0, j=0, base_count=0, readtemp;
	long int   Key_offset,offset, PreKeyLength, IndexPosi_count, InitIndexPosi_count;
	char  basech, PreKeyEndchar;

	//open file
	srcFile=fopen(genFileName,"rb");
	if(!srcFile)
	{
		printf("error input genome file :%s\n",genFileName);
		return NULL;
	}
	//read file
	genSize=InitsrcFile(srcFile);
	if(genSize>0)
	{
		gen_buf=(char *)malloc(genSize);
		BaseBuf=(char *) calloc( genSize, sizeof(char) );
		if(!gen_buf || ! BaseBuf) { printf("lack of memory to load the genomes!\n"); return NULL;  }
		readtemp=fread(gen_buf,sizeof(char),genSize,srcFile); //if readtemp < genSize  ERROR
	}

	else { printf("error genome file size\nMax Genome Size :4G\n"); return NULL;	}
	//test readtemp & genSize
	BaseBuf=BaseBuf+5;  //the start is 5 zero! for mapping

	PreKeyLength=strlen(PreKey);
	PreKeyEndchar=Caps(*(PreKey+PreKeyLength-1));
	Key_offset=PreKeyLength-2;
	offset=2;
	IndexPosi_count=0;   InitIndexPosi_count=128*1024;
	IndexPosition=(int*)calloc( InitIndexPosi_count, sizeof(int) );
	for(j=0; j<readtemp; j++)
	{
		if( ( basech=Caps(*(gen_buf+j)) )>64 )  //del the LF
		{
			*(BaseBuf+base_count)=basech;  base_count++;  //record base			
			//get index depend on PreKey
			if( basech == PreKeyEndchar ) //request PreKeyLength>1 
			{
				while( *(BaseBuf+base_count-offset)==*(PreKey+PreKeyLength-offset) )
				{					
					if(PreKeyLength-offset==0)
					{
						*(IndexPosition+IndexPosi_count)=base_count; //+1 means the next char of "CG..."

						IndexPosi_count++;
						if(IndexPosi_count==InitIndexPosi_count)  //request more memory
						{
							InitIndexPosi_count+=InitIndexPosi_count;
							IndexPosition=(int*)realloc(IndexPosition, InitIndexPosi_count*sizeof(int));
						}
						break;
					}						
					offset++;
				}
				offset=2;
			}

		}
	}
	while( *(IndexPosition+IndexPosi_count-1) > base_count-20 ) IndexPosi_count--;  //delete the last 

	//record the information of struct Ind
	BaseCount = base_count;       //base_num[0]=0 
	IndexPosiCount = IndexPosi_count; 
	//HashStrlen = 8;   //8 or 10 14
	strcpy( IndexPreKey , PreKey );

	fclose(srcFile);  
	if(gen_buf!=NULL) {	free(gen_buf); gen_buf=NULL;  }
	//get index from IndexPosition
	BaseBuf=getIndex(IndexPosi_count, IndexPosition, BaseBuf);
	return BaseBuf;
}
