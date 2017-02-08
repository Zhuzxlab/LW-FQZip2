//align function:  align the read short to the genomes like the BWA progress
#include "head.h"

//#define DEBUG

void mainUsage(void);


int  *base_num=NULL, *KeyofPosiCount=NULL, file_count=0; 
char *title=NULL,    *genBaseBuf=NULL,   IndexPreKey[16]={0};  
int  HashStrlen=8,  IndexPosiCount=0,  BaseCount=0, *IndexPosition=NULL,VALUELEN=12,cp=1; //Index information
float E=0.05;
//int CONSIDERLEN;
int main(int argc,char *argv[]) 
{	
	//CONSIDERLEN=E*200;
	int i,j, next=0, str_len=0,  Readmode=1;
	int	d=0,  f=0; //d, f is parameters for map function
	char genFileDir[256]={0}, readFileDir[256]={0}, outFile[256]={0},INDEXKEY[10]={0}; //index num <= 10;
	char *testargv[3]={"cmp", "srcFile", "readFile" };

	//for(j=0; j<10; j++)  //init the INDEXKEY
	//	INDEXKEY[j]=(char *)calloc(256,sizeof(char));


#ifdef DEBUG //if debug had defined
	//genomeRead("db.fa", BaseBuffer);
	//ind=getIndex(BaseBuffer, "CG");
	genBaseBuf=genProcess("NC_017634.1.fasta",genBaseBuf,"CG"); //15829254.fasta  NC_000913.3.fasta  NC 001136.10
	// OutputIndex(ind);  //ind is a struct
	map("SRR1063349.fastq", genBaseBuf,"result.txt",testargv); //map read to db, output the result to read.sam  ERR231645.fastq

#else
	//getParameter( argc,argv );
	if(argc<3) { printf("wrong parameter!\n");  mainUsage(); return -1;}
	strcpy(genFileDir,argv[1]);
	strcpy(readFileDir,argv[2]);
	strcpy(outFile,argv[2]);
	strcat(outFile,".map.txt");
	strcpy(INDEXKEY,"CG");
	i=3;
	while(argc>i) 
	{					
		if(*argv[i]=='-')  //option: -P -k -h 
		{
			switch ( *(argv[i]+1) )
			{
				case 'P':
				case 'p': i++;   //printf("P : %s\n",argv[i]); //test
					  strcpy(INDEXKEY,argv[i]);   
					  for(j=0; INDEXKEY[j]!='\0'; j++)
						  INDEXKEY[j]=Caps(INDEXKEY[j]);  //convert to caps
					  break;  

				case 'K':
				case 'k': i++;   //printf("k: %s\n",argv[i]); //test
					  HashStrlen = atoi(argv[i]); 
					  break;

				case 'E':
				case 'e': i++; E = atof(argv[i]); 

					  break;

				case 'L':
				case 'l': i++; //printf("L: %s\n",argv[i]);
					  VALUELEN = atoi(argv[i]);					
					  break;
				
				case 'O':
				case 'o': i++; //printf("L: %s\n",argv[i]);
					  cp = atoi(argv[i]);
					   break;
			}
			i++;
		}	
		else break;
	}	
	//genomeRead(genFileDir,BaseBuffer);  //genRead and getIndex should in the same time
	if( (genBaseBuf=genProcess(genFileDir,genBaseBuf,INDEXKEY))==NULL)
		return -1;
	map(readFileDir,genBaseBuf,outFile, argv);


#endif



#ifdef DEBUG //if debug had defined-- temp output 


	//#elif


#endif
	return 0;
}

void mainUsage(void)
{
	printf("  cmd <db.fa> <read.fq>\n");
	printf("\n  get help by 'cmd -h'\n");
	printf("Mapping Options Options:\n");
	printf("  -b, --the number of mapping thread(Default: 10)\n");
	printf("  -p, --specify the kmer prefixes, e.g.,'CG', 'AT', and 'TAG' (Default: '-p CG'). 'AA' is not recommended as a prefix.\n");
	printf("  -k, --length of a kmer used in locate local alignment. (Default: '-k 8')\n");		
	printf("  -e, --the tolerance rate of mismatches.(Default: '-e 0.05')\n");
	printf("  -L, --the mini length of a legal alignment.(Default: '-l 12')\n");
	printf("  -o, --open the complementart palindrome mode.(Default: '-o 1' means start,otherwise'-o 0')\n");
}
