#include "head.h"

//map the short read to the buf_sequence
#define READLEN     300000  //1024*80
#define RBUFSIZE    104857600  //100M
#define block_size 2000000
//#define VALUELEN    30       //value length for match result
//#define CONSIDERLEN 10     //consider length when a wrong match&& insert or delete length
#define INFINITY    0x7fffffff  //max int value
char    *ReadBuf=NULL;    //block buffer: read use fread();

int  ReadPreKeyarray[READLEN]={0},RvP[READLEN]={0},Word[READLEN]={0},RWord[READLEN]={0},*pWord[READLEN]={NULL},*pEnd[READLEN]={NULL};  //store the PreKey position&HashWord in the short read
int  outFLAG,tempCiagr[10000],getLengthout[10000],out[2][READLEN]={0},Rout[2][READLEN]={0},exactMatch=0;//out[]&Rout[] store the map position info
long long mappedNum,doubleMap,read_count,Wsize,RWsize;
char outCiagr[10000], readTitle[10000]={0};
int  setNum,ReadLength,repeat2[4], setArray[READLEN]={0},setOrder[READLEN]={0}; //maybe the could be a struct datatype
unsigned long long totalbase=0, mappedbase=0;
int CONSIDERLEN,cp;
//pos.txt,cigar.txt,add.txt,cor.txt
FILE *pos=NULL,*cigar=NULL,*add=NULL,*cor=NULL,*ref_g=NULL,*qs=NULL;
char pos_name[500]={'\0'},cigar_name[500]={'\0'},cor_name[500]={'\0'},add_name[500]={'\0'},qs_name[500]={'\0'};
char add_temp[READLEN]={'\0'},clip_seq[READLEN]={'\0'},base_temp[READLEN]={'\0'},match_temp[READLEN]={'\0'},base_[READLEN]={'\0'};
int a=1;

void* m_malloc(size_t n)
{
	void* p = malloc(n);
	if (p == NULL) {
		fprintf(stderr, "Can not allocate %zu bytes.\n", n);
		exit(EXIT_FAILURE);
	}
	return p;
}

char *m_strcat(char *dst, const char* src)
{
	while(*dst)   dst++;
	while((*dst ++ = *src++));
	return --dst ;
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

long int get_file_size(char *filename){

	struct stat buf;
	stat (filename, &buf);
	long int filesize;
	filesize = buf.st_size;
	return filesize;
}

char m_toupper(char ch)
{
	if (ch>96)	return ch-32;
	else return ch;
}

void MemeryArray(int a[], int n, int b[], int m, int c[])  //merge sort
{
	int i, j, k, anum, bnum;
	i = j = k = 0;
	anum=Word[n];  bnum=Word[m];
	while (i < anum && j < bnum)
	{
		if (a[i]-ReadPreKeyarray[n] < b[j]-ReadPreKeyarray[m])
			c[k++] = a[i++]-ReadPreKeyarray[n];
		else
			c[k++] = b[j++]-ReadPreKeyarray[m];
	}

	while (i < anum)
		c[k++] = a[i++]-ReadPreKeyarray[n];

	while (j < bnum)
		c[k++] = b[j++]-ReadPreKeyarray[m];
}
void Fintersection(int a[], int aNum, int b[], int bNum)//find intersection
{
	int i,j;
	for(i=0; i<aNum; i++)
	{

	}
}
void InitCmpOrder(int *a, int n, int *b, int v)  //may use the struct data!!!////////
{
	int i,j, min, temp;
	for(j=0; j<n; j++)  //init cmpOrder
	{
		min=j;
		for( i=j+1;i<n; i++)
		{
			if( a[b[min]]>a[ b[i] ] )
				min=i;
		}
		//printf("%d\n",a[b[min]]);
		temp=b[j];
		b[j]=b[min];
		b[min]=temp;
	}

	temp=1; j=1; b[n+2]=1;
	for(i=1; i<v; i++) //find same
	{
		if(a[ b[i] ] - a[ b[i-1] ] < CONSIDERLEN )
			temp++;
		else //not same
		{
			min=i;
			if(j<temp)
			{	b[n]=temp;  b[n+1]=min; b[n+2]=min-temp; j=temp;}
			temp=1;
		}
	}
	if(j<temp || j==1)
	{	b[n]=temp;  b[n+1]=v; b[n+2]=v-temp; }
}

int partCompare(int *a,int n, int *b, int*p ) //find the mini value & return repeatCount
{
	int i=0, p0, p1, valueNum;
	int count, temp, maxcount;
	//if( (p=b[i])!=b[n+2] ) //b[0]==b[n+2]
	p0=b[0];
	valueNum=b[n+3];
	if( pWord[p0 ]<pEnd[p0 ])//insert a value to a[] and updata the b[]
	{
		a[p0]=*pWord[p0]-p[p0]; //offset
		pWord[p0]++;
		//InsertValue(a,n,b);
		for(i=1; i<valueNum; i++)  //compare the valueNum, not the n--contain INFINITY
		{
			p1=b[i];
			if(a[p0]>a[p1])
				b[i-1]=p1;
			else {  b[i-1]=p0; break; } //insert
		}
		if(i==valueNum) b[valueNum-1]=p0;
		//updata n++
		count=1; maxcount=1;
		for(i=1;i<valueNum; i++)
		{
			if(a[b[i] ]-a[b[i-1] ]<CONSIDERLEN) count++;
			else
			{
				temp=i;
				if(maxcount<count)
				{ b[n]=count; b[n+1]=temp; b[n+2]=temp-count; maxcount=count;}
				count=1;
			}
		}
		if(maxcount<count || maxcount==1)
		{
			b[n]=count; b[n+1]=valueNum; b[n+2]=valueNum-count;
		}
		//-------efficient ways to updata n++---

		//------------------------
	}
	else //no value
	{
		for(i=1; i<valueNum; i++)
			b[i-1]=b[i];
		--valueNum;
		b[valueNum]=p0;
		a[p0] =	INFINITY;
		b[n+3]=valueNum;

		if(b[n+2]>0)
		{	b[n+2]--; b[n+1]--;  }
		else //==0
		{
			if(b[n]>1)
				b[n]--;
			//else b[n+1]=1;
		}
	}
	return valueNum;
}

int	getSet( int n, int *p, int *W)  //map position
{
	int i, temp, maxrepeatTime=1;
	//find intersection
	temp=0;
	for(i=0; i<n; i++) //init cmpArray
	{
		pEnd[i]=(*(Indexarray+W[i])+ *(IndexarrayNum+W[i]));
		pWord[i]=*(Indexarray+W[i]);
		*(setOrder+i)=i;
		if( pWord[i]<pEnd[i]) //  *(IndexarrayNum+Word[i]) >0
		{
			*(setArray+i)=*pWord[i]-p[i];  //pop the first value
			pWord[i]++; //point next position
			temp++;
		}
		else { *(setArray+i)=INFINITY;	} //no value
	}
	setNum=temp; //init value set num
	if(temp<2) return temp;  //temp=0 or 1  means the repeat time!

	*(setOrder+n+3)=temp;
	InitCmpOrder(setArray,n,setOrder, setNum);//init cmpOrder
	while(temp>1)
	{
		if( setOrder[n]>2 ) return setOrder[n]; //repeat time 3 4 5

		if(maxrepeatTime<setOrder[n]) //record the result
		{
			maxrepeatTime=setOrder[n];
			i=setOrder[n+2];   //start cmpArray
			repeat2[0]=setArray[ setOrder[i] ]; //src start posi
			repeat2[1]=p[ setOrder[i] ];  //offset for read
			repeat2[2]=setArray[ setOrder[i+1] ];
			repeat2[3]=p[ setOrder[i+1] ];
		}
		if( temp==maxrepeatTime ) return temp; //2 ...
		temp=partCompare(setArray, n, setOrder,p);
	}
	return maxrepeatTime; //could be 1 or 2
}
int writeOutArray(int repeat, int n, int *Rp,int status)
{
	int s,l,temp,minOrder,maxOrder; //the size of o[] is 10
	int *o, *p;
	if(status)
	{	o=out[0];  p=out[1];   }
	else
	{   o=Rout[0]; p=Rout[1];  }
	//printf("repeat:%d\n",repeat);
	if( repeat>2)  //3 4 5
	{
		l=setOrder[n+1];
		s=setOrder[n+2];
		////-----sort Order[]; write o[1] and o[ setOrder[n] ]
		minOrder=setOrder[s]; maxOrder=setOrder[l-1];
		if(setArray[maxOrder]-setArray[minOrder] > CONSIDERLEN) //InDel Consider length
		{ o[0]=1;	 return 1;		}
		if(maxOrder<minOrder) //ensure the max>min
		{
			temp=maxOrder;
			maxOrder=minOrder;
			minOrder=temp;
		}
		o[1]=setArray[minOrder]+Rp[minOrder];
		p[1]=Rp[minOrder];
		o[repeat]=setArray[maxOrder]+Rp[maxOrder];
		p[repeat]=Rp[maxOrder];
	}
	//consider o[n]=1 or 0
	else if(repeat==2)
	{
		if(repeat2[1]<repeat2[3]) //ensure p1<p2
		{
			o[1]=repeat2[0]+repeat2[1];
			p[1]=repeat2[1];
			o[2]=repeat2[2]+repeat2[3];
			p[2]=repeat2[3];
		}
		else
		{
			o[1]=repeat2[2]+repeat2[3];
			p[1]=repeat2[3];
			o[2]=repeat2[0]+repeat2[1];
			p[2]=repeat2[1];
		}
	}

	o[0]=1;

	return repeat;
}
int getWord(char *R,int Rp[],int *W, int pCount)
{
	int i,j,temp, tempCount, HashWord;
	tempCount=pCount-1;
	temp=ReadLength-HashStrlen;
	while(Rp[tempCount]>temp) //delete the unvalue prekey
		tempCount--;
	tempCount++;
	for(i=0; i<tempCount; i++) //get all HashWord--tempCount
	{
		HashWord=0;
		temp    =Rp[i]+HashStrlen;

		for(j=Rp[i]; j<temp; j++)
			HashWord= base_hash( R[j],HashWord );

		W[i]=HashWord;
	}
	return tempCount;
}
int getMapPosi(int pCount, char *R, char *Rv)  //normal and reverse
{
	int repeatNum, RrepeatNum, i,j; //variaty for set
	if(pCount==0)
		return -2;
	//normal--get HashWord depend on ReadPreKeyarray[]
	Wsize=getWord(R, ReadPreKeyarray,Word,pCount);
	repeatNum=getSet(Wsize,ReadPreKeyarray, Word);	//get map position
	repeatNum=writeOutArray(repeatNum,Wsize,ReadPreKeyarray,1); //1 means normal

	//if(repeatNum>2) return 0; //return a FLAG

	//reverse--get HashWord depend on ReadPreKeyarray[]
	j=0;
	for(i=pCount-1;i>=0;i--)
		RvP[j++]=ReadLength-ReadPreKeyarray[i]+2; //+2 means the strlen(PreKey)
	RWsize=getWord(Rv,RvP,RWord,pCount);
	RrepeatNum=getSet(RWsize,RvP,RWord);
	RrepeatNum=writeOutArray(RrepeatNum,RWsize,RvP,0);
	//printf("R:%d  Rout00 Rout01 03:%d  %d  %d\n",RrepeatNum,*Rout[0],Rout[0][1],Rout[3]);
	if(cp == 0){
		RrepeatNum =1;
	}
	if(RrepeatNum>1 && repeatNum>1) //two match
		return 2;
	else if(RrepeatNum>1 && repeatNum==1) //reverse match only
		return 1;
	else if(RrepeatNum==1 && repeatNum>1) // normal match only
		return 0;
	else		//no match
		return -1;
}

int getLength(int p,int r,  char *R) //return match length
{//p is posi of src, r is posi of R, Readlen is the length of R
	int  Rstart, genstart, matchlen;
	//------------------------
	//int i,j; //test prekey
	//genstart=p-2;
	//Rstart=r-2;
	//for(i=0; i<10; i++)
	//{
	//	if( R[Rstart+i]!=genBaseBuf[genstart+i] &&R[Rstart+i]!='N' )
	//	{
	//		printf("i:%d NO:%5d %d %d:",i,read_count,outFLAG,out[0][0]);
	//		printf("%c%c%c %c%c%c\n",R[Rstart+i],R[Rstart+i+1],R[Rstart+i+2],genBaseBuf[genstart+i],genBaseBuf[genstart+i+1],genBaseBuf[genstart+i+2] );
	//		break;
	//	}
	//}
	//------------------------------
	genstart=p+HashStrlen;
	for(Rstart=r+HashStrlen; Rstart<ReadLength; Rstart++, genstart++) //forward
	{
		if(R[Rstart]!=genBaseBuf[genstart])  //unmatch: judge the next char
		{
			if(R[Rstart+1]!=genBaseBuf[genstart+1]) //contain for N && R[Rstart]!='N'
				break;
			else { Rstart++; genstart++; }
		}
	}
	//if(Rstart==Readlen-1) Rstart++; //the last char unmatch!!!
	p=p-2;  //backward
	for(r=r-2; r>0; r--, p--)// r-3 request strlen(INDEXPRE)=2
	{
		if(R[r]!=genBaseBuf[p])  //unmatch: judge the next char
		{
			if(R[r-1]!=genBaseBuf[p-1])  //r-1 means r>=1
				break;
			//else { r--; p--; }
		}
	}
	if(r>0) { r++; p++; } //from r++ to match---if r==0 need not r++
	getLengthout[0]=p;
	getLengthout[1]=r;
	getLengthout[2]=Rstart;
	if( (matchlen=Rstart-r)==ReadLength )
		exactMatch++;

	return Rstart-r;	//match length
}
void storetempCiagr(int p, int r, int Rstart)
{
	int Dvalue=0;
	if(tempCiagr[0]<0) //no position
	{
		tempCiagr[0]=p;
		tempCiagr[1]=r;
		tempCiagr[2]=Rstart;  //read match end position
		tempCiagr[4]=-1; //means there is no next match
	}
	else //cmpare the two result: maybe the is a I Or D
	{//tempCiagr[3] between match and match is 3 probability 0:unmatch  -:I  +:D
		tempCiagr[4]=r;
		tempCiagr[5]=Rstart;

		Dvalue=tempCiagr[0]-tempCiagr[1]; //small posi
		tempCiagr[3]=p-r-Dvalue; //small+Dvalue-big
		//if(r>tempCiagr[1])
		//{
		//	tempCiagr[4]=r;
		//	tempCiagr[5]=Rstart;

		//	Dvalue=tempCiagr[0]-tempCiagr[1]; //small posi
		//	tempCiagr[3]=p-r-Dvalue; //small+Dvalue-big
		//}
		//else
		//{
		//	//test
		//	printf("getlength NO:%d ,1:%d, 4:%d\n",read_count,tempCiagr[1],tempCiagr[4]);
		//	//------------
		//	Dvalue=tempCiagr[0]-tempCiagr[1]; //big posi
		//	tempCiagr[0]=p;
		//	tempCiagr[4]=tempCiagr[1];
		//	tempCiagr[5]=tempCiagr[2];
		//	tempCiagr[1]=r;
		//	tempCiagr[2]=Rstart;
		//
		//	tempCiagr[3]=Dvalue-p+r; //big posi-small posi==if(t[3]>0) D; else I;
		//}
		//tempCiagr[6]=0; -------------------
		//if(tempCiagr[3]>CONSIDERLEN) //test setOrder: delete tempCiagr[4]
		//{
		//	printf("read_count:%d  tempC[1-5] %d--%d--%d--%d--%d\n",read_count,tempCiagr[1],tempCiagr[2],tempCiagr[3],tempCiagr[4],tempCiagr[5]);
		//	tempCiagr[4]=-1;
		//}
	}

}
void writeCIAGR(int Readlen) //from tmepCiagr to outCiagr
{
	int outCiagrlen=0,temp, temp2;
	//temp=tempCiagr[4]-tempCiagr[2];
	if(tempCiagr[4]<0)
	{
		if((temp=Readlen-tempCiagr[2])<2)
		{
			//tempCiagr[2]=Readlen;
			if (tempCiagr[1] < 2)//c[1]>0
			{
				outCiagrlen = sprintf(outCiagr, "%dM", Readlen);
				mappedbase += Readlen;
			}
			else
			{
				outCiagrlen = sprintf(outCiagr, "%dS%dM", tempCiagr[1], Readlen - tempCiagr[1]);
				mappedbase += (Readlen-tempCiagr[1]);
			}
		}
		else
		{
			if (tempCiagr[1] < 2)//c[1]>0
			{
				outCiagrlen = sprintf(outCiagr, "%dM%dS", tempCiagr[2], temp);
				mappedbase += tempCiagr[2];
			}
			else
			{
				outCiagrlen = sprintf(outCiagr, "%dS%dM%dS", tempCiagr[1], tempCiagr[2] - tempCiagr[1], temp);
				mappedbase += (tempCiagr[2] - tempCiagr[1]);
			}
		}
	}
	else//c[4]>0 means there is two
	{
		if (tempCiagr[1] < 2)//c[1]>0
		{
			outCiagrlen = sprintf(outCiagr, "%dM", tempCiagr[2]);
			mappedbase += tempCiagr[2];
		}
		else
		{
			outCiagrlen = sprintf(outCiagr, "%dS%dM", tempCiagr[1], tempCiagr[2] - tempCiagr[1]);
			mappedbase += (tempCiagr[2] - tempCiagr[1]);
		}

		if( (temp=Readlen-tempCiagr[5])<2)
			tempCiagr[5]=Readlen;

		if(tempCiagr[3]==0)
		{
			outCiagrlen+=sprintf(outCiagr+outCiagrlen,"%dS%dM",tempCiagr[4]-tempCiagr[2],tempCiagr[5]-tempCiagr[4]);
			mappedbase += (tempCiagr[5] - tempCiagr[4]);
		}
		else if(tempCiagr[3]>0) //D
		{
			if ((temp = tempCiagr[4] - tempCiagr[2]) > 0)
			{
				outCiagrlen += sprintf(outCiagr + outCiagrlen, "%dD%dS%dM", tempCiagr[3], temp, tempCiagr[5] - tempCiagr[4]);
				mappedbase += (tempCiagr[5] - tempCiagr[4]);
			}
			else
			{
				outCiagrlen += sprintf(outCiagr + outCiagrlen, "%dD%dM", tempCiagr[3], tempCiagr[5] - tempCiagr[2]);
				mappedbase += (tempCiagr[5] - tempCiagr[2]);
			}
		}
		else //<0 : I
		{
			temp2=-tempCiagr[3];
			if ((temp = tempCiagr[4] - tempCiagr[2] - temp2) > 0)
			{
				outCiagrlen += sprintf(outCiagr + outCiagrlen, "%dI%dS%dM", temp2, temp, tempCiagr[5] - tempCiagr[4]);
				mappedbase += (tempCiagr[5] - tempCiagr[4]);
			}
			else
			{
				outCiagrlen += sprintf(outCiagr + outCiagrlen, "%dI%dM", temp2, tempCiagr[5] - tempCiagr[2] - temp2);
				mappedbase += (tempCiagr[5] - tempCiagr[2] - temp2);
			}
		}
		if( Readlen-tempCiagr[5]>1)
			outCiagrlen+=sprintf(outCiagr+outCiagrlen,"%dS", Readlen-tempCiagr[5]);
	}
}
void RwriteCIAGR(int Readlen) //from tmepCiagr to outCiagr
{
	int outCiagrlen=0,temp, temp2;
	//temp=tempCiagr[4]-tempCiagr[2];
	if(tempCiagr[4]<0)
	{
		if((temp=Readlen-tempCiagr[2])<2)
		{
			//tempCiagr[2]=Readlen;
			if (tempCiagr[1] < 2)//c[1]>0  ????????should change the tempCiagr[0]---------
			{
				outCiagrlen = sprintf(outCiagr, "%dM", Readlen);
				mappedbase += Readlen;
			}
			else
			{
				outCiagrlen = sprintf(outCiagr, "%dS%dM", tempCiagr[1], Readlen - tempCiagr[1]);
				mappedbase += (Readlen - tempCiagr[1]);
			}
		}
		else
		{
			if (tempCiagr[1] < 2)//c[1]>0
			{
				outCiagrlen = sprintf(outCiagr, "%dM%dS", tempCiagr[2], temp);
				mappedbase += tempCiagr[2];
			}
			else
			{
				outCiagrlen = sprintf(outCiagr, "%dS%dM%dS", tempCiagr[1], tempCiagr[2] - tempCiagr[1], temp);
				mappedbase += (tempCiagr[2] - tempCiagr[1]);
			}
		}
	}
	else//c[4]>0 means there is two
	{
		if (tempCiagr[1] < 2)//c[1]>0
		{
			outCiagrlen = sprintf(outCiagr, "%dM", tempCiagr[2]);
			mappedbase += tempCiagr[2];
		}
		else
		{
			outCiagrlen = sprintf(outCiagr, "%dS%dM", tempCiagr[1], tempCiagr[2] - tempCiagr[1]);
			mappedbase += (tempCiagr[2] - tempCiagr[1]);
		}

		if( (temp=Readlen-tempCiagr[5])<2)
			tempCiagr[5]=Readlen;

		if(tempCiagr[3]==0)
		{
			outCiagrlen+=sprintf(outCiagr+outCiagrlen,"%dS%dM",tempCiagr[4]-tempCiagr[2],tempCiagr[5]-tempCiagr[4]);
			mappedbase += (tempCiagr[5] - tempCiagr[4]);
		}
		else if(tempCiagr[3]>0) //D
		{
			if ((temp = tempCiagr[4] - tempCiagr[2]) > 0)
			{
				outCiagrlen += sprintf(outCiagr + outCiagrlen, "%dD%dS%dM", tempCiagr[3], temp, tempCiagr[5] - tempCiagr[4]);
				mappedbase += (tempCiagr[5] - tempCiagr[4]);
			}
			else
			{
				outCiagrlen += sprintf(outCiagr + outCiagrlen, "%dD%dM", tempCiagr[3], tempCiagr[5] - tempCiagr[2]);
				mappedbase += (tempCiagr[5] - tempCiagr[2]);
			}
		}
		else //<0 : I
		{
			temp2=-tempCiagr[3];
			if ((temp = tempCiagr[4] - tempCiagr[2] - temp2) > 0)
			{
				outCiagrlen += sprintf(outCiagr + outCiagrlen, "%dI%dS%dM", temp2, temp, tempCiagr[5] - tempCiagr[4]);
				mappedbase += (tempCiagr[5] - tempCiagr[4]);
			}
			else
			{
				outCiagrlen += sprintf(outCiagr + outCiagrlen, "%dI%dM", temp2, tempCiagr[5] - tempCiagr[2] - temp2);
				mappedbase += (tempCiagr[5] - tempCiagr[2] - temp2);
			}
		}
		if( Readlen-tempCiagr[5]>1)
			outCiagrlen+=sprintf(outCiagr+outCiagrlen,"%dS", Readlen-tempCiagr[5]);
	}
}
//try fun: write the two CIAGR
void tryRwriteCIAGR(int s,int l)
{
	int outCiagrlen=0, temp, temp2;
	if(tempCiagr[4]<0)
	{
		if((temp=l-tempCiagr[2])<2)
		{
			//tempCiagr[2]=Readlen;
			outCiagrlen=sprintf(outCiagr,"%dM",l-s);
		}
		else
		{
			outCiagrlen=sprintf(outCiagr,"%dM%dS",tempCiagr[2]-s,temp );
		}
	}
	else//c[4]>0 means there is two
	{
		outCiagrlen=sprintf(outCiagr,"%dM",tempCiagr[2]-s );

		if( l-tempCiagr[5]<2)
			tempCiagr[5]=l;

		if(tempCiagr[3]==0)
		{
			outCiagrlen+=sprintf(outCiagr+outCiagrlen,"%dS%dM",tempCiagr[4]-tempCiagr[2],tempCiagr[5]-tempCiagr[4]);
		}
		else if(tempCiagr[3]>0) //D
		{
			if( (temp=tempCiagr[4]-tempCiagr[2])>0 )
				outCiagrlen+=sprintf(outCiagr+outCiagrlen,"%dD%dS%dM",tempCiagr[3],temp,tempCiagr[5]-tempCiagr[4]);
			else
				outCiagrlen+=sprintf(outCiagr+outCiagrlen,"%dD%dM",tempCiagr[3],tempCiagr[5]-tempCiagr[2]);
		}
		else //<0 : I
		{
			temp2=-tempCiagr[3];
			if( (temp=tempCiagr[4]-tempCiagr[2]-temp2)>0 )
				outCiagrlen+=sprintf(outCiagr+outCiagrlen,"%dI%dS%dM",temp2,temp,tempCiagr[5]-tempCiagr[4]);
			else
				outCiagrlen+=sprintf(outCiagr+outCiagrlen,"%dI%dM",temp2,tempCiagr[5]-tempCiagr[2]-temp2);
		}
		if( l-tempCiagr[5]>1)
			outCiagrlen+=sprintf(outCiagr+outCiagrlen,"%dS",l-tempCiagr[5]);
	}
}
void trywriteCIAGR(int s, int l)
{
	int outCiagrlen=0,temp, temp2;
	//temp=tempCiagr[4]-tempCiagr[2];
	if(tempCiagr[4]<0)
	{
		if((temp=l-tempCiagr[2])<2)
		{
			//outCiagrlen=sprintf(outCiagr,"%dM",l-s );
			if(tempCiagr[1]==s)
				outCiagrlen=sprintf(outCiagr,"%dM",l-s );
			else
				outCiagrlen=sprintf(outCiagr,"%dS%dM",tempCiagr[1]-s,l-tempCiagr[1]);
		}
		else
		{
			//outCiagrlen=sprintf(outCiagr,"%dM%dS",tempCiagr[2]-s,temp );
			if(tempCiagr[1]==s)//c[1]>0
				outCiagrlen=sprintf(outCiagr,"%dM%dS",tempCiagr[2]-s,temp );
			else
				outCiagrlen=sprintf(outCiagr,"%dS%dM%dS",tempCiagr[1]-s,tempCiagr[2]-tempCiagr[1],temp);
		}
	}
	else//c[4]>0 means there is two
	{
		//	outCiagrlen=sprintf(outCiagr,"%dM",tempCiagr[2]-s );
		if(tempCiagr[1]==s)
			outCiagrlen=sprintf(outCiagr,"%dM",tempCiagr[2]-s );
		else
			outCiagrlen=sprintf(outCiagr,"%dS%dM",tempCiagr[1]-s,tempCiagr[2]-tempCiagr[1]);

		if( l-tempCiagr[5]<2)
			tempCiagr[5]=l;

		if(tempCiagr[3]==0)
		{
			outCiagrlen+=sprintf(outCiagr+outCiagrlen,"%dS%dM",tempCiagr[4]-tempCiagr[2],tempCiagr[5]-tempCiagr[4]);
		}
		else if(tempCiagr[3]>0) //D
		{
			if( (temp=tempCiagr[4]-tempCiagr[2])>0 )
				outCiagrlen+=sprintf(outCiagr+outCiagrlen,"%dD%dS%dM",tempCiagr[3],temp,tempCiagr[5]-tempCiagr[4]);
			else
				outCiagrlen+=sprintf(outCiagr+outCiagrlen,"%dD%dM",tempCiagr[3],tempCiagr[5]-tempCiagr[2]);
		}
		else //<0 : I
		{
			temp2=-tempCiagr[3];
			if( (temp=tempCiagr[4]-tempCiagr[2]-temp2)>0 )
				outCiagrlen+=sprintf(outCiagr+outCiagrlen,"%dI%dS%dM",temp2,temp,tempCiagr[5]-tempCiagr[4]);
			else
				outCiagrlen+=sprintf(outCiagr+outCiagrlen,"%dI%dM",temp2,tempCiagr[5]-tempCiagr[2]-temp2);
		}
		if( l-tempCiagr[5]>1)
			outCiagrlen+=sprintf(outCiagr+outCiagrlen,"%dS",l-tempCiagr[5]);
	}
}
int  tradMatch( char *R, char *Rv) //tradition match
{
	long i,j, num, maxlen,len,flag;
	int word, offset;

	//init the getLenghout
	getLengthout[1]=0;
	getLengthout[2]=0;
	maxlen=VALUELEN;
	flag=4;
	if(Wsize>0) //normal
	{
		word=Word[0];
		num=*(IndexarrayNum+Word[0]);
		offset=ReadPreKeyarray[0];
		for(i=0; i<num; i++)
		{
			len=getLength( *(*(Indexarray+word)+i),offset,  R);

			if(maxlen<len)
			{
				tempCiagr[0]=-1;
				storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
				maxlen=len;
				flag=0;
			}
		}
	}
	if(maxlen==ReadLength) return flag;

	if(cp ==1){
		if( RWsize>0 )  //Reverse
		{
			word  =RWord[0];
			offset=RvP[0];
			num   =*(IndexarrayNum+word);
			for(i=0; i<num; i++)
			{
				len=getLength( *(*(Indexarray+word)+i),offset,Rv);
				if(maxlen<len)
				{
					tempCiagr[0]=-1;
					storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
					maxlen=len;
					flag  =16;
				}
			}
		}
	}
	return flag;
}

char *new_split_cigar(char *add_temp,char *seq,char *samcigar,char *flag,char *clip_seq)
{
    /*puts("Red54-0");
    puts(add_temp);
    puts(seq);
    puts(samcigar);
    puts(clip_seq);
    puts("Red54-1");*/
	size_t i=0,j=0,k=0;
	char NUM[30]={'\0'},unmap[]={"*"},char_cig[100]={'\0'};
	int num_cig[100]={0},num_cig_sum[100]={0};

	add_temp[0]='\0';
	clip_seq[0]='\0';

	if(strcmp(samcigar,unmap)==0){
		strcpy(add_temp,seq);
    //    printf("red54:new_split_cigar-0:add_temp:%s\n", add_temp);
		return add_temp;
	}
	else
	{
		/*
		for(i=0;i<strlen(samcigar);i++){
			num_cig[i]=0;
			num_cig_sum[i]=0;
		}*/
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
    //    printf("red54:new_split_cigar-1:num_cig_sum[%d]:%d\n", i, num_cig_sum[i]);
    //    printf("red54:new_split_cigar-1:num_cig[%d]:%d\n", i, num_cig[i]);
						 for(j=num_cig_sum[i];j<(num_cig[i]+num_cig_sum[i]);j++){
							 *(clip_seq+k)=*(seq+j);
                             clip_seq[k+1] = '\0';
      //  printf("red54:new_split_cigar-s:clip_seq:%s\n", clip_seq);
							 k++;
						 }
        //printf("red54:new_split_cigar-1:clip_seq:%s\n", clip_seq);
						 strcat(add_temp,clip_seq);
    //    printf("red54:new_split_cigar-1:add_temp:%s\n", add_temp);
						 clip_seq[0]='\0';
						 break;
					 }
				case 'I':{
      //  printf("red54:new_split_cigar-2:num_cig_sum[%d]:%d\n", i, num_cig_sum[i]);
        //printf("red54:new_split_cigar-2:num_cig[%d]:%d\n", i, num_cig[i]);
						 for(j=num_cig_sum[i];j<(num_cig[i]+num_cig_sum[i]);j++){
							 *(clip_seq+k)=*(seq+j);
                             clip_seq[k+1] = '\0';
							 k++;
						 }
    //    printf("red54:new_split_cigar-2:clip_seq:%s\n", clip_seq);
						 strcat(add_temp,clip_seq);
      //  printf("red54:new_split_cigar-2:add_temp:%s\n", add_temp);
						 clip_seq[0]='\0';
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
        //printf("red54:new_split_cigar-3:add_temp:%s\n", add_temp);
		return add_temp;
	}
}

char *get_mismatch(char *match_temp,const char *seq,const char *base,char *pos)
{
	match_temp[0]='\0';
	//printf("red54:match_temp1:%s\n",match_temp);
	//printf("red54:seq:%s\n",seq);
	//printf("red54:base:%s\n",base);
	//printf("red54:pos:%s\n",pos);

	if(pos[0]==0)
	{
		strcpy(match_temp,"\n");
		return match_temp;
	}
	size_t i=0,temp=0,first=0;
	char num[30]={'\0'};
	char match[READLEN]={'\0'};
	for(i=0;i<strlen(seq);i++)
	{
 //       printf("red54:*(seq+%d)=%c,*(base+%d)=%c\n", i, *(seq+i), i, *(base+i));
		if(*(seq+i)!=*(base+i))
		{
			temp=i-first;
			sprintf(num,"%zu",temp);
			first=i;
			strcat(match,num);
			match[strlen(match)]=*(seq+i);
          //  printf("red54: match: %s\n", match);
			memset(num,'\0',strlen(num));
		}
	}
	match[strlen(match)]='\n';
	strcpy(match_temp,match);
//	printf("red54:match_temp2:%s\n",match_temp);
	return match_temp;
}

char *new_get_base(char *base,char *add_,char *samcigar,char *ref_gen,char *pos,char *base_)
{
/*	printf("red54:base1:%s\n",base);
	printf("red54:add_:%s\n",add_);
	printf("red54:samcigar:%s\n",samcigar);
	//puts(ref_gen);
	printf("red54:pos:%s\n",pos);
	printf("red54:base_:%s\n",base_);*/
    memset(base,'\0',READLEN);
	base_[0]='\0';
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
//	printf("red54:base2:%s\n",base);
	return base;
}

int match(char *readStr, char *srcBuf, char *RStr, FILE*fp, char *ref_gen) //start(main) match--return mapped number!
{
	//printf("considerlen=%d\n", CONSIDERLEN);
	//getchar();
	int  i, j, PreKeyCount, result,temp, matchend,matchstart;
	unsigned int HashWord=0;
	char *R,mappos[10];
	//add_temp=m_malloc(READLEN);
	//clip_seq=m_malloc(READLEN);

	//memset(add_temp,'\0',sizeof(char)*READLEN);
	//memset(clip_seq,'\0',sizeof(char)*READLEN);

	i=0;  //or i=10
	PreKeyCount=0;
	while(  (*RStr=REV( *(readStr+i) ))>64 ) //find all PreKey
	{
		//( *(readStr+i)=Caps( *(readStr+i) ) )
		RStr--;
		if( *(readStr+i)==IndexPreKey[0] )
		{
			j=1;
			while( *(readStr+i+j)==IndexPreKey[j] ) j++;  //test j=end
			if( IndexPreKey[j]==0 )
			{
				ReadPreKeyarray[PreKeyCount++]=i+2;
			}
		}
		i++;
	}
	*(readStr+i)='\0';
	ReadLength=i;
	RStr++;
	if(PreKeyCount==0)  //no prekey
	{
		//fprintf(fp,"%d\t4\t0\t*\t%s\n",read_count,readStr);
		fprintf(cigar,"%c\n",'*');
		fprintf(pos,"%c\n",'0');
		//new_split_cigar(add_temp,readStr,"*","4",clip_seq);
		fprintf(add,"%s\n",readStr);
		//new_get_base(base_temp,readStr,"*",ref_gen,"0",base_);
		//get_mismatch(match_temp,readStr,base_temp,'\0');
		fputc('\n',cor);
		return mappedNum;
	}
	//get map position and return FLAG
	result=getMapPosi(PreKeyCount,readStr,RStr);  //normal and reverse
	//if(outFLAG) //0 means readStr; 1 means RStr
	//	R=RStr;
	//else R=readStr;
	//write the detail to file  normal:

	if( result==2)  //match two ...
	{
		i=*out[0]; j=*Rout[0]; temp=ReadLength/2;  //temp is half of readlength
		if(out[1][i]+HashStrlen<temp && Rout[1][j]+HashStrlen<temp)  //match to start:normal & reverse
		{
			//printf("result2_one\n");
			getLength(out[0][1], out[1][1], readStr);//out[0] is src posi; out[1] is read offset
			tempCiagr[0]=-1;
			//------------
			//printf("nor 0: %d;  1:%d;  2:%d\nout1i:%d\n",getLengthout[0],getLengthout[1],getLengthout[2],out[1][i]);
			//------------
			storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
			if(tempCiagr[2]< out[1][ i ])  //modifier the func:writeOutArray :sort o[]
			{
				getLength(out[0][ i ],out[1][ i ],readStr);
				if(getLengthout[1]<tempCiagr[1])  //stop negative I or D
					tempCiagr[0]=-1;

				storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
			}
			//printf("ReadLength-Rout[1][1]<getLengthout[2]:%d %d %d\n",Rout[1][1],Rout[1][j],getLengthout[2]);
			if(ReadLength-Rout[1][j]<getLengthout[2])  //correct...5-16
			{
				writeCIAGR(ReadLength);
				//fprintf(fp,"%d\t0\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,readStr);
				fprintf(cigar,"%s\n",outCiagr);
				fprintf(pos,"%d\n",tempCiagr[0]+1);
				new_split_cigar(add_temp,readStr,outCiagr,"0",clip_seq);
				fprintf(add,"%s\n",add_temp);
				sprintf(mappos, "%d", tempCiagr[0]+1);
				new_get_base(base_temp,add_temp,outCiagr,ref_gen,mappos,base_);
				get_mismatch(match_temp,readStr,base_temp,mappos);
				fprintf(cor,"%s",match_temp);
				//
				mappedNum++;
				return mappedNum;
			}
			else
				matchend=getLengthout[2];

			writeCIAGR(matchend); //part of read
			readStr[matchend]='\0';
			//fprintf(fp,"%d\t0\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,readStr);
			//
			fprintf(cigar,"%s\n",outCiagr);
			fprintf(pos,"%d\n",tempCiagr[0]+1);
			new_split_cigar(add_temp,readStr,outCiagr,"0",clip_seq);
			fprintf(add,"%s\n",add_temp);
			sprintf(mappos, "%d", tempCiagr[0]+1);
			new_get_base(base_temp,add_temp,outCiagr,ref_gen,mappos,base_);
			get_mismatch(match_temp,readStr,base_temp,mappos);
			fprintf(cor,"%s",match_temp);
			//fputs(outCiagr,cigar);
			//fputc('\n',cigar);
			//-------------reverse  part-----------------
			getLength(Rout[0][j],Rout[1][j],RStr);//out[0] is src posi; out[1] is read offset
			tempCiagr[0]=-1;
			//---------test---------
			//printf("rev 0: %d;  1:%d;  2:%d\n",getLengthout[0],getLengthout[1],getLengthout[2]);
			//--------------------
			//			if(getLengthout[2]< Rout[1][ j ])  //modifier the func:writeOutArray :sort o[]
			//			{
			//				storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
			//				getLength(Rout[0][ j ],Rout[1][ j ],RStr);
			//
			//				if(getLengthout[1]<tempCiagr[1])  //stop negative I or D
			//					tempCiagr[0]=-1;
			//			}
			if(getLengthout[2]>ReadLength-matchend) //overlap
			{
				getLengthout[2]=ReadLength-matchend;
				//if(getLengthout[2]<getLengthout[1])
				//tempCiagr[0]=-1;
			}
			storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);

			matchend=ReadLength-matchend;
			RwriteCIAGR(matchend);
			RStr[matchend]='\0';
			//fprintf(fp,"%d\t16\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,RStr);
			fprintf(cigar,"0%s\n",outCiagr);
			fprintf(pos,"%d\n",tempCiagr[0]+1);
			new_split_cigar(add_temp,RStr,outCiagr,"16",clip_seq);
			fprintf(add,"%s\n",add_temp);
			sprintf(mappos, "%d", tempCiagr[0]+1);
			new_get_base(base_temp,add_temp,outCiagr,ref_gen,mappos,base_);
		//	printf("match:%s\nstr:%s\nbase:%s\nmappos:%s\n",match_temp,RStr,base_temp,mappos);
			get_mismatch(match_temp,RStr,base_temp,mappos);
			fprintf(cor,"%s",match_temp);

			mappedNum++; doubleMap++;
		}
		else if(out[1][1]>temp && Rout[1][1]>temp)  // match to end:reverse & normal
		{
			//printf("result2_two\n");
			//-------------reverse  part-----------------
			getLength(Rout[0][1],Rout[1][1],RStr);//out[0] is src posi; out[1] is read offset
			tempCiagr[0]=-1;
			//--------
			//printf("out11: %d;  Rout11:%d;  temp:%d\n",out[1][1],Rout[1][1],temp);
			//printf("out1j: %d;  Rout1j:%d;  i j:%d %d\n",out[1][i],Rout[1][j],i,j);
			//---------
			storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
			if(tempCiagr[2]< Rout[1][ j ])  //modifier the func:writeOutArray :sort o[]
			{
				getLength(Rout[0][ j ],Rout[1][ j ],RStr);
				if(getLengthout[1]<tempCiagr[1])  //stop negative I or D
					tempCiagr[0]=-1;

				storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
			}
			if(out[1][i]<ReadLength-tempCiagr[1])  //correct...5-16
			{
				RwriteCIAGR(ReadLength);
				//fprintf(fp,"%d\t16\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,RStr);
				fprintf(cigar,"0%s\n",outCiagr);
				fprintf(pos,"%d\n",tempCiagr[0]+1);
				new_split_cigar(add_temp,RStr,outCiagr,"16",clip_seq);
				fprintf(add,"%s\n",add_temp);
				sprintf(mappos, "%d", tempCiagr[0]+1);
				new_get_base(base_temp,add_temp,outCiagr,ref_gen,mappos,base_);
			//	printf("match:%s\nstr:%s\nbase:%s\nmappos:%s\n",match_temp,RStr,base_temp,mappos);
				get_mismatch(match_temp,RStr,base_temp,mappos);
				fprintf(cor,"%s",match_temp);
				mappedNum++;
				return mappedNum;
			}
			else
				matchstart=tempCiagr[1];

			tryRwriteCIAGR(matchstart,ReadLength);
			R=RStr+matchstart;
			//fprintf(fp,"%d\t16\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,R);
			fprintf(cigar,"0%s\n",outCiagr);
			fprintf(pos,"%d\n",tempCiagr[0]+1);
			new_split_cigar(add_temp,R,outCiagr,"16",clip_seq);
			fprintf(add,"%s\n",add_temp);
			sprintf(mappos, "%d", tempCiagr[0]+1);
			new_get_base(base_temp,add_temp,outCiagr,ref_gen,mappos,base_);
		//	printf("match:%s\nstr:%s\nbase:%s\nmappos:%s\n",match_temp,R,base_temp,mappos);
			get_mismatch(match_temp,R,base_temp,mappos);
			fprintf(cor,"%s",match_temp);

			//---------normal match----------------
			getLength(out[0][i],out[1][i],readStr);//out[0] is src posi; out[1] is read offset
			tempCiagr[0]=-1;
			//---------test---------
			//printf("nor 0: %d;  1:%d;  2:%d\n",getLengthout[0],getLengthout[1],getLengthout[2]);
			//--------------------
			if(getLengthout[1]<ReadLength-matchstart) //overlap
			{	getLengthout[0]+=ReadLength-matchstart-getLengthout[1]; getLengthout[1]=ReadLength-matchstart;  }
			//			if(getLengthout[1]<ReadLength-matchstart) //overlap
			//			{
			//				getLengthout[0]+=ReadLength-matchstart-getLengthout[1]; getLengthout[1]=ReadLength-matchstart;
			//
			//			}
			//
			//			if(getLengthout[2]< out[1][ i ])  //modifier the func:writeOutArray :sort o[]
			//			{
			//				if(getLengthout[2]>getLengthout[1])
			//					storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
			//				getLength(out[0][ i ],out[1][ i ],readStr);
			//				if(getLengthout[1]<tempCiagr[1])  //stop negative I or D
			//				{
			//				   tempCiagr[0]=-1;
			//				   if(getLengthout[1]<ReadLength-matchstart) //overlap
			//			     {	getLengthout[0]+=ReadLength-matchstart-getLengthout[1]; getLengthout[1]=ReadLength-matchstart;  }
			//		    }
			//
			//			}

			storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);

			matchstart=ReadLength-matchstart;
			trywriteCIAGR(matchstart,ReadLength);
			R=readStr+matchstart;
			//fprintf(fp,"%d\t0\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,R);
			fprintf(cigar,"%s\n",outCiagr);
			fprintf(pos,"%d\n",tempCiagr[0]+1);
			new_split_cigar(add_temp,R,outCiagr,"0",clip_seq);
			fprintf(add,"%s\n",add_temp);
			sprintf(mappos, "%d", tempCiagr[0]+1);
			new_get_base(base_temp,add_temp,outCiagr,ref_gen,mappos,base_);
			get_mismatch(match_temp,R,base_temp,mappos);
			fprintf(cor,"%s",match_temp);
			mappedNum++;  doubleMap++;
		}
		else
		{
			//printf("result2_three\n");
			if(i>j)  //compare the more matched kmer
			{
				temp=*out[0];  //try Normal
				getLength(out[0][1],out[1][1],readStr);//out[0] is src posi; out[1] is read offset
				tempCiagr[0]=-1;
				storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
				if(tempCiagr[2]< out[1][ temp ])  //modifier the func:writeOutArray :sort o[]
				{
					getLength(out[0][ temp ],out[1][ temp ],readStr);
					if(getLengthout[1]<tempCiagr[1])  //stop negative I or D
						tempCiagr[0]=-1;
					storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
				}
				//if(outFLAG==0)
				writeCIAGR(ReadLength);
				//else
				//RwriteCIAGR(ReadLength);

				//fprintf(fp,"%d\t0\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,readStr);
				fprintf(cigar,"%s\n",outCiagr);
				fprintf(pos,"%d\n",tempCiagr[0]+1);
				new_split_cigar(add_temp,readStr,outCiagr,"0",clip_seq);
				fprintf(add,"%s\n",add_temp);
				sprintf(mappos, "%d", tempCiagr[0]+1);
				new_get_base(base_temp,add_temp,outCiagr,ref_gen,mappos,base_);
				get_mismatch(match_temp,readStr,base_temp,mappos);
				fprintf(cor,"%s",match_temp);
				mappedNum++;
			}
			else
			{
				temp=*Rout[0];  //try reverse
				getLength(Rout[0][1],Rout[1][1],RStr);//out[0] is src posi; out[1] is read offset
				tempCiagr[0]=-1;
				storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
				if(tempCiagr[2]< Rout[1][ temp ])  //modifier the func:writeOutArray :sort o[]
				{
					getLength(Rout[0][ temp ],Rout[1][ temp ],RStr);
					if(getLengthout[1]<tempCiagr[1])  //stop negative I or D
						tempCiagr[0]=-1;
					storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
				}
				RwriteCIAGR(ReadLength);

				//fprintf(fp,"%d\t16\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,RStr);
				fprintf(cigar,"0%s\n",outCiagr);
				fprintf(pos,"%d\n",tempCiagr[0]+1);
				new_split_cigar(add_temp,RStr,outCiagr,"16",clip_seq);
				fprintf(add,"%s\n",add_temp);
				sprintf(mappos, "%d", tempCiagr[0]+1);
				new_get_base(base_temp,add_temp,outCiagr,ref_gen,mappos,base_);
			//	printf("match:%s\nstr:%s\nbase:%s\nmappos:%s\n",match_temp,RStr,base_temp,mappos);
				get_mismatch(match_temp,RStr,base_temp,mappos);
				fprintf(cor,"%s",match_temp);
				mappedNum++;
			}
		}
	}
	else if( result==0 )  //normal match
	{
		////-------contain 2 or 3------
		//getLength(out[0][1],out[1][1],ShortReadLength,R);//out[0] is src posi; out[1] is read offset
		//tempCiagr[0]=-1;
		//storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
		//if(tempCiagr[1]>out[1][2] || tempCiagr[2]<out[1][2])
		//{
		//	getLength(out[0][2],out[1][2],ShortReadLength,R);
		//	storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
		//}
		//
		//if(*out[0]>3)
		//{
		//	if(tempCiagr[4]>0) //two match
		//	{
		//		if(out[1][3]<tempCiagr[1] || out[1][3]>tempCiagr[5] || (out[1][3]>tempCiagr[2] && out[1][3]<tempCiagr[4]))
		//		{
		//			getLength(out[0][3],out[1][3],ShortReadLength,R);
		//			storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
		//		}
		//	}
		//	else //one match
		//		if(tempCiagr[1]>out[1][3] || tempCiagr[2]<out[1][3])
		//		{
		//			getLength(out[0][3],out[1][3],ShortReadLength,R);
		//			storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
		//		}
		//}
		////----consider more ---------//
		temp=*out[0];
		getLength(out[0][1],out[1][1],readStr);//out[0] is src posi; out[1] is read offset
		tempCiagr[0]=-1;
		storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
		if(tempCiagr[2]< out[1][ temp ])  //modifier the func:writeOutArray :sort o[]
		{
			getLength(out[0][ temp ],out[1][ temp ],readStr);
			if(getLengthout[1]<tempCiagr[1])  //stop negative I or D
				tempCiagr[0]=-1;
			storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
		}
		//if(outFLAG==0)
		writeCIAGR(ReadLength);
		//else
		//RwriteCIAGR(ReadLength);

		//fprintf(fp,"%d\t0\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,readStr);
		fprintf(cigar,"%s\n",outCiagr);
		fprintf(pos,"%d\n",tempCiagr[0]+1);
		new_split_cigar(add_temp,readStr,outCiagr,"0",clip_seq);
		fprintf(add,"%s\n",add_temp);
		sprintf(mappos, "%d", tempCiagr[0]+1);
		new_get_base(base_temp,add_temp,outCiagr,ref_gen,mappos,base_);
		get_mismatch(match_temp,readStr,base_temp,mappos);
		fprintf(cor,"%s",match_temp);
		mappedNum++;
	}
	else if(result==1) //reverse match
	{
		temp=*Rout[0];
		i=getLength(Rout[0][1],Rout[1][1],RStr);//out[0] is src posi; out[1] is read offset
		tempCiagr[0]=-1;
		storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
		//printf("match:%d;temp:%d  0 1 2:%d  %d  %d\n",i,temp,getLengthout[0],getLengthout[1],getLengthout[2]);
		if(tempCiagr[2]< Rout[1][ temp ])  //modifier the func:writeOutArray :sort o[]
		{
			getLength(Rout[0][ temp ],Rout[1][ temp ],RStr);
			if(getLengthout[1]<tempCiagr[1])  //stop negative I or D
				tempCiagr[0]=-1;
			storetempCiagr(getLengthout[0],getLengthout[1],getLengthout[2]);
		}
		RwriteCIAGR(ReadLength);

		//fprintf(fp,"%d\t16\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,RStr);
		fprintf(cigar,"0%s\n",outCiagr);
		fprintf(pos,"%d\n",tempCiagr[0]+1);
		new_split_cigar(add_temp,RStr,outCiagr,"16",clip_seq);
		fprintf(add,"%s\n",add_temp);
		sprintf(mappos, "%d", tempCiagr[0]+1);
		new_get_base(base_temp,add_temp,outCiagr,ref_gen,mappos,base_);
	//	printf("match:%s\nstr:%s\nbase:%s\nmappos:%s\n",match_temp,RStr,base_temp,mappos);
		get_mismatch(match_temp,RStr,base_temp,mappos);
		fprintf(cor,"%s",match_temp);
		mappedNum++;
	}
	else //-1  or -2
	{
		outFLAG=tradMatch(readStr,RStr);
		if(outFLAG==4){ //value match
			//fprintf(fp,"%d\t4\t0\t*\t%s\n",read_count,readStr);
			fprintf(cigar,"%c\n",'*');
			fprintf(pos,"%c\n",'0');
			//new_split_cigar(add_temp,readStr,"*","4",clip_seq);
			fprintf(add,"%s\n",readStr);
			//new_get_base(base_temp,readStr,"*",ref_gen,"0",base_);
			//get_mismatch(match_temp,readStr,base_temp,'\0');
			//fprintf(cor,"%s",match_temp);
			fputc('\n',cor);
		}
		else
		{
			mappedNum++;
			if(outFLAG==0)
			{
				writeCIAGR(ReadLength);
				//fprintf(fp,"%d\t%d\t%d\t%s\t%s\n",read_count,outFLAG,tempCiagr[0]+1,outCiagr,readStr);
				fprintf(cigar,"%s\n",outCiagr);
				fprintf(pos,"%d\n",tempCiagr[0]+1);
				new_split_cigar(add_temp,readStr,outCiagr,"0",clip_seq);
				fprintf(add,"%s\n",add_temp);
				sprintf(mappos, "%d", tempCiagr[0]+1);
				new_get_base(base_temp,add_temp,outCiagr,ref_gen,mappos,base_);
				get_mismatch(match_temp,readStr,base_temp,mappos);
				fprintf(cor,"%s",match_temp);
			}
			else if(cp ==1)//outFLAG==16 reverse
			{
				RwriteCIAGR(ReadLength);
				//fprintf(fp,"%d\t%d\t%d\t%s\t%s\n",read_count,outFLAG,tempCiagr[0]+1,outCiagr,RStr);
				fprintf(cigar,"0%s\n",outCiagr);
				fprintf(pos,"%d\n",tempCiagr[0]+1);
				new_split_cigar(add_temp,RStr,outCiagr,"16",clip_seq);
				fprintf(add,"%s\n",add_temp);
				sprintf(mappos, "%d", tempCiagr[0]+1);
				new_get_base(base_temp,add_temp,outCiagr,ref_gen,mappos,base_);
			//	printf("match:%s\nstr:%s\nbase:%s\nmappos:%s\n",match_temp,RStr,base_temp,mappos);
				get_mismatch(match_temp,RStr,base_temp,mappos);
				fprintf(cor,"%s",match_temp);
			}
		}
	}
	//free(add_temp); add_temp=NULL;
	//free(clip_seq); clip_seq=NULL;
	return mappedNum;
}

void get_qscore(char *file,char *qsread, FILE *fastq)
{
	/*size_t buffer_size=get_file_size(file),buf_count=0;
	if(buffer_size<block_size)
		buf_count = buffer_size;
	else buf_count = block_size;*/

	//char *_buffer=m_malloc(buf_count+READLEN);
	//char *read_temp=m_malloc(READLEN);
	//char *meta_temp=m_malloc(READLEN);
	char *qs_temp=m_malloc(READLEN);
	//char *meta_mid=m_malloc(READLEN);
	//char *meta_part=m_malloc(READLEN);

	//memset(_buffer,'\0',sizeof(char)*(buf_count+READLEN));
	//memset(read_temp,'\0',sizeof(char)*READLEN);
	//memset(meta_temp,'\0',sizeof(char)*READLEN);
	memset(qs_temp,'\0',sizeof(char)*READLEN);
	/*memset(meta_mid,'\0',sizeof(char)*READLEN);
	memset(meta_part,'\0',sizeof(char)*READLEN);

	size_t fread_num=0,read_num=1,i=0,j=0,k=0,l=0,codelen=0;
	char num[20]={'\0'};
	//Copy the first metadata
	strcpy(_buffer,metaread);
	//get_meta(meta_mid,_buffer);
	sscanf(_buffer,"%*s%s",meta_mid);
	char *ptr=NULL;
	ptr=strstr(meta_mid,"length=");

	//memset(meta_temp,'\0',strlen(meta_temp));
	//memset(meta_mid,'\0',strlen(meta_mid));

	//Incremental coding of metadata
	//get_meta(meta_temp,read_temp);
    long pos, now;
    pos = ftell(fastq);
    fgets(metaread,READLEN,fastq);
    now = ftell(fastq);
    if (pos == now) goto qs;
    fseek(fastq, pos, SEEK_SET);
	sscanf(metaread,"%*s%s",meta_temp);
	if(ptr==NULL){
        //printf("Red54: meta_temp=%s\n", meta_temp);
        //printf("Red54: meta_mid=%s\n", meta_mid);
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
	//codelen is the coded length of qualityScore
	//codelen=rll_quality_score(qs_temp, read_temp);
	//codelen=fwrite(qs_temp,1,codelen,qs);
	//get_qs(qs_temp,read_temp);
qs:*/
	sscanf(qsread,"%s",qs_temp);
	fwrite(qs_temp,1,strlen(qs_temp),qs);
	fputc('\n',qs);

	//memset(read_temp,'\0',strlen(read_temp));

	//fprintf(stderr,"read_num=%d\n",read_num/4);
	//free(_buffer);  _buffer=NULL;
	//free(meta_temp); meta_temp=NULL;
	//free(meta_mid);  meta_mid=NULL;
	//free(meta_part);  meta_part=NULL;
	free(qs_temp);  qs_temp=NULL;
	//free(read_temp);  read_temp=NULL;
}

int map(char *readFileName, char *buf, char *outFileName, char *argv[])  //main call-func
{
	FILE *readFile, *outFile;
	int  readStrSize;
	char readStr[READLEN]={0}, tempStr[READLEN]={0}, metaStr[READLEN]={0},qsStr[READLEN]={0},RevStr[READLEN]={0};
	char *p=&RevStr[80000];
	long long fread_num=0,buffer_size=0,buf_count=0,j=0,i=0,l=0,k=0;
	//puts("0");
	buffer_size=get_file_size(argv[1]);
	//printf("fasta file size:%ld\n",buffer_size);
	buf_count = buffer_size+10000;
	//printf("buf_count:%ld\n",buf_count);
	char *ref_temp=m_malloc(buf_count);
	char *ref_gen=m_malloc(buf_count);

	memset(ref_temp,'\0',sizeof(char)*buf_count);
	memset(ref_gen,'\0',sizeof(char)*buf_count);
	//puts("1");
	char ref_ch;

	long long buffer_sizes=get_file_size(argv[2]),buf_counts=0,read_num=0,fread_nums=0;
	if(buffer_sizes<block_size)
		buf_counts = buffer_sizes;
	else buf_counts = block_size;
	char *_buffer=m_malloc(buf_counts+READLEN);
	char *read_temp=m_malloc(READLEN);
	//char *meta_temp=m_malloc(READLEN);
	memset(_buffer,'\0',sizeof(char)*(buf_counts+READLEN));
	memset(read_temp,'\0',sizeof(char)*READLEN);
	//memset(meta_temp,'\0',sizeof(char)*READLEN);


	ref_g=fopen(argv[1],"rb");
	fgets(ref_temp,2000,ref_g);
	j=0;
	while((fread_num=fread(ref_temp,sizeof(char),buf_count, ref_g))){
		//printf("fread_num=%ld\n",fread_num);
		for(i=0;i<fread_num;i++){
			if( ( ref_ch=m_toupper(*(ref_temp+i)) )>64){
				*(ref_gen+j)=ref_ch;
				//printf("%c",ref_ch);
				j++;
			}
		}
	}
	//puts("2");
	free(ref_temp);ref_temp=NULL;
	fclose(ref_g);

	outFile=fopen(outFileName,"wb");
	readFile=fopen(readFileName, "rb");
	if( !readFile || !outFile )
	{
		printf("error readfile or outFile !\n");
		return -1;
	}
	ReadBuf=(char *)malloc(RBUFSIZE);
	if(!ReadBuf) return -2;
	fprintf(outFile,"@SQ\tSN:%s\tLN:%d\n",title,BaseCount);//write infomation for outFile
	fprintf(outFile,"@PG\tID:map\tPN:map\tVN:0.1\tCL:%s %s %s\n",argv[0],argv[1],argv[2]);

	read_count=0;
	mappedNum =0;
	doubleMap =0;
	//fread(ReadBuf, sizeof(char), RBUFSIZE, readFile);  //block read:depend on @ & +

	strcpy(pos_name,argv[2]);
	strcat(pos_name,".pos.txt");
	strcpy(cigar_name,argv[2]);
	strcat(cigar_name,".cigar.txt");
	strcpy(cor_name,argv[2]);
	strcat(cor_name,".cor.txt");
	strcpy(add_name,argv[2]);
	strcat(add_name,".add.txt");
	//strcpy(meta_name,argv[2]);
	//m_strcat(meta_name,".meta.txt");
	strcpy(qs_name,argv[2]);
	m_strcat(qs_name,".qs.txt");

	pos=fopen(pos_name,"wb");
	cor=fopen(cor_name,"wb");
	cigar=fopen(cigar_name,"wb");
	add=fopen(add_name,"wb");
	//meta=fopen(meta_name,"wb");
	qs=fopen(qs_name,"wb");
	//puts("3");
	/*fgets(metaStr,READLEN,readFile);
	strcpy(_buffer,metaStr);
	//get_meta(meta_mid,_buffer);
	//sscanf(_buffer,"%*s%s",meta_mid);
	for(i=0,l=0;i<strlen(_buffer);i++,l++)
	{
		if(_buffer[i]==' '&&_buffer[i+1]=='l'&&_buffer[i+2]=='e')//remove "length="
			break;
		else *(meta_temp+l)=*(_buffer+i);
	}
	fputs(meta_temp,meta);
	fputc('\n',meta);
	fseek(readFile,0,SEEK_SET);*/
	//frequency IO read
	//puts("4");
	while((fread_nums=fread(_buffer,sizeof(char),buf_counts,readFile)))
	{	//puts("5");
		_buffer[fread_nums]='\0';
		//puts("6");
		for(i=0;i<fread_nums;i++){
			*(read_temp+k)=*(_buffer+i);
			k++;
			if(*(_buffer+i)=='\n')
			{	
				read_num++;
				if(read_num%4==2){
					long lens = 0;
					lens = strlen(read_temp) - 1;
					CONSIDERLEN = lens *E;
					//printf("considerlen=%d,length=%d    \n", CONSIDERLEN, lens); getchar();
					totalbase += lens;
					//if(read_temp) //there is a short rea
					//{
						read_count++;
					//printf("read_temp=%s\n",read_temp);
						strcpy(readStr,read_temp);
						match(readStr, buf, p, outFile,ref_gen); //p as the rev
					//}
				}
				//puts("6-4");
				//fgets(tempStr,READLEN,readFile);  //quality data
				//fgets(qsStr,READLEN,readFile);
				if(read_num%4==0){
					strcpy(qsStr,read_temp);
					get_qscore(argv[2],qsStr,readFile);
				}
				//puts("6-5");
				memset(read_temp,'\0',strlen(read_temp));
				k=0;
			}
		}
		//puts("7");
	}
	
	/*printf("Unmapped read= %d\n",read_count-mappedNum);
	  printf("All Mapped reads= %d\n",mappedNum);
	  printf("Exact Match= %d\n",exactMatch);
	  printf("Inexact Match= %d\n",mappedNum-exactMatch);
	  printf("N*= %d\n\n",doubleMap);
	  printf("Total readCounts: %d\nMappedread rate: %f\n",read_count,(float)mappedNum/read_count);
	  printf("Total base:%d\nMappedbase:%d\nMappedbase rate:%f\n", totalbase, mappedbase,(float)mappedbase / totalbase);*/
	FILE *Out=NULL;
	Out=fopen("Output.txt","a");
	fprintf(Out,"%lld,%lld,%lld,%lld,%lld,%lld\n",read_count,mappedNum,exactMatch,doubleMap,totalbase,mappedbase);
	fclose(Out);
	fclose(readFile);
	fclose(outFile);
	fclose(pos);
	fclose(cigar);
	fclose(add);
	fclose(cor);
	fclose(qs);
	//fclose(meta);
	free(ref_gen); ref_gen=NULL;
	free(_buffer);  _buffer=NULL;
	free(read_temp);  read_temp=NULL;
	return 0;
}


//normal forword
int NorFor(int Mposi, char *pstr, int p)
{
	return 0;
}
//normal backword
int NorBack()
{
	return 0;
}
//reverse forword
int RevFor()
{
	return 0;
}
//reverse backword
int RevBac()
{
	return 0;
}

int ErrNorFor(char *Read, int a, char *Gen, int b) //CONSIDERLEN=10 to ERR rate: 20%
{
	int i=1, unmatchBase=0, ErrRate=CONSIDERLEN/5;
	while(i<=CONSIDERLEN)  //no consider the Insert and Delete
	{
		if( *(Read+a+i)!=*(Gen+b+i) )
			unmatchBase++;
		i++;
	}
	if(unmatchBase<ErrRate) return 1;
	else return 0;
}
int ErrNorBack(char *Read, int a,char *Gen,int b)  //CONSIDERLEN=10 to ERR rate: 20%
{
	return ErrNorFor(Read,a-11,Gen,b-11);  //offset = Consider + 1
}
