#include "head.h"

//map the short read to the buf_sequence
#define READLEN     300000  //1024*80
#define RBUFSIZE    104857600  //100M
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
int match(char *readStr, char *srcBuf, char *RStr, FILE*fp) //start(main) match--return mapped number!
{
	//printf("considerlen=%d\n", CONSIDERLEN);
	//getchar();
	int  i, j, PreKeyCount, result,temp, matchend,matchstart;
	unsigned int HashWord=0;
	char *R;
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
		fprintf(fp,"%d\t4\t0\t*\t%s\n",read_count,readStr);
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
				fprintf(fp,"%d\t0\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,readStr);
				mappedNum++;
				return mappedNum;
			}
			else
				matchend=getLengthout[2];

			writeCIAGR(matchend); //part of read
			readStr[matchend]='\0';
			fprintf(fp,"%d\t0\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,readStr);

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
			fprintf(fp,"%d\t16\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,RStr);
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
				fprintf(fp,"%d\t16\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,RStr);
				mappedNum++;
				return mappedNum;
			}
			else
				matchstart=tempCiagr[1];

			tryRwriteCIAGR(matchstart,ReadLength);
			R=RStr+matchstart;
			fprintf(fp,"%d\t16\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,R);

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
			fprintf(fp,"%d\t0\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,R);
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

				fprintf(fp,"%d\t0\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,readStr);
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

				fprintf(fp,"%d\t16\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,RStr);
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

		fprintf(fp,"%d\t0\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,readStr);
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

		fprintf(fp,"%d\t16\t%d\t%s\t%s\n",read_count,tempCiagr[0]+1,outCiagr,RStr);
		mappedNum++;
	}
	else //-1  or -2
	{
		outFLAG=tradMatch(readStr,RStr);		
		if(outFLAG==4) //value match
			fprintf(fp,"%d\t4\t0\t*\t%s\n",read_count,readStr);		
		else 
		{
			mappedNum++;						
			if(outFLAG==0)
			{
				writeCIAGR(ReadLength);
				fprintf(fp,"%d\t%d\t%d\t%s\t%s\n",read_count,outFLAG,tempCiagr[0]+1,outCiagr,readStr);
			}
			else if(cp ==1)//outFLAG==16 reverse
			{
				RwriteCIAGR(ReadLength);
				fprintf(fp,"%d\t%d\t%d\t%s\t%s\n",read_count,outFLAG,tempCiagr[0]+1,outCiagr,RStr);
			}
		}
	}
	return mappedNum;
}
int map(char *readFileName, char *buf, char *outFileName, char *argv[])  //main call-func
{
	FILE *readFile, *outFile;
	int  readStrSize;
	char readStr[READLEN]={0}, tempStr[READLEN]={0}, RevStr[READLEN]={0};
	char *p=&RevStr[80000];
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

	//frequency IO read 
	while( fgets(tempStr,READLEN,readFile)!=NULL )
	{
		int lens = 0;
		fgets(readStr, READLEN, readFile);
		lens = strlen(readStr) - 1;
		CONSIDERLEN = lens *E;
		//printf("considerlen=%d,length=%d    \n", CONSIDERLEN, lens); getchar();
		totalbase += lens;

		if(readStr) //there is a short rea
		{
			read_count++;		
			match(readStr, buf, p, outFile); //p as the rev
		}
		fgets(tempStr,READLEN,readFile);  //quality data
		fgets(tempStr,READLEN,readFile);
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
	fclose(readFile);  fclose(outFile);
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
