#include "head.h"

//there are some assistant function

int base_hash(unsigned char a, int b) 
{
	/*
convert:
W M R N : A
S Y : C
K : T	  
other : A
	 */
	int base_status=0;
	switch(a)
	{
		case 'A':
		case 'M':
		case 'R':base_status=1;
		case 'W':a=0; break;
		case 'C':
		case 'Y':base_status=1;
		case 'S':a=1; break;
		case 'G':a=2; break;
		case 'T':base_status=1;
		case 'K':a=3; break;
		case 'N':base_status=1; a=0; break;    
		default :base_status=2; a=0; 
	}
	b=b<<2;
	b=b+a;
	return b;
}
char REV(char a)
{
	switch(a)
	{
		case 'A':
			//	case 'M':
			//	case 'R':
			//	case 'W':
			return 'T'; 
		case 'C':
			//case 'Y':
			//case 'S':
			return 'G'; 
		case 'G':
			return 'C'; 
		case 'T':
			//case 'K':
			return 'A'; 

		default : return a;
	}

}
int REVbase_hash(unsigned char a, int b)  //127-base_hash
{
	/*
convert:
W M R N : A
S Y : C
K : T	  
other : A
	 */
	int base_status=0;
	switch(a)
	{
		case 'A':
		case 'M':
		case 'R':base_status=1;
		case 'W':a=3; break;
		case 'C':
		case 'Y':base_status=1;
		case 'S':a=2; break;
		case 'G':a=1; break;
		case 'T':base_status=1;
		case 'K':a=0; break;
		case 'N':base_status=1; a=3; break;    
		default : //printf(" the base is error ---%c--%d---\n",a,a);  //test accident
			 base_status=2; a=3; 
	}
	b=b<<2;
	b=b+a;
	return b;
}
int stringcmp(char *str1, char *str2) //strlen(str2) based on str2
{
	int status=0, i;
	char c;
	for(i=0; ( c= *(str2+i) )!='\0'; i++ )
	{		
		if( *(str1+i)!=c) 
		{
			status=-1;
			break;
		}
	}
	return status; //0 means equal
}
char Caps(char a) //get Caps 
{
	if (a>96)	return a-32;
	else return a;
}
int _pow(int a, int b)
{
	int i, result=1;
	for(i=0; i<b; i++)
		result=a*result;
	return result;
}
