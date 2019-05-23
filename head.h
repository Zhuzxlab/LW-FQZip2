#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/stat.h>
//#include <math.h>

//linux head file 
//#include <dirent.h>
//#include <sys/stat.h>
//#include <unistd.h>
//#include <errno.h> 


//public para

extern char *title, *genBaseBuf;        //src file info

extern int **Indexarray, *IndexarrayNum, HashStrlen, IndexPosiCount,  BaseCount, *IndexPosition,VALUELEN;//index info 
extern char IndexPreKey[];
extern float E;
extern int *ref_index,  refunit, Readlen,CONSIDERLEN,cp;  

//associate function
int base_hash(unsigned char a, int b);
int stringcmp(char *str1, char *str2);
char Caps(char a);
char REV(char a);
//genProcess
char * genProcess(char *genFileName, char *genBaseBuf, char *PreKey);


//map function
int map(char *readFileName, char*buf, char *outFileName, char *argv[]);


//int output_test( chp filename, chp base_title );  //test
