#define NDEBUG

#include <algorithm>
#include <cassert>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <pthread.h>
#include "head.h"
#ifdef __GLIBC__
#define putc putc_unlocked
#define getc getc_unlocked
#endif

// same sort of thing for Windows:
#ifdef _MSC_VER
#define putc _fputc_nolock
#define getc _fgetc_nolock
#endif

#define ReadLen 300000
#define MaxNum 64

using namespace std;

typedef uint8_t  U8;
typedef uint16_t U16;
typedef uint32_t U32;
typedef uint64_t U64;
const U32 ARI_P_BITS  = 12;
const U32 DIVISOR_BITS     = 10;
const U32 PPM_P_BITS  = 22;
const U32 PPM_C_BITS  = 10;
const U32 RECIPROCAL_BITS  = 15;
const U32 PPM_C_SCALE = 32;
const U32 ARI_P_SCALE = 1 << ARI_P_BITS;  
const U32 DIVISOR_LIMIT    = 1 << DIVISOR_BITS;
const U32 RECIPROCAL_LIMIT = 1 << RECIPROCAL_BITS;
const U32 PPM_P_SCALE = 1 << PPM_P_BITS;
const U32 PPM_C_LIMIT = 1 << PPM_C_BITS;
const U32 PPM_P_MASK  = (PPM_P_SCALE  - 1) << PPM_C_BITS;
const U32 PPM_C_MASK  = PPM_C_LIMIT - 1;
const U32 PPM_P_START = PPM_P_SCALE / 2;
const U32 PPM_C_START = PPM_C_SCALE * 12;   
const U32 PPM_C_INH   = PPM_C_SCALE * 3 / 2;  
const U32 PPM_C_INC   = PPM_C_SCALE;

unsigned long long filesize, pos, control, maxsize, scanned, memo;
unsigned long rc[2], cc[2], r2, r1, conta, sce;  // Context: last 0-8 bits with a leading 1
unsigned char epix, mbit;  // 0 and 1 state count in context
unsigned short cuf[256][256];
unsigned short duf[256];
unsigned char *buffer;
unsigned char *buf;

pthread_t thread1,thread2,thread3,thread4,thread5,thread6,thread7;//thread_name
pthread_t thread[MaxNum];
//  options:
int command     = 0;   // 'c' or 'd'
int memoryLimit = 128; // memory limit in MiB
int orderLimit  = 4;   //  order limit in bytes
long totalread_threshold=20000000;
//PPM+arithmetic coding
class ReciprocalTable
{
	U16 t[DIVISOR_LIMIT];
	public:
	ReciprocalTable()
	{
		for (U32 n = 0; n < DIVISOR_LIMIT; ++n)
			t[n] = RECIPROCAL_LIMIT / (n + 2);
	}
	U32 operator[](U32 n)
	{
		assert(n < DIVISOR_LIMIT);
		return t[n];
	}
} reciprocals;

U32 Excess(U32 n, U32 m)
{
	return (n > m) ? n - m : 0;
}

struct F
{
	FILE* intput;
	FILE* output;
	bool progressBar;
}F1,F2,F3,F4,F6,F7,FQS[MaxNum]={NULL,NULL,0};

struct G
{
	char intput[50];
	char output[50];
}G1,G3,G4,G6,F5={{'\0'},{'\0'}};


U32 Divide(U32 x, U32 n, U32 y, U32 m)
{
	assert (x < (1 << n));
	assert (y < (1 << m));
	U32 dn = Excess(n, 32 - RECIPROCAL_BITS);//n-17:0
	U32 dm = Excess(m, DIVISOR_BITS);//m-15:0
	U32 dk = RECIPROCAL_BITS + dm - dn;
	return ((x >> dn) * reciprocals[y >> dm]) >> dk;
}
U32 Fit(U32 x, U32 n, U32 m)
{
	assert(x < (1 << n));
	return (n > m) ? x >> (n - m) : x << (m - n);
}

U32 Fit0(U32 x, U32 n, U32 m)
{
	assert(0 < x && x < (1 << n));
	return Fit(x, n, m) + 1 - (x >> (n - 1));
}
struct Node
{
	U32 ext0;
	U32 ext1;
	U32 sfx;
	U32 ctr;

	Node() {}

	Node(U32 ext0, U32 ext1, U32 sfx)
		: ext0(ext0),
		ext1(ext1),
		sfx(sfx),
		ctr((PPM_P_START << PPM_C_BITS) + PPM_C_START) {}

	Node(U32 sfx, const Node * sfxp)
		: ext0(0),
		ext1(0),
		sfx(sfx),
		ctr((sfxp->ctr & PPM_P_MASK) + PPM_C_INH) {}

	U32 Predict()
	{
		return ctr >> PPM_C_BITS;
	}

	template <bool bit> void Update()
	{
		U32 cnt = ctr  & PPM_C_MASK;
		U32 p1  = ctr >> PPM_C_BITS;

		if (cnt < PPM_C_LIMIT - PPM_C_INC)
			cnt += PPM_C_INC;
		else
			cnt = PPM_C_LIMIT-1;
		assert(cnt < PPM_C_LIMIT);

		if (bit)
			p1 += PPM_C_SCALE * Divide(PPM_P_SCALE - p1, PPM_P_BITS, cnt, PPM_C_BITS);
		else
			p1 -= PPM_C_SCALE * Divide(              p1, PPM_P_BITS, cnt, PPM_C_BITS);

		assert(0 < p1 && p1 < PPM_P_SCALE);

		ctr = (p1 << PPM_C_BITS) + cnt;
	}

	template <bool bit> U32 & Ext()
	{
		return bit ? ext1 : ext0;
	}
};
class PPM
{
	Node * nodes;
	Node * top;
	Node * end;

	Node * act;
	int order;

	const int nodesLimit;
	const int orderLimitBits;
	public:
	PPM()
		: nodesLimit(memoryLimit * (1 << 20) / sizeof(Node)),
		orderLimitBits(8 * orderLimit + 7)
	{
		nodes = new Node[nodesLimit];
		end = nodes + nodesLimit;
		top = nodes;
		act = nodes + 1;
		order = 0;

		*top++ = Node(1, 1, 0);                 // 1   root node
		for (int dst = 2; dst != 256; dst += 2)
			*top++ = Node(dst, dst+1, 0);       // 127 internal nodes
		for (int i = 0; i != 128; ++i)
			*top++ = Node(0, 0, 0);             // 128 leaf nodes
	}

	U32 Predict()
	{
		return Fit0(act->Predict(), PPM_P_BITS, ARI_P_BITS);
	}

	template <bool bit> void Update()
	{
		act->Update<bit>();

		Node * lst = act;
		while (act->Ext<bit>() == 0)
		{
			lst = act;
			act = nodes + act->sfx;
			order -= 8;
			act->Update<bit>();
		}

		if (act != lst && order+9 <= orderLimitBits && top < end)
		{
			Node * ext = nodes + act->Ext<bit>();
			lst->Ext<bit>() = top - nodes;
			*top = Node(ext - nodes, ext);
			act = top++;
			order += 9;
		}
		else
		{
			Node * ext = nodes + act->Ext<bit>();
			act = ext;
			order++;
		}
	}

	U32 GetUsedMemory()
	{
		return ((top - nodes) * sizeof(Node)) >> 20;
	}
};

class Encoder
{
	FILE * codeFile;
	U64 low;
	U32 range;
	U32 fluxLen;
	U8  fluxFst;
	public:
	Encoder(FILE * codeFile)
		: codeFile(codeFile),
		low(0),
		range(0xFFFFFFFF),
		fluxLen(1),
		fluxFst(0) {}

	template <bool bit> void Encode(U32 p1)
	{
		assert(0 < p1 && p1 < ARI_P_SCALE);
		U32 mid = range / ARI_P_SCALE * p1;
		if (bit)
			range = mid;
		else
			low += mid, range -= mid;
	}

	void Normalize()
	{
		while (range <= 0xFFFFFF)
		{
			U32 lo32 = low, hi32 = low >> 32;
			if (lo32 < 0xFF000000 || hi32 != 0)
			{
				putc(fluxFst + hi32, codeFile);
				while (--fluxLen)
					putc(0xFF + hi32, codeFile);
				fluxFst = lo32 >> 24;
			}
			++fluxLen;
			low = (lo32 << 8);
			range <<= 8;
		}
	}

	void FlushBuffer()
	{
		U32 lo32 = low, hi32 = low >> 32;
		putc(fluxFst + hi32, codeFile);
		while (--fluxLen)
			putc(0xFF + hi32, codeFile);
		putc(lo32 >> 24, codeFile);
		putc(lo32 >> 16, codeFile);
		putc(lo32 >>  8, codeFile);
		putc(lo32 >>  0, codeFile);
	}
};

class Decoder
{
	FILE * codeFile;
	U32 range;
	U32 cml; // code minus low
	public:
	Decoder(FILE * codeFile)
		: codeFile(codeFile),
		range(0xFFFFFFFF),
		cml(0) {}

	void FillBuffer()
	{
		for (int i = 0; i < 5; ++i)
			cml = (cml << 8) + getc(codeFile);
	}

	bool Decode(U32 p1)
	{
		assert(0 < p1 && p1 < ARI_P_SCALE);
		U32 mid = range / ARI_P_SCALE * p1;
		if (cml < mid)
		{
			range = mid;
			return 1;
		}
		else
		{
			cml -= mid, range -= mid;
			return 0;
		}
	}

	void Normalize()
	{
		while (range <= 0xFFFFFF)
		{
			cml = (cml << 8) + getc(codeFile);
			range <<= 8;
		}
	}
};

class ProgressBar
{
	static const int period = 1 << 18;
	clock_t start;
	public:
	ProgressBar()
	{
		start = clock();
	}

	void Display(U64 processed, U64 total, U64 memory)
	{
		if (processed % period != 0)
			return;

		int percentage = (100 * processed + total / 2) / total;
		printf("\r%3d%% ", percentage);

		const char blocks[] = "[########################################]";
		const char spaces[] = "[                                        ]";
		int maxBlocks = 40;
		int numBlocks = (maxBlocks * processed + total / 2) / total;
		int fromBlocks = numBlocks + 1;
		int fromSpaces = maxBlocks + 1 - numBlocks;
		fwrite(blocks, fromBlocks, 1, stdout);
		fwrite(spaces + fromBlocks, fromSpaces, 1, stdout);

		clock_t current = clock();
		int speed = processed / 1024 * CLOCKS_PER_SEC / (current - start + 1);
		//printf("%6d kiB/s %d/%d MiB", speed, memory, memoryLimit);

		fflush(stdout);
	}

	void Finish(U64 textLength, U64 codeLength)
	{
		double seconds = (double)(clock() - start) / CLOCKS_PER_SEC;
		double bpc = 8.0 * codeLength / textLength;

		if (command == 'd')
			swap(textLength, codeLength);

		putc('\n', stdout);
		//printf("%d -> %d, %.2f s, %.3f bpc.\n",textLength, codeLength, seconds, bpc);
	}
};





void Compress(FILE * textFile, FILE * codeFile, bool progressBar)
{
	fseek(textFile, 0, SEEK_END);
	U32 textLength = ftell(textFile);
	U64 textLength1 = ftell(textFile);
	fseek(textFile, 0, SEEK_SET);
	putc(textLength1 >> 56, codeFile);
	putc(textLength1 >> 48, codeFile);
	putc(textLength1 >> 40, codeFile);
	putc(textLength1 >> 32, codeFile);
	putc(textLength1 >> 24, codeFile);
	putc(textLength1 >> 16, codeFile);
	putc(textLength1 >> 8, codeFile);
	putc(textLength1 >> 0, codeFile);

	ProgressBar bar;
	Encoder rc(codeFile);
	PPM ppm;
	for (U64 processed = 0; processed != textLength1; ++processed)
	{
		if (progressBar == 1)bar.Display(processed, textLength1, ppm.GetUsedMemory());

		U32 c = getc(textFile);
		for (U32 mask = 1 << 7; mask != 0; mask >>= 1)
		{
			U32 p1 = ppm.Predict();
			if (c & mask)
			{
				rc.Encode<1>(p1);
				ppm.Update<1>();
			}
			else
			{
				rc.Encode<0>(p1);
				ppm.Update<0>();
			}
			rc.Normalize();
		}
	}
	rc.FlushBuffer();
	if (progressBar == 1)bar.Finish(textLength1, ftell(codeFile));
}

void Decompress(FILE * codeFile, FILE * textFile, bool progressBar)
{
	U32 textLength1 = 0, textLength2 = 0;
	U64 textLength = 0;
	textLength1 += getc(codeFile) << 24;
	textLength1 += getc(codeFile) << 16;
	textLength1 += getc(codeFile) << 8;
	textLength1 += getc(codeFile) << 0;
	textLength2 += getc(codeFile) << 24;
	textLength2 += getc(codeFile) << 16;
	textLength2 += getc(codeFile) << 8;
	textLength2 += getc(codeFile) << 0;
	textLength += textLength1;
	textLength = textLength << 32;
	textLength += textLength2;
	ProgressBar bar;
	Decoder rc(codeFile);
	PPM ppm;
	rc.FillBuffer();
	for (U64 processed = 0; processed != textLength; ++processed)
	{
		if (progressBar == 1)bar.Display(processed, textLength, ppm.GetUsedMemory());

		U32 c = 0;
		for (U32 mask = 1 << 7; mask != 0; mask >>= 1)
		{
			U32 p1 = ppm.Predict();
			if (rc.Decode(p1))
			{
				ppm.Update<1>();
				c |= mask;
			}
			else
			{
				ppm.Update<0>();
			}
			rc.Normalize();
		}
		putc(c, textFile);
	}
	if (progressBar == 1)bar.Finish(textLength, ftell(codeFile));
}


//addzip
unsigned long len(unsigned long number)
{
	unsigned long nbit;
	nbit = 1;
bscan:
	if ((number >> 1)>0)
	{
		number = number >> 1;
		nbit++;
		goto bscan;
	}
	return (nbit);
}
unsigned char Peekb(unsigned char *memory)
{
	unsigned char result;
	memcpy(&result, memory, 1);
	return (result);
}
void Pokeb(unsigned char *memory, unsigned char value)
{
	memcpy(memory, &value, 1);
}
int p()
{
	return 4096 * (Peekb(buf + rc[1]) + 1) / (Peekb(buf + rc[1]) + Peekb(buf + rc[0]) + 2);
}

void update(int y)
{
	if (Peekb(buf + rc[y]) + 8<255)
		Pokeb(buf + rc[y], Peekb(buf + rc[y]) + 8);
	else
	{
		Pokeb(buf + rc[y], Peekb(buf + rc[y]) >> 1);
		Pokeb(buf + rc[!y], Peekb(buf + rc[!y]) >> 1);
	}
	rc[y] >>= 1;          
	rc[y] += (1 << mbit);  
	rc[!y] >>= 1;         

}


typedef enum { COMPRESS, DECOMPRESS } Mode;
class Addzip {
	private:

		const Mode mode;       // Compress or decompress?
		FILE* archive;         // Compressed data file
		U32 x1, x2;            // Range, initially [0, 1), scaled by 2^32
		U32 x;                 // Last 4 input bytes of archive.
	public:
		Addzip(Mode m, FILE* f);
		void addzip(int y);    // Compress bit y
		int decode();          // Uncompress and return bit y
		void flush();          // Call when done compressing
};

// Constructor
Addzip::Addzip(Mode m, FILE* f) : mode(m), archive(f), x1(0),
	x2(0xffffffff), x(0) {

		// In DECOMPRESS mode, initialize x to the first 4 bytes of the archive
		if (mode == DECOMPRESS) {
			for (int i = 0; i<4; ++i) {
				int c = getc(archive);
				if (c == EOF) c = 0;
				x = (x << 8) + (c & 0xff);
			}
		}
	}

/* addzip(y) -- addzip bit y by splitting the range [x1, x2] in proportion
   to P(1) and P(0) as given by the predictor and narrowing to the appropriate
   subrange.  Output leading bytes of the range as they become known. */

inline void Addzip::addzip(int y) {

	// Update the range
	const U32 xmid = x1 + ((x2 - x1) >> 12) * p();
	assert(xmid >= x1 && xmid < x2);
	if (y)
		x2 = xmid;
	else
		x1 = xmid + 1;

	update(y);

	// Shift equal MSB's out
	while (((x1^x2) & 0xff000000) == 0) {
		putc(x2 >> 24, archive);
		x1 <<= 8;
		x2 = (x2 << 8) + 255;
	}
}

/* Decode one bit from the archive, splitting [x1, x2] as in the Addzip
   and returning 1 or 0 depending on which subrange the archive point x is in.
 */
inline int Addzip::decode() {

	// Update the range
	const U32 xmid = x1 + ((x2 - x1) >> 12) * p();
	assert(xmid >= x1 && xmid < x2);
	int y = 0;
	if (x <= xmid) {
		y = 1;
		x2 = xmid;
	}
	else
		x1 = xmid + 1;

	update(y);

	// Shift equal MSB's out
	while (((x1^x2) & 0xff000000) == 0) {
		x1 <<= 8;
		x2 = (x2 << 8) + 255;
		int c = getc(archive);
		if (c == EOF) c = 0;
		x = (x << 8) + c;
	}
	return y;
}

// Should be called when there is no more to compress
void Addzip::flush() {

	// In COMPRESS mode, write out the remaining bytes of x, x1 < x < x2
	if (mode == COMPRESS) {
		while (((x1^x2) & 0xff000000) == 0) {
			putc(x2 >> 24, archive);
			x1 <<= 8;
			x2 = (x2 << 8) + 255;
		}
		putc(x2 >> 24, archive);  // First unequal byte
	}
}
void endx()
{

	exit(1);
}


int AddEncode(int mode, char inputFile[], char outputFile[]) {
	//mode 0:compress , 1:decompress


	// Start timer
	clock_t start = clock();
	FILE *in;
	FILE *out;
	// Open files
	if (mode == 0)
	{
		mbit = 29;
		in = fopen(inputFile, "rb");
		if (!in) perror(inputFile), exit(1);
		out = fopen(outputFile, "wb");
		if (!out) perror(outputFile), exit(1);
	}
	else
	{
		in = fopen(inputFile, "rb");
		if (!in) perror(inputFile), exit(1);
		out = fopen(outputFile, "wb");
		if (!out) perror(outputFile), exit(1);
	}

	int c;
	maxsize = 1 << 20;
	buffer = (unsigned char *)malloc(maxsize);
	memset(cuf, 0, sizeof(cuf));
	memset(duf, 0, sizeof(duf));
	// Compress
	if (mode == 0)
	{
		memo = 1 << (mbit + 1);
		buf = (unsigned char *)malloc(memo);
		memset(buf, 0, sizeof(buf));
		while ((c = getc(in)) != EOF) // reduction of bits to encode
		{
			if (cuf[r1][c] == 0)
			{
				cuf[r1][c] = 1; duf[r1]++;
			}
			r1 = c;
		}
		rewind(in);
		fseek(in, 0, SEEK_END);
		filesize = ftell(in);
		rewind(in);
		fseek(out, 0, SEEK_SET);
		fwrite(&filesize, 1, 8, out);
		fwrite(&mbit, 1, 1, out);
		Addzip e(COMPRESS, out);
		for (r1 = 0; r1<256; r1++)
		{
			if (duf[r1]>0)
			{
				e.addzip(1);
				conta = 0;
				for (r2 = 0; r2<256; r2++)
				{
					e.addzip(cuf[r1][r2]);
					if (cuf[r1][r2] == 1){ cuf[r1][r2] = conta; conta++; }
				}
			}
			else
				e.addzip(0);
		}
		memset(buf, 0, sizeof(buf));
		r1 = 0;
newscan:
		if (scanned + maxsize>filesize) maxsize = filesize - scanned;
		pos = 0;
		fread(buffer, 1, maxsize, in);
		while (pos<maxsize)
		{
			c = Peekb(buffer + pos);
			for (int i = len(duf[r1] - 1) - 1; i >= 0; --i)
				e.addzip((cuf[r1][c] >> i) & 1);
			r1 = c;
			pos++;
			if (scanned + pos>control*(filesize / 79))
			{
				//printf(">");
				control++;
			}
		}
		scanned += maxsize;
		if (scanned<filesize) goto newscan;
		e.flush();
	}

	// Decompress
	else {
		fseek(in, 0, SEEK_SET);
		fread(&filesize, 1, 8, in);
		fread(&mbit, 1, 1, in);
		memo = 1 << (mbit + 1);
		buf = (unsigned char *)malloc(memo);
		memset(buf, 0, sizeof(buf));
		Addzip e(DECOMPRESS, in);
		for (r1 = 0; r1<256; r1++)
		{
			if (e.decode())
				for (r2 = 0; r2<256; r2++)
				{
					sce = e.decode();
					if (sce == 1)
					{
						cuf[r1][duf[r1]] = r2;
						duf[r1]++;
					}
				}
		}
		memset(buf, 0, sizeof(buf));
		r1 = 0;
newdscan:
		if (scanned + maxsize>filesize) maxsize = filesize - scanned;
		pos = 0;
		while (pos<maxsize)
		{
			c = 0;
			for (int i = len(duf[r1] - 1) - 1; i >= 0; --i)
				c += e.decode()*(1 << i);
			Pokeb(buffer + pos, cuf[r1][c]);
			r1 = cuf[r1][c];
			if (scanned + pos>control*(filesize / 79))
			{
				//printf(">");
				control++;
			}
			pos++;
		}
		fwrite(buffer, 1, maxsize, out);
		scanned += maxsize;
		if (scanned<filesize) goto newdscan;
	}
	// Print results
	/*
	   if (mode == 0)
	   {
	   printf("%ld -> %ld, in %1.2f s.\n"
	   ,ftell(in), ftell(out)
	   , ((double)clock() - start) / CLOCKS_PER_SEC);
	   }
	   else
	   {
	   printf("%ld -> %ld, in %1.2f s.\n"
	   ,  ftell(in),  ftell(out)
	   , ((double)clock() - start) / CLOCKS_PER_SEC);
	   }
	 */
	return 0;
}






char * optarg = NULL;
int    optind = 0;

int getopt(int argc, char ** argv, const char * spec)
{
	static char * next = NULL;
	static int end = 0;

	if (next == NULL || next[0] == '\0')
	{
		rotate(argv + optind, argv + end, argv + end + 1);
		++optind;
		++end;

		while (end != argc && (argv[end][0] != '-' || argv[end][1] == '\0'))
			++end;

		if (end == argc)
			return -1;

		next = &argv[end][1];
	}

	const char * opt = strchr(spec, next[0]);

	if (opt == NULL)
	{
		fprintf(stderr, "%s: invalid option '-%c'\n",
				argv[0], next[0]);
		optarg = NULL;
		++next;
		return '?';
	}
	else if (opt[1] == ':' && next[1] == '\0')
	{
		fprintf(stderr, "%s: missing argument for option '-%c'\n",
				argv[0], next[0]);
		optarg = NULL;
		next = NULL;
		return '?';
	}
	else if (opt[1] == ':' && next[1] != '\0')
	{
		optarg = next + 1;
		next = NULL;
		return opt[0];
	}
	else
	{
		optarg = NULL;
		++next;
		return opt[0];
	}
}



void*   Compress_encode_imp_thread(void * args) 
{
	F* arg=(F*)args;
	Compress(arg->intput, arg->output, arg->progressBar);
	return NULL;
}



void*   Compress1_encode_imp_thread(void * args) 
{
	G* arg=(G*)args;
	if (0 != (AddEncode(0, arg->intput, arg->output)))printf("add.txt compression occur problem!");
	return NULL;
}


void*   Compress2_encode_imp_thread(void * args) 
{
	G* arg=(G*)args;
	char order[500]={'\0'};
	sprintf(order, "./zpaq a %s %s -method 5 -threads 8", arg->output, arg->intput);
	if (-1 == (system(order)))printf("Please make sure you have assigned the executable permission to lzip tool correctly");
	return NULL;
}


void*   Compress3_encode_imp_thread(void * args) 
{
	G* arg=(G*)args;
	char order[500]={'\0'};
	sprintf(order, "./lpaq9m 6 %s %s ", arg->intput, arg->output);
	if (-1 == (system(order)))printf("Please make sure you have assigned the executable permission to lpaq8 tool correctly");
	return NULL;
}


void*   Decompress_encode_imp_thread(void * args)
{
	F* arg=(F*)args;
	Decompress(arg->intput, arg->output, arg->progressBar);
	return NULL;
}


void*   Decompress1_encode_imp_thread(void * args) 
{
	G* arg=(G*)args;
	if (0 != (AddEncode(1, arg->intput, arg->output)))printf("add.txt decompression occur problem!");
	return NULL;
}


void*   Decompress2_encode_imp_thread(void * args) 
{
	G* arg=(G*)args;
	char order[500]={'\0'};
	sprintf(order, "./zpaq e %s -threads 8 -force ", arg->intput);
	if (-1 == (system(order)))printf("Please make sure you have assigned the executable permission to lzip tool correctly");
	return NULL;
}


void*   Decompress3_encode_imp_thread(void * args) 
{
	G* arg=(G*)args;
	char order[500]={'\0'};
	sprintf(order, "./lpaq9m d %s %s ", arg->intput,arg->output);
	if (-1 == (system(order)))printf("Please make sure you have assigned the executable permission to lpaq8 tool correctly");
	return NULL;
}

char *get_fastq_name(char *fastqname, char *str)
{
	int i = 0, j = 0;
	for (i = 0; i<60; i++)
	{
		if (str[i] == '.'&&str[i + 1] == 'f' && str[i + 2] == 'a' )
			break;
		else fastqname[j++] = str[i];
	}
	return fastqname;
}


// MAIN
int main(int argc, char ** argv)
{

	bool help = false;
	bool version = false;

	int c;
	while ((c = getopt(argc, argv, "hVvqm:O:")) != -1)
	{
		if      (c == 'h') help    = true;
		else if (c == 'V') version = true;
		else if (c == 'm' || c == 'O')
		{
			errno = 0;
			char * rest;
			long val = strtol(optarg, &rest, 10);
			if (errno != 0 || *rest != '\0' || val < 0)
			{
				fprintf(stderr,
						"%s: invalid argument '%s' for option '%c'\n",
						argv[0], optarg, c);
				return 1;
			}
			if (c == 'm') memoryLimit = val;
			else          orderLimit  = val;
		}
		else return 1;
	}

	help = help || (!version && optind == argc);

	if (help)
	{
		puts("To compress a file invoke the program like this\n"
				"  FQZip c INPUT OUTPUT\n"
				"To decompress\n"
				"  FQZip d INPUT OUTPUT\n"
				"Existing output files are overwritten.\n"
				"\n"
		    );
	}


	if ((argv[optind][0] != 'c' && argv[optind][0] != 'd') ||
			argv[optind][1] != 0)
	{
		fprintf(stderr, "%s: unrecognized command '%s'\n",
				argv[0], argv[optind]);
		return 1;
	}
	
	int highestCom_flag;
	command = argv[optind][0];
	FILE * input1 = fopen(argv[optind+1], "rb");
	FILE * input2 = fopen(argv[optind+2], "rb");
	FILE * input3 = fopen(argv[optind+3], "rb");
	FILE * input4 = fopen(argv[optind+4], "rb");
	FILE * input6 = fopen(argv[optind+6], "rb");
	FILE * output1 = fopen(argv[optind+7], "wb");
	FILE * output2 = fopen(argv[optind+8], "wb");
	FILE * output3 = fopen(argv[optind+9], "wb");
	FILE * output4 = fopen(argv[optind+10], "wb");
	FILE * output6 = fopen(argv[optind+12], "wb");
	char fasta_intput[50]="\0", fasta_output[50]="\0";
	int file_number;

	F1.intput=input1;
	F1.output=output1;
	F2.intput=input2;
	F2.output=output2;
	F3.intput=input3;
	F3.output=output3;
	F4.intput=input4;
	F4.output=output4;
	strcpy(F5.intput,argv[optind+5]);
	strcpy(F5.output,argv[optind+11]);
	F6.intput=input6;
	F6.output=output6;
	
	strcpy(G1.intput,argv[optind+1]);
	strcpy(G1.output,argv[optind+7]);
	strcpy(G3.intput,argv[optind+3]);
	strcpy(G3.output,argv[optind+9]);
	strcpy(G4.intput,argv[optind+4]);
	strcpy(G4.output,argv[optind+10]);
	strcpy(G6.intput,argv[optind+6]);
	strcpy(G6.output,argv[optind+12]);	


	if (command == 'c')
	{		
		unsigned long totalread;
		if(argc==18)//assemble-based model
		highestCom_flag=atoi(argv[optind+16]);
		if(argc==17)//reference-based model
		highestCom_flag=atoi(argv[optind+15]);
		file_number=atoi(argv[optind+14]);
		if(atol(argv[optind+13])==0)
		printf("This Fastq file is too big to condcut the compression!");
		totalread=atol(argv[optind+13]);//obtain total read count
		if(totalread<=totalread_threshold||highestCom_flag==1)//20 milllion=~5GB all file size,<20 million don't cread multithread for qs compression
		{
		if(highestCom_flag==0)// normal model
		{
		F6.progressBar = 1;
		pthread_create(&thread1, 0, &Compress_encode_imp_thread, &F1);
		pthread_create(&thread2, 0, &Compress_encode_imp_thread, &F2);
		pthread_create(&thread3, 0, &Compress_encode_imp_thread, &F3);
		pthread_create(&thread4, 0, &Compress_encode_imp_thread, &F4);
		pthread_create(&thread5, 0, &Compress1_encode_imp_thread, &F5);
		pthread_create(&thread6, 0, &Compress_encode_imp_thread, &F6);
		}
		else	//highest compression model
		{
		pthread_create(&thread1, 0, &Compress3_encode_imp_thread, &G1);
		pthread_create(&thread2, 0, &Compress_encode_imp_thread, &F2);
		pthread_create(&thread3, 0, &Compress3_encode_imp_thread, &G3);
		pthread_create(&thread4, 0, &Compress3_encode_imp_thread, &G4);
		pthread_create(&thread5, 0, &Compress2_encode_imp_thread, &F5);
		pthread_create(&thread6, 0, &Compress2_encode_imp_thread, &G6);
		}
		if(argc==18)//assemble-based model
		{
			strcpy(fasta_intput,argv[optind+15]);
			strcpy(fasta_output,argv[optind+15]);
			strcat(fasta_output,".lz");
			FILE * input7=fopen(fasta_intput, "rb");
			FILE * output7=fopen(fasta_output, "wb");
			F7.intput=input7;
			F7.output=output7;
			pthread_create(&thread7, 0, &Compress_encode_imp_thread, &F7);
		}
		pthread_join(thread1, 0);
		pthread_join(thread2, 0);
		pthread_join(thread3, 0);
		pthread_join(thread4, 0);
		pthread_join(thread5, 0);
		pthread_join(thread6, 0);
		if(argc==18)
		pthread_join(thread7, 0);
		}
		else	//normal model
		{
			pthread_create(&thread1, 0, &Compress_encode_imp_thread, &F1);
			pthread_create(&thread2, 0, &Compress_encode_imp_thread, &F2);
			pthread_create(&thread3, 0, &Compress_encode_imp_thread, &F3);
			pthread_create(&thread4, 0, &Compress_encode_imp_thread, &F4);
			pthread_create(&thread5, 0, &Compress1_encode_imp_thread, &F5);
			if(argc==18)//assemble-based model
			{
			strcpy(fasta_intput,argv[optind+15]);
			strcpy(fasta_output,argv[optind+15]);
			strcat(fasta_output,".lz");
			FILE * input7=fopen(fasta_intput, "rb");
			FILE * output7=fopen(fasta_output, "wb");
			F7.intput=input7;
			F7.output=output7;
			pthread_create(&thread7, 0, &Compress_encode_imp_thread, &F7);
			}		
			int i;
			FILE *fp_qs = fopen(argv[optind + 6], "r");

			char tempStr[ReadLen] = {'\0'};
			long subreadcounts = totalread / file_number;//file_number thread to compress qs.txt
			long lastblockcounts = 0;
			long subline;
			char index[10];
			char input_name[500] = { '\0' };
			char output_name[500]= {'\0'};
			if (totalread%file_number == 0){ lastblockcounts = subreadcounts; }
			else{ lastblockcounts = subreadcounts + (totalread - subreadcounts*file_number); }
			FILE *subfastq[file_number];
			for (i = 0; i < file_number - 1; i++)
			{
				subline = 0;
				sprintf(input_name, "%s_%d", argv[optind + 6], i+1);
				subfastq[i] = fopen(input_name, "wb");
				while (fgets(tempStr, ReadLen, fp_qs) != NULL)
				{
					subline++;
					fputs(tempStr, subfastq[i]);
					if (subline == subreadcounts)
						break;
				}
						
				fclose(subfastq[i]);
		
				sprintf(output_name, "%s_%d.lz", argv[optind + 6],i+1);
				FILE * input_qs = fopen(input_name, "rb");
				FILE * output_qs = fopen(output_name, "wb");
				FQS[i].intput = input_qs;
				FQS[i].output = output_qs;
				pthread_create(&thread[i], 0, &Compress_encode_imp_thread, &FQS[i]);
				
			}
			//To process the last subfile
			subline = 0;
			sprintf(input_name, "%s_%d", argv[optind + 6], file_number);
			subfastq[file_number-1] = fopen(input_name, "wb");
			while (fgets(tempStr, ReadLen, fp_qs) != NULL)
			{
				subline++;
				fputs(tempStr, subfastq[file_number - 1]);
			}
			fclose(subfastq[file_number - 1]);
			sprintf(output_name, "%s_%d.lz", argv[optind + 6], file_number);
			FILE * input_qs = fopen(input_name, "rb");
			FILE * output_qs = fopen(output_name, "wb");
			FQS[file_number-1].intput = input_qs;
			FQS[file_number-1].output = output_qs;
			FQS[file_number-1].progressBar = 1;
			pthread_create(&thread[file_number - 1], 0, &Compress_encode_imp_thread, &FQS[file_number-1]);
			pthread_join(thread1, 0);
			pthread_join(thread2, 0);
			pthread_join(thread3, 0);
			pthread_join(thread4, 0);
			if(argc==18)
			pthread_join(thread7, 0);
			for (i = 0; i < file_number;i++)
			pthread_join(thread[i], 0);	
			pthread_join(thread5, 0);
			
		}
	}
	if (command == 'd')
	{

		char input_name[500] = { '\0' };
		char output_name[500] = { '\0' };
		char tempStr[ReadLen];
		char fastq_name[500] = { '\0' };
		int thread_num = 0;
		get_fastq_name(fastq_name, argv[optind + 6]);
		bool ismulti = false;
		if(argc==16)//assemble-based model
		highestCom_flag=atoi(argv[optind+14]);
		if(argc==15)//reference-based model
		highestCom_flag=atoi(argv[optind+13]);
		strcat(tempStr, fastq_name);
		strcat(tempStr, ".fastq.qs.txt_1.lz");
		FILE *fp_qs;
		fp_qs = fopen(tempStr, "r");
		if (fp_qs != NULL)
			ismulti = true;
		if(ismulti==false||highestCom_flag==1)
		{
		if(highestCom_flag==0)
		{
		F6.progressBar = 1;
		pthread_create(&thread1, 0, &Decompress_encode_imp_thread, &F1);
		pthread_create(&thread2, 0, &Decompress_encode_imp_thread, &F2);
		pthread_create(&thread3, 0, &Decompress_encode_imp_thread, &F3);
		pthread_create(&thread4, 0, &Decompress_encode_imp_thread, &F4);
		pthread_create(&thread5, 0, &Decompress1_encode_imp_thread, &F5);
		pthread_create(&thread6, 0, &Decompress_encode_imp_thread, &F6);
		}
		else
		{
		pthread_create(&thread1, 0, &Decompress3_encode_imp_thread, &G1);
		pthread_create(&thread2, 0, &Decompress_encode_imp_thread, &F2);
		pthread_create(&thread3, 0, &Decompress3_encode_imp_thread, &G3);
		pthread_create(&thread4, 0, &Decompress3_encode_imp_thread, &G4);
		pthread_create(&thread5, 0, &Decompress2_encode_imp_thread, &F5);
		pthread_create(&thread6, 0, &Decompress2_encode_imp_thread, &G6);
		}
		if(argc==16)//assemble-based model
		{
			strcpy(fasta_intput,argv[optind+13]);
			strcat(fasta_intput,".lz");
			strcpy(fasta_output,argv[optind+13]);
			FILE * input7=fopen(fasta_intput, "rb");
			FILE * output7=fopen(fasta_output, "wb");
			F7.intput=input7;
			F7.output=output7;
			pthread_create(&thread7, 0, &Decompress_encode_imp_thread, &F7);
		}
		pthread_join(thread1, 0);
		pthread_join(thread2, 0);
		pthread_join(thread3, 0);
		pthread_join(thread4, 0);
		pthread_join(thread5, 0);
		pthread_join(thread6, 0);
		if(argc==16)
		pthread_join(thread7, 0);
		}
		else
		{	
			pthread_create(&thread1, 0, &Decompress_encode_imp_thread, &F1);
			pthread_create(&thread2, 0, &Decompress_encode_imp_thread, &F2);
			pthread_create(&thread3, 0, &Decompress_encode_imp_thread, &F3);
			pthread_create(&thread4, 0, &Decompress_encode_imp_thread, &F4);
			pthread_create(&thread5, 0, &Decompress1_encode_imp_thread, &F5);
			for (int i = 2; i <= MaxNum; i++)
			{
				
				sprintf(tempStr, "%s.fastq.qs.txt_%d.lz", fastq_name, i);
				fp_qs = fopen(tempStr, "r");
				if (fp_qs == NULL)
				{
					break;
				}
				thread_num = i;
				fclose(fp_qs);
			}
			
			printf("Detect this decompress file has %d block, start multithread decompression.\n", thread_num);
			for (int i = 1; i <= thread_num; i++)
			{
				sprintf(input_name, "%s.fastq.qs.txt_%d.lz", fastq_name, i);
				sprintf(output_name, "%s.fastq.qs.txt_%d", fastq_name, i);
				FILE * input_qs = fopen(input_name, "rb");
				FILE * output_qs = fopen(output_name, "wb");
				FQS[i - 1].intput = input_qs;
				FQS[i - 1].output = output_qs;
				if (i==thread_num)
					FQS[i - 1].progressBar = 1;
				pthread_create(&thread[i-1], 0, &Decompress_encode_imp_thread, &FQS[i-1]);
				
			}
			pthread_join(thread1, 0);
			pthread_join(thread2, 0);
			pthread_join(thread3, 0);
			pthread_join(thread4, 0);
			for (int i = 0; i <thread_num; i++)
			{
			pthread_join(thread[i], 0);
			fclose(FQS[i].intput);
			fclose(FQS[i].output);
			}
			char order[500];
			stpcpy(order,"cat ");

			for (int i = 1; i <= thread_num; i++)
			{
				sprintf(tempStr, " %s.fastq.qs.txt_%d ", fastq_name, i);
				strcat(order,tempStr);
			}
			sprintf(tempStr, " >> %s.fastq.qs.txt", fastq_name);
			strcat(order,tempStr);
			if (-1 == (system(order)))printf("can't merge multi qs blocks!\n");
			//fclose(fp_qs);
			pthread_join(thread5, 0);
			if(argc==16)//assemble-based model
			{
			strcpy(fasta_intput,argv[optind+13]);
			strcat(fasta_intput,".lz");
			strcpy(fasta_output,argv[optind+13]);
			FILE * input7=fopen(fasta_intput, "rb");
			FILE * output7=fopen(fasta_output, "wb");
			F7.intput=input7;
			F7.output=output7;
			pthread_create(&thread7, 0, &Decompress_encode_imp_thread, &F7);
			pthread_join(thread7, 0);
			}			
			sprintf(order, "rm %s.fastq.qs.txt_* ", fastq_name);
			if (-1 == (system(order)))printf("can't delete multi qs block intermediate files\n");
			
			
		}
		
	}
}
