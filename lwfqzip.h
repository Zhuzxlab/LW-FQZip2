#ifndef _MAIN_H
#define _MAIN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include <time.h>
#include <stddef.h>
#include <pthread.h>

#include "compress.h"
#include "decompress.h"

size_t read_len=300000;

void print_help();
void print_version();

#endif // _MAIN_H
