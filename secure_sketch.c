/*
 * Original from: https://github.com/Parrot-Developers/bch/blob/master/Documentation/bch/tu_correct.c
 */
 
/*
 * Copyright (C) 2011 Parrot S.A.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

// static int rev8(int bit)
// {
// 	return (bit & ~7)|(7-(bit & 7));
// }

// static int compar(const void *a, const void *b) {
//     return *(int*)a - *(int*)b;
// }

/**
 * Introduce up to t errors, generates a vector containing positions at which errors will occur.
 *
 * @param len Length of Vector.
 * @param bit Reference to the vector that the result will be written to.
 * @param size Number of errors.
 * @param seed Seed for randomization.
 */
// static void generate_error_vector(int len, int *bit, int size, unsigned int seed)
// {
// 	int i, done = 0;

// 	/* corrupt data */
// 	srand48(seed);

// 	while (!done) {
// 		for (i = 0; i < size; i++) {
// 			bit[i] = lrand48() % len;
// 		}
// 		qsort(bit, size, sizeof(int), compar);
// 		done = 1;
// 		for (i = 0; i < size-1; i++) {
// 			if (bit[i] == bit[i+1]) {
// 				if ((i > 0) && (bit[i-1]+1 < bit[i])) {
// 					bit[i]--;
// 					continue;
// 				}
// 				if ((i+1 < size-1) && (bit[i+1]+1 < bit[i+2])) {
// 					bit[i+1]++;
// 					continue;
// 				}
// 				done = 0;
// 				break;
// 			}
// 		}
// 	}
// 	for (i = 0; i < size; i++) {
// 		bit[i] = rev8(bit[i]);
// 	}
// }



// char *bin2hex(const unsigned char *bin, size_t len)
// {
// 	char   *out;
// 	size_t  i;

// 	if (bin == NULL || len == 0)
// 		return NULL;

// 	out = malloc(len*2+1);
// 	for (i=0; i<len; i++) {
// 		out[i*2]   = "0123456789ABCDEF"[bin[i] >> 4];
// 		out[i*2+1] = "0123456789ABCDEF"[bin[i] & 0x0F];
// 	}
// 	out[len*2] = '\0';

// 	return out;
// }
