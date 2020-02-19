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

/**
 * Generates Vector with random binary entries.
 *
 * @param len Length of Vector.
 * @param data Reference to the vector that the result will be written to.
 */
static void generate_random_vector(int len, uint8_t *data)
{
	int i;
	for (i = 0; i < len; i++)
    {
        data[i] = rand()%2;
    }
}

/**
 * Bitwise XOR message with given vector
 *
 * @param len Length of message (also vector and data).
 * @param message Reference to the message vector.
 * @param vector Reference to the vector containing random bits.
 * @param data Reference to the vector that the result will be written to.
 */
static void generate_xor_vector(int len, uint8_t *message, uint8_t *vector, uint8_t *data)
{
	int i;
	for (i = 0;	i < len; i++) {
        data[i] = message[i] ^ vector[i];
    }
}
/**
 * Bitwise Multiplication data1 with data2
 *
 * @param len Length of message (also vector and data).
 * @param data1 Reference to the message vector.
 * @param data2 Reference to the vector containing random bits.
 * @param result Reference to the vector that the result will be written to.
 */
static void generate_mul_vector(int len, uint8_t *data1, uint8_t *data2, uint8_t *result)
{
	int i;
	for (i = 0;	i < len; i++) {
        result[i] = data1[i] * data2[i];
    }
}
/**
 * Flipping Bits in the given data vector.
 *
 * @param bitflip Vector containing the positions at which the bits will be flipped.
 * @param ncorrupt Number of Errors (Length of bitflip).
 * @param data Reference to the vector that the result will be written to.
 */
static void corrupt_data(int *bitflip, uint8_t *data, int ncorrupt)
{
	int i;
	for (i = 0; i < ncorrupt; i++) {
        if(data[bitflip[i]] == 1){
            data[bitflip[i]] = 0;
        } else {
            data[bitflip[i]] = 1;
        }
	}
}

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
