#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <inttypes.h>

#include "lib/bch/bch_codec.h"
#include "parameter_estimator.h"
#include "alice_bob_phase.c"
// #include "2019-07-10_16-03-24_csi_log.c"
// #include "2019-07-10_15-12-45_csi_log.c"
//~ #include "2020-02-19_14-05-25_csi_log.c"

#define MAX_ERRORS   2048
#define WINDOW_SIZE     4
#define DEBUG_KEYGEN    1
#define	POLY_DEGREE     7   /* BCH polynomial degree */
#define ERROR_CORR_CAP  9   /* erroc correction capability */

/**
 * Data is an array containing 2 * sketch_len + ecc_len bytes.
 * It can be considered a struct with minimal alignment of the following form:
 *
 * struct data {
 *  uint8_t sketch[sketch_len];
 *  uint8_t rand_x[rand_x_len];
 *  uint8_t ecc[ecc_len];
 * }__attribute__((packed));
 *
 * TODO: Use bits instead of bytes.
 */
struct key_recovery_data {
    uint16_t sketch_len;
    uint16_t rand_x_len;
    uint16_t ecc_len;
    uint8_t data[];
};

/**
 * Helper function to print the various uint8_t arrays.
 *
 * \param array The array of uint8_t to print.
 * \param length The length of the array to print.
 */
inline static void print_array(uint8_t array[], int length)
{
    for (int i = 0; i < length; ++i) {
        printf("%" PRIu8, array[i]);
    }
}

inline static unsigned int calculate_nr_block_bits(struct bch_control *bch) { return (bch->n + 1) / (2 * ERROR_CORR_CAP - 2); }
inline static unsigned int calculate_nr_blocks(struct bch_control *bch) { return (bch->n + 1) / calculate_nr_block_bits(bch); }

/**
 * Generates Vector with random binary entries.
 *
 * @param len Length of Vector.
 * @param data Reference to the vector that the result will be written to.
 */
inline static void generate_random_vector(uint8_t *data, size_t len)
{
	for (size_t i = 0; i < len; i++) {
        data[i] = rand() % 2;
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
inline static void generate_xor_vector(int len, uint8_t *message, uint8_t *vector, uint8_t *data)
{
	for (int i = 0;	i < len; ++i) {
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
inline static void generate_mul_vector(int len, uint8_t *data1, uint8_t *data2, uint8_t *result)
{
	for (int i = 0;	i < len; ++i) {
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
inline static void corrupt_data(unsigned int *bitflip, uint8_t *data, unsigned int ncorrupt)
{
	for (unsigned int i = 0; i < ncorrupt; ++i) {
        if(data[bitflip[i]] == 1){
            data[bitflip[i]] = 0;
        } else {
            data[bitflip[i]] = 1;
        }
	}
}

/**
 *
 */
static double calculate_metric_entropy(uint8_t message[], unsigned int length)
{
    /* H(X) = [-[(freq_0)log2(freq_0))+((freq_1)log2(freq_1))]] / number of bits */
	double freq_0 = 0.0;
	double freq_1 = 0.0;
	double c = 0.0;

    /* Calculate frequency of each alphabet */
	for (unsigned int i = 0; i < length; ++i) {
		if (message[i] == 0)
			++c;
	}
	freq_0 = c / length;
	printf("freq_0: %f\n", freq_0);

	c = 0.0;
	for (unsigned int i = 0; i < length; ++i) {
		if (message[i] == 1)
			++c;
	}
	freq_1 = c / length;
	printf("freq_1: %f\n", freq_1);

	return (- ((freq_0 * log2(freq_0)) + (freq_1 * log2(freq_1)))) / length;
}

/**
 *
 */
static void parameter_quantization(uint8_t message[], int rx, int tx, int nr_packets, double csi_phase_data[][rx][tx][56])
{
    double params[5 * WINDOW_SIZE];
    double sliding_average;

    /* TODO: Calculate sliding window average of parameters for all antenna combinations. */
  	for (int i = 0; i < nr_packets / WINDOW_SIZE; ++i) {
        /* Caculate window average. */
        sliding_average = 0;
		for(int j = 0; j < WINDOW_SIZE; ++j) {
            estimate_csi_parameters(&csi_phase_data[i * WINDOW_SIZE + j][0][0][0], &params[j * 5]);
            sliding_average += params[j * 5];
		}
        sliding_average /= WINDOW_SIZE;

        /* Quantize. */
    	for(int j = 0; j < WINDOW_SIZE; ++j) {
			message[i * WINDOW_SIZE + j] = (params[j * 5] > sliding_average) ? 1 : 0;
    	}
	}
}

void csi_quantization(uint8_t message[], int rx, int tx, int nr_packets, double csi_phase_data[][rx][tx][56])
{
	double sliding_average;
	/* TODO: Calculate sliding window average of parameters for all antenna combinations. */
	for (int i = 0; i < nr_packets / WINDOW_SIZE; ++i) {
		/* Caculate window average. */
		sliding_average = 0;
		for(int j = 0; j < WINDOW_SIZE; ++j) {
			sliding_average += csi_phase_data[i * WINDOW_SIZE + j][0][0][0];
		}
		sliding_average /= WINDOW_SIZE;
		/* Quantize. */
		for(int j = 0; j < WINDOW_SIZE; ++j) {
			message[i * WINDOW_SIZE + j] = (csi_phase_data[j][0][0][0] > sliding_average) ? 1 : 0;
		}
	}
}
/**
 *
 * TODO: Use bits instead of bytes to use AND instead of multiplication for better space and time efficiency.
 */
struct key_recovery_data* encode(uint8_t message[], unsigned int msg_len, struct bch_control *bch)
{
    unsigned int nr_block_bits = calculate_nr_block_bits(bch);
    unsigned int nr_blocks = calculate_nr_blocks(bch);
    uint16_t vector_length = bch->n - bch->ecc_bits;

    /* Lengths for key recovery data. */
    uint16_t sketch_len = vector_length * nr_blocks;
    uint16_t ecc_len = bch->ecc_bits * nr_blocks;

    /* Create key recovery data struct on heap. */
    uint8_t *key_r_data = malloc(sizeof(struct key_recovery_data) + 2 * sketch_len + ecc_len);
    struct key_recovery_data *key_r = (struct key_recovery_data *) key_r_data;
    key_r->sketch_len = sketch_len;
    key_r->rand_x_len = sketch_len;
    key_r->ecc_len = ecc_len;
    uint8_t *sketch_s = key_r_data + sizeof(struct key_recovery_data);
    uint8_t *rand_x = sketch_s + sketch_len;
    uint8_t *ecc = rand_x + sketch_len;

    uint8_t rand_k[vector_length];
    uint8_t bit_mul[vector_length];
    uint8_t data_block[vector_length];

    assert(nr_block_bits <= vector_length);

    /* Add zero padding to end of data block, so it can be XOR'd with the sketch block. */
    memset(data_block + nr_block_bits, 0, vector_length - nr_block_bits);

    /* TODO: Check if this is really necessary. */
    memset(sketch_s, 0, sketch_len);

    generate_random_vector(rand_x, sketch_len);

	for (unsigned int i = 0; i < nr_blocks; ++i) {
		/* Multiplication of random vectors x and k will reduce randomness. */
        /* TODO: No function calls to avoid iterating twice. */
        generate_random_vector(rand_k, vector_length);
        generate_mul_vector(vector_length, rand_x + i * (vector_length), rand_k, bit_mul);

		/* Encode bit_mul */
        encodebits_bch(bch, bit_mul, ecc + i * (bch->ecc_bits));

        /* Copy block of data from message to data_block with zero padding. */
        memcpy(data_block, message + i * nr_block_bits, nr_block_bits);

        /* XOR data blocks with bit_mul to produce sketch_s block */
        generate_xor_vector(vector_length, data_block, bit_mul, sketch_s + i * (vector_length));
	}

    return (struct key_recovery_data*) key_r_data;
}

/**
 *
 * TODO: Use bits instead of bytes to use AND instead of multiplication for better space and time efficiency.
 */
void decode(uint8_t message[], unsigned int msg_len, struct bch_control *bch, struct key_recovery_data *key_r_data)
{
    unsigned int nr_block_bits = calculate_nr_block_bits(bch);
    unsigned int nr_blocks = calculate_nr_blocks(bch);
    unsigned int nr_errors;
    unsigned int errloc[ERROR_CORR_CAP];
    uint16_t vector_length = bch->n - bch->ecc_bits;

    /* Access data for key recovery. */
    uint8_t *sketch_s = key_r_data->data;
    uint8_t *rand_x = sketch_s + key_r_data->sketch_len;
    uint8_t *ecc_a = rand_x + key_r_data->sketch_len;

    uint8_t data_block[vector_length];
    assert(nr_block_bits <= vector_length);
    /* Add zero padding to end of data block, so it can be XOR'd with the sketch block. */
    memset(data_block + nr_block_bits, 0, vector_length - nr_block_bits);

	uint8_t recov_data_block[vector_length];
    uint8_t r1_b[vector_length];
    uint8_t bit_mul_b[vector_length];
    uint8_t ecc_b[bch->ecc_bits];

	for (unsigned int i = 0; i < nr_blocks; ++i) {
        /* Copy block of data from message to data_block with zero padding. */
        memcpy(data_block, message + i * nr_block_bits, nr_block_bits);

		/* Generate r1 by XOR data_block and sketch_a_s_block. */
		generate_xor_vector(vector_length, data_block, sketch_s + i * vector_length, r1_b);

		/* BCH decode r1 with the recieved ecc_a */
		memset(errloc, 0, ERROR_CORR_CAP * sizeof(unsigned int));
		nr_errors = decodebits_bch(bch, r1_b, ecc_a + i * (bch->ecc_bits), errloc);
		printf("\nr_errors detected: %d\n", nr_errors);

		/* Correct  r1, it becomes r2, variable name DOES NOT change! */
		corrupt_data(errloc, r1_b, nr_errors); // After correction r1_b is as same as bit_mul

		/* BCH encode r2,  it becomes r3 , variable name DOES NOT change! */
		memset(ecc_b, 0, bch->ecc_bits);
		encodebits_bch(bch, r1_b, ecc_b);

		/* Multiplication of r3 and rand_a_x_block */
		generate_mul_vector(vector_length, r1_b, rand_x + i * vector_length, bit_mul_b);

		/* XOR sketch_a_s_block with bit_mul_b to genrate recov_data_block */
		generate_xor_vector(vector_length, sketch_s + i * vector_length, bit_mul_b, recov_data_block);

		/* Recovering data_block */
		for (unsigned int j =  (i * nr_block_bits) ; j < nr_block_bits * (i + 1); ++j) {
			message[j] = recov_data_block[j - (i * nr_block_bits)];
		}
	}
}

int main(void)
{	/* Common */
	time_t begin = time(NULL);  // to calculate total execution time;
	unsigned int nrPkt = 128;   // to produce 128 bits;
	unsigned int nr_block_bits;   // Nr of bits in each blocks;
	unsigned int nr_blocks;      // Nr of blocks;
    unsigned int generator = 0;
    struct bch_control *bch;
    struct key_recovery_data *key_r_data;

    #if (DEBUG_KEYGEN == 0)
	srand(time(NULL)); // Seed for the random number generator
    #else
    printf("RNG seed will always be 1 in DEBUG mode!\n");
    #endif

	/* common */
    {
        int len = 0;
        assert((POLY_DEGREE >= 5) && (POLY_DEGREE <= 15));
        if (len == 0) {
            len = 1 << (POLY_DEGREE - 4);
        }
        assert((ERROR_CORR_CAP > 0) && (ERROR_CORR_CAP <= MAX_ERRORS));
        assert((len > 0) && (8*len + POLY_DEGREE * ERROR_CORR_CAP <= (1 << POLY_DEGREE) - 1));
    }

	bch = init_bch(POLY_DEGREE, ERROR_CORR_CAP, generator);
	nr_block_bits = calculate_nr_block_bits(bch);
	nr_blocks = calculate_nr_blocks(bch);
	printf("nr_block_bits		: %d\n", nr_block_bits);
	printf ("nr_blocks		: %d\n", nr_blocks);

    printf("\n**** Quantization (Moving Window) ****\n");
	uint8_t *message_a = malloc(nrPkt * sizeof(uint8_t)); // quantized bits
	uint8_t *message_b = malloc(nrPkt * sizeof(uint8_t)); // Bob: for quantized bits
    // parameter_quantization(message_a, 3, 3, nrPkt, alice_phase);
    // parameter_quantization(message_b, 3, 3, nrPkt, bob_phase);

	csi_quantization(message_a, 3, 3, nrPkt, alice_phase);
	csi_quantization(message_b, 3, 3, nrPkt, bob_phase);

	#if (DEBUG_KEYGEN != 0)
    {
        printf("\nmessage_a = ");
        print_array(message_a, nrPkt);
        printf("\nmessage_b = ");
        print_array(message_b, nrPkt);

        /* Check bit mismatch between A,B after quantization (DEBUG)*/
        printf("\nBit mismatch message_a, message_b: ");
        int k = 0;
        for (unsigned int i = 0; i < nrPkt; ++i) {
            if (message_a[i] != message_b[i]) {
                printf("%d ", i);
                k++;
            }
        }
        printf("\nNr. of errors		:%d", k);
        printf("\n");
    }
	#endif

	printf("\n**** Node A encodes  ****\n");
    key_r_data = encode(message_a, nrPkt, bch);

    printf("\n\n/* Node A sends sketch_a_s, rand_a_x and ecc_a to Node B */\n\n");

	printf("\n**** Node B decodes  ****\n");
    decode(message_b, nrPkt, bch, key_r_data);

	/* Check if correctly reconstructed. */
    #if (DEBUG_KEYGEN != 0)
	{
        int k = 0;
        printf("\n");
        printf("\nBit mismatch between message_a and message_b after reconciliating all blocks : \n");
        for (unsigned int i = 0; i < nrPkt; ++i) {
            if (message_a[i] != message_b[i]) {
                printf("%d ", i);
                k++;
            }
        }
        if (k != 0 ){
            printf("\n0, Successfully reconliced and reconstructed at Node B\n");
        }
    }
    #endif

    /* Calculate Metric Entropy */
    {
        double entropy;
        printf("\n\n**** Metric entropy of quantized bits ****\n");
        entropy = calculate_metric_entropy(message_a, nrPkt);
        printf("Metric entropy of quantized bits: %f\n\n", entropy);
    }

    free(message_a);
    free(message_b);

	time_t end = time(NULL);
	printf("\nTotal time elapsed: %ld sec.", (end - begin));
	printf("\n");

	return 0;
}
