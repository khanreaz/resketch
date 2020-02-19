#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "lib/bch/bch_codec.c"
#include "secure_sketch.c"
#include "parameter_estimator.cpp"
// #include "alice_bob_phase.c"
// #include "2019-07-10_16-03-24_csi_log.c"
// #include "2019-07-10_15-12-45_csi_log.c"
#include "2020-02-19_14-05-25_csi_log.c"

#define MAX_ERRORS 2048


int main()
{	/* Common */
	// time_t begin = time(NULL); // to calculate total execution time;
	int m = 7; // BCH polynomial degree;
	int t = 9; // erroc correction capability;
	int len = 0;
	int nrPkt = 128; // to produce 128 bits;
	int i = 0, j = 0, k = 0, nerrors = 0;
	int nrBlockBits = 0;// Nr of bits in each blocks;
	int nrBlocks = 0; // Nr of blocks;
	double *pe_sum;
	struct bch_control *bch;
	unsigned int *pattern = NULL, generator = 0, errloc[t];

	/* for Node A */
	uint8_t *message_a, *data_a_block, *rand_a_x, *rand_a_x_block, *rand_a_k, *rand_a_k_block, *ecc_a, *ecc_a_block,  *bit_mul_a, *sketch_a_R, *sketch_a_s_block, *sketch_a_s;
	double *pe_data_a;

	/* for Node B*/
	uint8_t *message_b, *data_b_block, *ecc_b, *recov_sketch_R_block, *bit_mul_b, *r1_b, *recov_data_a_block;
	double *pe_data_b;


	// srand(time(NULL)); // Seed for the random number generator

	/* common */
	assert((m >= 5) && (m <= 15));
	if (len == 0) {
		len = 1 << (m-4);
	}
	assert((t > 0) && (t <= MAX_ERRORS));
	assert((len > 0) && (8*len+m*t <= (1 << m)-1));

	bch = init_bch(m, t, generator);
	assert(bch);


	nrBlockBits = (bch->n + 1) / (2 * t - 2);
	nrBlocks = (bch->n + 1) / nrBlockBits;
	printf("\nnrBlockBits		: %d", nrBlockBits);
	printf ("\nnrBlocks		: %d", nrBlocks);


	/* for Node A */
	message_a = malloc( nrPkt * sizeof(uint8_t));; // quantized bits
	assert(message_a);

	data_a_block = malloc( nrPkt * sizeof(uint8_t));
	assert(data_a_block);

	rand_a_x = malloc((bch->n - bch->ecc_bits) * nrBlocks * sizeof(uint8_t));
	assert(rand_a_x);

	rand_a_x_block = malloc((bch->n - bch->ecc_bits) * sizeof(uint8_t));
	assert(rand_a_x_block);

	rand_a_k = malloc((bch->n - bch->ecc_bits) * sizeof(uint8_t));
	assert(rand_a_k);

	rand_a_k_block = malloc((bch->n - bch->ecc_bits) * sizeof(uint8_t));
	assert(rand_a_k_block);

	ecc_a = malloc(bch->ecc_bits * nrBlocks * sizeof(uint8_t));
	assert(ecc_a);

	ecc_a_block = malloc(bch->ecc_bits * sizeof(uint8_t));
	assert(ecc_a_block);

	sketch_a_R = malloc((bch->n - bch->ecc_bits) * sizeof(uint8_t));
	assert(sketch_a_R);

	sketch_a_s = malloc((bch->n - bch->ecc_bits) * nrBlocks * sizeof(uint8_t));
	assert(sketch_a_s);

	sketch_a_s_block = malloc((bch->n - bch->ecc_bits) * sizeof(uint8_t));
	assert(sketch_a_s_block);

	bit_mul_a = malloc((bch->n - bch->ecc_bits) * sizeof(uint8_t));
	assert(bit_mul_a);



	/* for Node B */
	message_b = malloc(nrPkt * sizeof(uint8_t)); // Bob: for quantized bits
	assert(message_b);

	data_b_block = malloc(nrPkt * sizeof(uint8_t));
	assert(data_b_block);

	r1_b = malloc((bch->n - bch->ecc_bits) * sizeof(uint8_t));
	assert(r1_b);

	bit_mul_b =  malloc((bch->n - bch->ecc_bits) * sizeof(uint8_t));
	assert(bit_mul_b);

	ecc_b = malloc((bch->n - bch->ecc_bits) * sizeof(uint8_t));
	assert(ecc_b);

	recov_sketch_R_block = malloc((bch->n - bch->ecc_bits) * sizeof(uint8_t));
	assert(recov_sketch_R_block);

	recov_data_a_block = malloc((bch->n - bch->ecc_bits) * sizeof(uint8_t));
	assert(data_b_block);

	/* EOVariables*/

	/* Quantization method (1): simple mean based thresholding */

	// pe_data_a = malloc(5*sizeof(double));// Alice:placeholder for the each of the 5 values of the parameter_estimator
	// pe_data_b = malloc(5*sizeof(double));//Bob: placeholder for the each of the 5 values of the parameter_estimator
	// pe_sum = malloc(2*sizeof(double)); // placeholder for threshold, 2 for 2 nodes.
    // for(i=0; i<nrPkt; i++){

	// for(j=0; j<1; j++){ // only considering Rx0Tx0 at the moment

	// 		estimate_csi_parameters(alice_phase[i][j][j], pe_data_a);
	// 		pe_sum[0] += pe_data_a[0];
	// 		printf("a = %f\n ", pe_data_a[0]);
	// 		estimate_csi_parameters(bob_phase[i][j][j], pe_data_b);
	// 		pe_sum[1] += pe_data_b[0];
	// 		printf("b = %f\n ", pe_data_b[0]);
	// 	}
	// }


	// pe_sum[0] /= nrPkt;
	// pe_sum[1] /= nrPkt;

	// 	printf("\nAlice threshold: %f",pe_sum[0]);
	// 	printf("\nBob threshold	: %f\n",pe_sum[1]);

	// for(i=0; i<nrPkt; i++) //Finally quantize!
	// {
	// 	message_a[i] = (pe_data_a[i] > pe_sum[0]) ? 1 : 0;
	// 	message_b[i] = (pe_data_b[i] > pe_sum[1]) ? 1 : 0;
	// }
	/* EOQunatization (1) */

	/* Quantization method (2): "x" point moving window  based */
	printf("\n**** Quantization (Moving Window) ****\n");

	int wSize = 3; // Set window size, x: nrPkt > wSize > 1.

  	pe_data_a = malloc(5 * wSize * sizeof(double));
  	pe_data_b = malloc(5 * wSize * sizeof(double));
	pe_sum = malloc(2*sizeof(double)); // placeholder for threshold, 2 for 2 nodes.

  	for (i = 0; i < nrPkt / wSize; i++){
		for(j = 0; j < wSize; j++){
        	for(k = 0; k < 1; k++){
          		estimate_csi_parameters(alice_phase[i * wSize + j][k][k], pe_data_a, j * 5); // For Node A
          		pe_sum[0] += pe_data_a[j * 5];
	          	estimate_csi_parameters(bob_phase[i * wSize + j][k][k], pe_data_b, j * 5); // for Node B
          		pe_sum[1] += pe_data_b[j * 5];
        	}
		}
    	pe_sum[0] /= wSize;
    	pe_sum[1] /= wSize;
    	for(j = 0; j < wSize; j++){
			message_a[i * wSize + j] = (pe_data_a[j * 5] > pe_sum[0]) ? 1 : 0;
        	message_b[i * wSize + j] = (pe_data_b[j * 5] > pe_sum[1]) ? 1 : 0;
    	}
      	pe_sum[0] = 0;
      	pe_sum[1] = 0;
	}

	free(pe_data_a);
	free(pe_data_b);
	free(pe_sum);
	// EOQunatization(2)

	/* copy message_a to a new location recov_message_a */
	// memcpy(recov_message_a, message_a, bch->n);

	/* copy message_b to a new location org_message_b */
	// memcpy(org_message_b, message_b, bch->n);

	/* DEBUG Info */
	printf("\nmessage_a = ");
    for (i = 0;	i < nrPkt; i++) {
        printf("%d", message_a[i]);
    }

	printf("\nmessage_b = ");
    for (i = 0;	i < nrPkt; i++) {
        printf("%d", message_b[i]);
    }

	/* Check bit mismatch between A,B after quantization (DEBUG)*/
	k = 0;
	printf("\nBit mismatch message_a, message_b: ");
	for (i = 0;	i < nrPkt; i++){
		if (message_a[i] != message_b[i]){
			printf("%d ", i);
			k++;

		}
	}
	printf("\nNr. of errors		:%d", k-1);
	printf("\n");
	/* EODebug info */

	printf("\n**** Node A encodes  ****\n");


	/* Generate  Random  Vector as same length as data block */
	generate_random_vector((bch->n - bch->ecc_bits) * nrBlocks  , rand_a_x);
	// printf("\nrand_a_x		= ");
	// for (j = 0;	j < (bch->n - bch->ecc_bits) * nrBlocks; j++) {
	// 	printf("%d", rand_a_x[j]);
	// }
	// printf("\n");

	/* Generate second Random  Vector as same length as data block */
	generate_random_vector((bch->n - bch->ecc_bits) * nrBlocks  , rand_a_k);
	// printf("\nrand_a_k		= ");
	// for (j = 0;	j < (bch->n - bch->ecc_bits) * nrBlocks; j++) {
	// 	printf("%d", rand_a_k[j]);
	// }
	// printf("\n");

	/* Check bit mismatch between rand_a_x, rand_a_k  (DEBUG)*/
	k = 0;
	// printf("\nBit mismatch rand_a_x, rand_a_k: ");
	for (i = 0;	i < (bch->n - bch->ecc_bits) * nrBlocks; i++){
		if (rand_a_x[i] != rand_a_k[i]){
			// printf("%d ", i);
			k++;

		}
	}
	// printf("\nNr. of  rand diff	:%d", k-1);
	// printf("\n");

	/* create  blocks from message */
	memset(data_a_block, 0, (bch->n - bch->ecc_bits)); // set all block bits to 0 at first;

	for (i = 0; i < nrBlocks; i++){
		for (j = i * nrBlockBits; j < (i + 1) * nrBlockBits; j++){
			data_a_block[j - i * nrBlockBits] = message_a[j];
		}
		// for (j = 0;	j < bch->n - bch->ecc_bits; j++) {
		// 	printf("%d", data_a_block[j]);
		// }
		// printf("\n");

		/* create blocks from the rand_a_x  */
		for (j = i * (bch->n - bch->ecc_bits); j < (i + 1) * (bch->n - bch->ecc_bits); j++){
			rand_a_x_block[j - i * (bch->n - bch->ecc_bits)] = rand_a_x[j];
		}

		// printf("\nrand_a_x_block		= ");
		// for (j = 0;	j < bch->n - bch->ecc_bits; j++) {
		// 	printf("%d", rand_a_x_block[j]);
		// }
		// printf("\n");

		/* XOR data blocks with the  first random_vector rand_a_x  */
		generate_xor_vector(bch->n - bch->ecc_bits, data_a_block, rand_a_x_block, sketch_a_R);
		// printf("\nsketch_a_R		= ");
		// for (j = 0;	j < bch->n - bch->ecc_bits; j++) {
		// 	printf("%d", sketch_a_R[j]);
		// }
		// printf("\n");

		/* create blocks from the second random vector rand_a_k  */
		for (j = i * (bch->n - bch->ecc_bits); j < (i + 1) * (bch->n - bch->ecc_bits); j++){
			rand_a_k_block[j - i * (bch->n - bch->ecc_bits)] = rand_a_k[j];
		}

		// printf("\nrand_a_k_block		= ");
		// for (j = 0;	j < bch->n - bch->ecc_bits; j++) {
		// 	printf("%d", rand_a_k_block[j]);
		// }
		// printf("\n");

		/*  multiplication of x and k */
		// printf("\nbit_mul_a 		= ");
		generate_mul_vector(bch->n - bch->ecc_bits, rand_a_x_block, rand_a_k_block, bit_mul_a);
		// for (j = 0;	j < bch->n - bch->ecc_bits; j++) {
		// 	printf("%d", bit_mul_a[j]);
		// }
		// printf("\n");

		/* encoding bit_mul_a  */

		encodebits_bch(bch, bit_mul_a, ecc_a_block);
		// printf("\necc_a_block 		= "); // for DEBUG
		// for (j = 0; j < bch->ecc_bits ; j++){
		// 	printf("%d", ecc_a_block[j]);
		// }
		// printf("\n");

		/* create ecc_a by concatenating all ecc_a_block */
		for (j = i * (bch->ecc_bits); j < (bch->ecc_bits) * (i + 1); j++){
			ecc_a[j] = ecc_a_block[j - i * (bch->ecc_bits)];
			// printf("%d", ecc_a[j]);
		}

		/* XOR data blocks with encoded bit_mul_a to produce sketch_a_s_block  */
		generate_xor_vector(bch->n - bch->ecc_bits, data_a_block, bit_mul_a, sketch_a_s_block);
		// printf("\nsketch_a_s_block		= "); // for DEBUG
		// for (j = 0; j < bch->n - bch->ecc_bits ; j++){
		// 	printf("%d", sketch_a_s_block[j]);
		// }
		// printf("\n");

		/* create sketch_a_s by concatenating all sketch_a_s_block */
		for (j = i * (bch->n - bch->ecc_bits); j < (bch->n - bch->ecc_bits) * (i + 1); j++){
			sketch_a_s[j] = sketch_a_s_block[j - i * (bch->n - bch->ecc_bits)];
			printf("%d", sketch_a_s[j]);
		}
		printf("\n");
	}
	// EORobustScheme


	printf("\n\n/* Node A sends sketch_a_s, rand_a_x and ecc_a to Node B */\n\n");

	printf("\n**** Node B decodes  ****\n");


	/* create  blocks from message */

	memset(data_b_block, 0, (bch->n - bch->ecc_bits)); // set all block bits to 0 at first;
	memset(sketch_a_s_block, 0, (bch->n - bch->ecc_bits));

	for (i = 0; i < nrBlocks; i++){
		for (j = i * nrBlockBits; j < (i + 1) * nrBlockBits; j++){
			data_b_block[j - i * nrBlockBits] = message_b[j];
		}

		/* create blocks from the received sketch_a_s to process it with each data_b_block */
		for (j = i * (bch->n - bch->ecc_bits); j < (i + 1) * (bch->n - bch->ecc_bits); j++){
				sketch_a_s_block[j - i * (bch->n - bch->ecc_bits)] = sketch_a_s[j];
		}

		// for (j = 0; j < bch->n - bch->ecc_bits ; j++){
		// printf("%d", sketch_a_s_block[j]);
		// }
		// printf("\n");

		/*Generate r1 by XOR data_b_block and sketch_a_s_block */
		generate_xor_vector(bch->n - bch->ecc_bits, data_b_block, sketch_a_s_block , r1_b);
		// for (j = 0; j< (bch->n - bch->ecc_bits); j++){ // Print for DEBUG
			// printf("%d", r1_b[j]);
		// }
		// printf("\n");

		/* BCH decode r1 with the recieved ecc_a */
		/* create blocks from the received ecc_a  */
		for (j = i * (bch->ecc_bits); j < (i + 1) * (bch->ecc_bits); j++){
			ecc_a_block[j - i * (bch->ecc_bits)] = ecc_a[j];
		}

		memset(errloc, 0, t);
		nerrors = decodebits_bch(bch, r1_b, ecc_a_block, errloc);
		// printf("\nErrors detected: %d", nerrors); // for DEBUG;

		/* Correct  r1, it becomes r2, variable name DOES NOT change! */
		corrupt_data(errloc, r1_b, nerrors); // After correction r1_b is as same as bit_mul_a
		for (j = 0; j< (bch->n - bch->ecc_bits); j++){ // Print for DEBUG
			// printf("%d", r1_b[j]);
		}
		// printf("\n");

		/* BCH encode r2,  it becomes r3 , variable name DOES NOT change! */
		memset(ecc_b, 0, bch->ecc_bits);
		encodebits_bch(bch, r1_b, ecc_b);

		// for (j = 0; j< (bch->ecc_bits); j++){ // Print for DEBUG
		// 	printf("%d", ecc_b[j]);
		// }
		// printf("\n");

		/* create blocks from the received rand_a_x  */
		for (j = i * (bch->n - bch->ecc_bits); j < (i + 1) * (bch->n - bch->ecc_bits); j++){
			rand_a_x_block[j - i * (bch->n - bch->ecc_bits)] = rand_a_x[j];
		}
		// for (j = 0; j< (bch->n - bch->ecc_bits); j++){ // Print for DEBUG
		// 	printf("%d", rand_a_x_block[j]);
		// }
		// printf("\n");

		/* Multiplication of r3 and rand_a_x_block */
		generate_mul_vector(bch->n - bch->ecc_bits, r1_b, rand_a_x_block, bit_mul_b);

		/* XOR sketch_a_s_block with bit_mul_b to genrate recov_data_a_block */
		generate_xor_vector(bch->n - bch->ecc_bits, sketch_a_s_block, bit_mul_b , recov_data_a_block);

		/* Recovering data_a_block */
		// printf("\n");
		for (j =  (i * nrBlockBits) ; j < nrBlockBits * (i + 1); j++) {
			message_b[j] = recov_data_a_block[j - (i * nrBlockBits)];
			// printf("%d", message_b[j]);
		}
	}
	/*EOReconciliaition*/



	/* Check if correctly reconstructed (DEBUG) */
	// printf("\nmessage_a 	= ");
	printf("\n");
	// for (i = 0; i < nrPkt ; i++) {
		// printf ("%d" , message_a[i]);
	// }

	printf("\nBit mismatch between message_a and message_b after reconciliating all blocks : \n");
	for (i = 0;	i < nrPkt; i++){
		if (message_a[i] != message_b[i]){
			printf("%d ", i);
			k =0;
		}
	}
	if (k != 0 ){
		printf("\n0, Successfully reconliced and reconstructed at Node B\n");
	}





	/* Calculate Metric Entropy */
	// printf("\n\n**** Metric entropy of quantized bits ****\n");
	 double entropy_H;
	// H(X) = [-[(freq_0)log2(freq_0))+((freq_1)log2(freq_1))]] / number of bits ;


	/* Calculate frequency of each alphabet */
	double freq_0 = 0.0;
	double freq_1 = 0.0;
	double c = 0.0;

	for (int i = 0; i < nrPkt; i++) {
		if (message_a[i] == 0)
			++c;
	}
	freq_0 = c / nrPkt;
	printf("\nfreq_0: %f", freq_0);

	c =0.0;
	for (int i = 0; i < nrPkt; i++) {
		if (message_a[i] == 1)
			++c;
	}
	freq_1 = c / nrPkt;
	printf("\nfreq_1: %f", freq_1);

	entropy_H = (- ((freq_0 * log2(freq_0)) + (freq_1 * log2(freq_1)))) / nrPkt;
	printf("\nMetric entropy of quantized bits: %f", entropy_H);
	printf("\n\n");

	// printf("\n\n**** Metric entropy of PRNG bits ****\n");
	double entropy_PRNG;
	// H(X) = [-[(freq_0)log2(freq_0))+((freq_1)log2(freq_1))]] / number of bits ;


	/* Calculate frequency of each alphabet */
	freq_0 = 0.0;
	freq_1 = 0.0;
	c = 0.0;

	for (int i = 0; i < (nrBlockBits * nrBlocks); i++) {
		if (rand_a_x[i] == 0)
			++c;
	}
	freq_0 = c / (nrBlockBits * nrBlocks);
	printf("\nfreq_0: %f", freq_0);

	c =0.0;
	for (int i = 0; i < (nrBlockBits * nrBlocks); i++) {
		if (rand_a_x[i] == 1)
			++c;
	}
	freq_1 = c / (nrBlockBits * nrBlocks);
	printf("\nfreq_1: %f", freq_1);

	entropy_PRNG = (- ((freq_0 * log2(freq_0)) + (freq_1 * log2(freq_1)))) / nrPkt;
	printf("\nMetric entropy of PRNG bits: %f", entropy_PRNG);
	printf("\n\n");



	/* convert to HEX */
	// const char *a = "This is Message";
	// char *hex;

	// hex = bin2hex((unsigned char *)a, strlen(a));





	// // strrev(hex);
	// printf("HEX: %s", hex );



	/* Common */
	free(pattern);
	free_bch(bch);

	/* Node A */

	free(data_a_block);
	free(rand_a_k);
	// free(rand_a_k_block);
	free(rand_a_x);
	free(rand_a_x_block);

	free(sketch_a_s);
	free(sketch_a_s_block);
	free(bit_mul_a);


	/* Node B */

	free(data_b_block);
	free(r1_b);
	free(bit_mul_b);
	free(recov_data_a_block);
	free(recov_sketch_R_block);
	free(ecc_b);


	// free(rand_a_x); // Node B also has it so needs to be freed!



	// time_t end = time(NULL);
	// printf("\nTotal time elapsed: %d sec.", (end - begin));
	// printf("\n");


	//return message_a;
	//return message_b;
	return 0;
}