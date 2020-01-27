#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <unistd.h>
#include "lib/bch/bch_codec.c"
#include "secure_sketch.c"
#include "parameter_estimator.cpp"
#include "alice_bob_phase.c"

#define MAX_ERRORS 2048
static int bitflip[MAX_ERRORS];

int main()
{
	int m = 6;
	int t = 5;
	int len = 0;
	int nrPkt = 128.0; // to produce 128 bits it needs to be set to 128
	int c, i, j, k, l, nerrors, nerrors_b;
	int opt_decode = 0;
	int opt_cache_encode = 0;
	int niterations = 1;
	uint8_t  *data_a, *data_b, *rand_data_a_x, *rand_data_a_k,*message_a, *message_b, *recov_message_a, *org_message_b, *ecc_a, *ecc_b, *bit_mul_a, *bit_mul_b, *sketch_a_s,*sketch_b, *r1_b, *r2_b, *r3_b, *wclean_b;
	struct bch_control *bch;
	unsigned int *pattern = NULL;
	unsigned int tmax;
	unsigned int generator = 0;
	double *pe_data_a, *pe_data_b, *pe_sum;

	srand(time(NULL));
	// srand48(123);
	assert((m >= 5) && (m <= 15));
	if (len == 0) {
		len = 1 << (m-4);
	}
	assert((t > 0) && (t <= MAX_ERRORS));
	assert((len > 0) && (8*len+m*t <= (1 << m)-1));

	bch = init_bch(m, t, generator);
	assert(bch);

	message_a = malloc(bch->n+bch->ecc_bits); // Alice: for quantized bits
	assert(message_a);
	message_b = malloc(bch->n+bch->ecc_bits); // Bob: for quantized bits
	assert(message_b);

	recov_message_a = malloc(bch->n);
	assert(recov_message_a);
	// org_message_b = malloc(bch->n);
	// assert(org_message_b);


	rand_data_a_x = malloc(bch->n);
	assert(rand_data_a_x);
	rand_data_a_k = malloc(bch->n);
	assert(rand_data_a_k);

	data_a = malloc(bch->n+bch->ecc_bits); // for secure sketch bits
	assert(data_a);
	data_b = malloc(bch->n+bch->ecc_bits);
	assert(data_b);

	ecc_a = malloc(bch->ecc_bits);
	assert(ecc_a);
	ecc_b = malloc(bch->ecc_bits);
	assert(ecc_b);

	// bit_mul_a = malloc(bch->n+bch->ecc_bits);
	// assert(bit_mul_a);
	// bit_mul_b = malloc(bch->n+bch->ecc_bits);
	// assert(bit_mul_b);

	// sketch_a_s = malloc(bch->n+bch->ecc_bits);
	// assert(sketch_a_s);
	// sketch_b = malloc(bch->n+bch->ecc_bits);
	// assert(sketch_b);

	// r1_b = malloc(bch->n+bch->ecc_bits);
	// assert(r1_b);
	// r2_b = malloc(bch->n+bch->ecc_bits);
	// assert(r2_b);
	// r3_b = malloc(bch->n+bch->ecc_bits);
	// assert(r3_b);
	// wclean_b = malloc(bch->n+bch->ecc_bits);
	// assert(wclean_b);

	/* Quantization method 1: simple mean based thresholding */

	// pe_data_a = malloc(5*sizeof(double));// Alice:placeholder for the each of the 5 values of the parameter_estimator
	// pe_data_b = malloc(5*sizeof(double));//Bob: placeholder for the each of the 5 values of the parameter_estimator
	// pe_sum = malloc(2*sizeof(double)); // placeholder for threshold, 2 for 2 nodes.
    // for(i=0; i<nrPkt; i++){
	//
	// for(j=0; j<1; j++){ // only considering Rx0Tx0 at the moment
	//
	// 		estimate_csi_parameters(alice_phase[i][j][j], pe_data_a);
	// 		pe_sum[0] += pe_data_a[0];
	// 		// printf("a = %f\n ", pe_data_a[0]);
	// 		estimate_csi_parameters(bob_phase[i][j][j], pe_data_b);
	// 		pe_sum[1] += pe_data_b[0];
	// 		// printf("b = %f\n ", pe_data_b[0]);
	// 	}
	// }
	// 	printf("\nAlice threshold: %f",pe_sum[0]);
	// 	printf("\nBob threshold: %f\n",pe_sum[1]);

	// pe_sum[0] /= nrPkt;
	// pe_sum[1] /= nrPkt;

	// for(i=0; i<nrPkt; i++) //Finally quantize!
	// {
	// 	message_a[i] = (pe_data_a[i] > pe_sum[0]) ? 1 : 0;
	// 	message_b[i] = (pe_data_b[i] > pe_sum[1]) ? 1 : 0;
	// }
	// EOS

	/* Quantization method (2): three point moving mean based */

	int wSize = 8;
  	int wNr = 0;

  	pe_data_a = malloc(5 * wSize * sizeof(double));
  	pe_data_b = malloc(5 * wSize * sizeof(double));
	pe_sum = malloc(2*sizeof(double)); // placeholder for threshold, 2 for 2 nodes.

  	for (i = 0; i < nrPkt / wSize; i++){
		for(j = 0; j < wSize; j++){
        	for(k = 0; k < 1; k++){
          		estimate_csi_parameters(alice_phase[i * wSize + j][k][k], pe_data_a, j * 5);
          		pe_sum[0] += pe_data_a[j * 5];
	          	estimate_csi_parameters(bob_phase[i * wSize + j][k][k], pe_data_b, j * 5);
          		pe_sum[1] += pe_data_b[j * 5];
        	}
		}
    	pe_sum[0] /= wSize;
    	pe_sum[1] /= wSize;
    	for(j = 0; j < wSize; j++){
			message_a[i * wSize + j] = (pe_data_a[j * wSize] > pe_sum[0]) ? 1 : 0;
        	message_b[i * wSize + j] = (pe_data_b[j * wSize] > pe_sum[1]) ? 1 : 0;
    	}
      	pe_sum[0] = 0;
      	pe_sum[1] = 0;
	}

  	for (i = (nrPkt / wSize) * wSize; i < nrPkt; i++){
    	for (k = 0; k < 1; k++){
      		estimate_csi_parameters(alice_phase[i][k][k], pe_data_a);
      		estimate_csi_parameters(bob_phase[i][k][k], pe_data_b);
      		message_a[i] = (pe_data_a[0] > pe_sum[0]) ? 1 : 0;
      		message_b[i] = (pe_data_b[0] > pe_sum[1]) ? 1 : 0;
    	}
  	}  // EOS


	/* copy message_a to a new location recov_message_a */
	// memcpy(recov_message_a, message_a, bch->n);

	/* copy message_b to a new location org_message_b */
	// memcpy(org_message_b, message_b, bch->n);


	printf("\nquantized_a(message_a) 	= ");
    for (i = 0;	i < bch->n; i++) {
        printf("%d", message_a[i]);
    }

	printf("\nquantized_b(message_b) 	= ");
    for (i = 0;	i < bch->n; i++) {
        printf("%d", message_b[i]);
    }

	/* Check bit mismatch between A,B after quantization*/

	printf("\nBit mismatch message_a, message_b: \n");
	for (i = 0;	i < bch->n; i++){
		if (message_a[i] != message_b[i]){
			printf("%d, ", i);
		}
	}

	/* encode after quantization */
	for (i = 0; i < niterations; i++) {
		memset(message_a+bch->n, 0, bch->ecc_bits); /* since parity bits add at the end of source bits, initially setting those to 0 */
		encodebits_bch(bch, message_a, message_a+bch->n);
		// encodebits_bch(bch, rand_data_a_k, ecc_a);
		printf("\necc_a		= ");
			for (j = 0;	j < bch->ecc_bits; j++) {
			printf("%d", message_a[bch->n+j]);
			}
			printf("\n");
	}

	printf("bch->ecc: %d", bch->ecc_bits);
	/* decode after quantization */

	unsigned int errloc[t];
    memset(errloc, 0xff, t);
    nerrors = decodebits_bch(bch, message_b, message_a+bch->n, errloc);

    printf("\nNr. Errors in message_b: %d\nat bit-Position: ", nerrors);
    for(i = 0; i < nerrors; i++){
        printf("%d ", errloc[i]);
    }
    // printf("\n");

    /*  correcting errors using corrupt_data function */
	correctbits_bch(bch, message_b, errloc, nerrors);

	// unsigned int errloc[t];
    // memset(errloc, 0xff, t);
    // nerrors = decodebits_bch(bch, message_b, ecc_a, errloc);

	// correctbits_bch(bch, message_b, errloc, nerrors);
	// printf("data_b			= ");
    // for (i = 0;	i < bch->n; i++) {
    //     printf("%d", data_b[i]);
    // }
	/* Check if data_b has been successfully decoded*/
	printf("\nBit mismatch after decoding message: \n");
	for (i = 0;	i < bch->n; i++){
		if (message_a[i] != message_b[i]){
			printf("%d ", i);
		}
	}

	printf("\nmessage_b(after): ");
	for (i = 0;	i < bch->n; i++){
			printf("%d", message_b[i]);
		}
	printf("\n");
	printf("\nmessage_a: ");
	for (i = 0;	i < bch->n; i++){
			printf("%d", message_a[i]);
		}
	printf("\n");


	/* Generate Random  Vector as same length as message_a */
    // generate_random_vector(bch->n, rand_data_a_x);
	// // printf("\nrand_data_a_x		= ");
    // for (i = 0;	i < bch->n; i++) {
    //     // printf("%d", rand_data_a_x[i]);
    // }
	// printf("\n");



	/* XOR "quantized bits with the  random_vector  */
	// generate_xor_vector(bch->n, message_a, rand_data_a_x, data_a);
	// printf("\ndata_a			= ");
    // for (i = 0;	i < bch->n; i++) {
    //     printf("%d", data_a[i]);
    // }
	// printf("\n");



	/* Generate additional Random  Vector  */
	// generate_random_vector(bch->n, rand_data_a_k);
	// printf("\nrand_data_a_k(init)	= ");
    // for (i = 0;	i < bch->n; i++) {
    //     printf("%d", rand_data_a_k[i]);
    // }
	// printf("\n");

	/* Encode with BCH */
	// printf("\n bch->n = %d \n bch->ecc_bits = %d",bch->n,bch->ecc_bits);
	/* encode_a to generate ecc from bch encoder */
	// for (i = 0; i < niterations; i++) {
	// 	memset(rand_data_a_k+bch->n, 0, bch->ecc_bits);
    //     encodebits_bch(bch, rand_data_a_k, rand_data_a_k+bch->n);
    //     printf("\nbch(rand_data_a_k)	= ");
    //     for (j = 0;	j < bch->ecc_bits; j++) {
    //         printf("%d", rand_data_a_k[bch->n+j]);
    //     }
    //     printf("\n");
	// }

	// alternate test with  ecc output

		// for (i = 0; i < niterations; i++) {
		// memset(rand_data_a_k+bch->n, 0, bch->ecc_bits); /* since parity bits add at the end of source bits, initially setting those to 0 */
		// 	encodebits_bch(bch, rand_data_a_k, rand_data_a_k+bch->n);
		// 	// encodebits_bch(bch, rand_data_a_k, ecc_a);
		// 	printf("\nbch(rand_data_a_k)	= ");
		// 	for (j = 0;	j < bch->n+bch->ecc_bits; j++) {
		// 		printf("%d", rand_data_a_k[j]);// print together with ecc bits
		// 	}
		// 	printf("\n");
		// }

	// alternate test with  bch encoding after XOR the random vector with message bits

		// for (i = 0; i < niterations; i++) {
		// 	memset(ecc_a, 0, bch->ecc_bits); /* since parity bits add at the end of source bits, initially setting those to 0 */
		// 	encodebits_bch(bch, data_a, ecc_a);
		// 	// encodebits_bch(bch, rand_data_a_k, ecc_a);
		// 	// printf("\nbch(data_a)		= ");
		// 	for (j = 0;	j < bch->n; j++) {
		// 		// printf("%d", data_a[j]);
		// 	}
		// 	// printf("\n");
		// }

	// bitwise multiplication of x and r
	//memset(rand_data_a_x+bch->n, 0, bch->ecc_bits); /* making rand_data_a_x same length as rand_data_a_k (which is now with parity bits) */
	// printf("r dot x         	= ");
	// for (i = 0; i < bch->n ; i++){
	// 	bit_mul_a[i] = rand_data_a_x[i] * rand_data_a_k[i];
	// 	printf("%d", bit_mul_a[i]);
	// }

	// printf("\n");

	/* XOR "message bits with x dot r  */
	// generate_xor_vector(bch->n+bch->ecc_bits, message_a, bit_mul_a, sketch_a_s);

	// printf("\nsketch_a_s 		= ");
    // for (i = 0;	i < bch->n; i++) {
    //     printf("%d", sketch_a_s[i]);
    // }
	// printf("\n");



	/* decode */

	// function [sketch_a_s, rand_data_a_x, data_a] = secure_sketch_generate(message_a, k)
	// function  data_a = secure_sketch_reproduce(message_b, sketch_a_s, rand_data_a_x, k)

	/* XOR message_b with sketch_a_s  */

	/*Generate r1 by XOR message_b and sketch_a_s */
	// generate_xor_vector(bch->n+bch->ecc_bits, message_b,sketch_a_s , r1_b);

	// printf("\nr1			= ");
    // for (i = 0;	i < bch->n; i++) {
    //     printf("%d", r1_b[i]);
    // }
	// printf("\n");

    /* BCH decode r1 */
	// unsigned int errloc[t];
    // memset(errloc, 0xff, t);
    // nerrors = decodebits_bch(bch, r1_b, r1_b+bch->n, errloc);
    // printf("Errors detected: %d\nat bit-Position: ", nerrors);
    // for(i = 0; i < nerrors; i++){
    //     printf("			%d ", errloc[i]);
    // }
    // printf("\n");


	/*correction*/
    // corrupt_data(errloc, r1_b, nerrors);
	/*after BCH decoding r1 becomes r2*/
	// printf("Decoded r1(r2)	= ");
	// for (i = 0;	i < bch->n; i++) {
    //     printf("%d", r1_b[i]);
    // }
    // printf("\n");

	// test correction using different function

	/* correctbits_bch - correct error locations as found in decodebits_bch
 	* @bch,@databits,@errloc: same as a previous call to decodebits_bch
 	* @nerr: returned from decodebits_bch
 	*/
	// correctbits_bch(bch, r1_b, errloc, nerrors);
	// printf("Decoded r1(r2)	= ");
	// for (i = 0;	i < bch->n; i++) {
    //     printf("%d", r1_b[i]);
    // }
    // printf("\n");


	/* bch encode r2,  it becomes r3 */
	// for (i = 0; i < niterations; i++) {
	// 	memset(r1_b+bch->n, 0, bch->ecc_bits);
    //     encodebits_bch(bch, r1_b, r1_b+bch->n);
    //     printf("\nr3         		= ");
    //     for (j = 0;	j < bch->n; j++) {
    //         printf("%d", r1_b[j]);
    //     }
    //     printf("\n");
	// }



	/*Bitwise Multiplication of r3 and rand_data_a_x*/

	// printf("r3 dot x         	= ");
	// for (i = 0; i < bch->n ; i++){
	// 	bit_mul_b[i] = r1_b[i] * rand_data_a_x[i];
	// 	printf("%d", bit_mul_b[i]);

	// }
	// printf("\n");

	/* XOR sketch_a_s with (r3 dot x) to genrate wclean_b */
	// generate_xor_vector(bch->n, sketch_a_s, bit_mul_b , wclean_b);
	// printf("Wclean         		= ");
	// for (i = 0; i < bch->n ; i++){
	// 	printf("%d", wclean_b[i]);
	// }
	// printf("\n");

	// for visual comparison printing again data_a
	// printf("\ndata_a			= ");
    // for (i = 0;	i < bch->n; i++) {
    //     printf("%d", data_a[i]);
    // }
	// printf("\n");


	/* XOR wclean_b with rand_data_a_x to genrate data_b */
	// generate_xor_vector(bch->n, wclean_b, rand_data_a_x , data_b);
	// printf("date_b         		= ");
	// for (i = 0; i < bch->n ; i++){
	// 	printf("%d", data_b[i]);
	// }
	// printf("\n");

// test steps with final BCH correction
		// for (i = 0; i < niterations; i++) {
		// memset(data_a+bch->n, 0, bch->ecc_bits); /* since parity bits add at the end of source bits, initially setting those to 0 */
		// 	encodebits_bch(bch, data_a, data_a+bch->n);
		// 	// encodebits_bch(bch, rand_data_a_k, ecc_a);
		// 	printf("\nbch(data_a)		= ");
		// 	for (j = 0;	j < bch->n+bch->ecc_bits; j++) {
		// 		printf("%d", data_a[j]);
		// 	}
		// 	printf("\n");
		// }

    /* BCH decode  */
	// unsigned int errloc_f[t];
    // memset(errloc_f, 0xff, t);
    // nerrors = decodebits_bch(bch, data_b, r1_b+bch->n, errloc_f);
    // printf("Errors detected: %d\nat bit-Position: ", nerrors);
    // for(i = 0; i < nerrors; i++){
    //     printf("			%d ", errloc_f[i]);
    // }
    // printf("\n");


	// /*correction*/
    // corrupt_data(errloc_f, data_b, nerrors);
	// /*after BCH decoding r1 becomes r2*/
	// printf("Decoded 	= ");
	// for (i = 0;	i < bch->n+bch->ecc_bits; i++) {
    //     printf("%d", data_b[i]);
    // }
    // printf("\n");


	/*To reproduce R*/
	// generate_xor_vector(bch->n, message_a, sketch_a_s, data_b);



    /* detect errors */
	// unsigned int errloc[t];
    // memset(errloc, 0xff, t);
    // nerrors = decodebits_bch(bch, data_b, data_a+bch->n, errloc);

    // printf("Errors detected: %d\nat bit-Position: ", nerrors);
    // for(i = 0; i < nerrors; i++){
    //     printf("%d ", errloc[i]);
    // }
    // printf("\n");

    /* alternate correcting errors using corrupt_data function */

    // corrupt_data(errloc, data_b, nerrors);
    // for (i = 0;	i < bch->n; i++) {
    //     printf("%d", data_b[i]);
    // }

	/* test recovering message_b from the data_b*/

	// generate_xor_vector(bch->n, message_b, rand_data_a_x, data_b); // XOR quantized bits with the random vector of A

	// printf("\ndata_b			= ");
    // for (i = 0;	i < bch->n; i++) {
    //     printf("%d", data_b[i]);
    // }
	// printf("\n");


	/* Decode data_b using the ecc_a*/
	/*
	unsigned int errloc[t];
    memset(errloc, 0xff, t);
    nerrors = decodebits_bch(bch, data_b, ecc_a, errloc);

    printf("\nNr. Errors in data_b: %d\nat bit-Position: ", nerrors);
    for(i = 0; i < nerrors; i++){
        printf("%d ", errloc[i]);
    }
    // printf("\n");

    /*  correcting errors using corrupt_data function */
    // corrupt_data(errloc, data_b, nerrors);
	// printf("data_b			= ");
    // for (i = 0;	i < bch->n; i++) {
    //     printf("%d", data_b[i]);
    // }
	/* Check if data_b has been successfully decoded*/
	// printf("\nBit mismatch data_a, data_b: \n");
	// for (i = 0;	i < bch->n; i++){
	// 	if (data_a[i] != data_b[i]){
	// 		printf("%d ", i);
	// 	}
	// }

	/* undo XOR for data_b  */
		// generate_xor_vector(bch->n, data_b, rand_data_a_x, recov_message_a);
		// printf("\ndata_b XOR rand_data_a(reconciled A)= ");
    	// for (i = 0;	i < bch->n; i++) {
        // 	printf("%d", recov_message_a[i]);
    	// }
	/* check if message_a is successfully recovered*/
	// for (i = 0; i < bch->n; i++){
	// 	if (message_a[i] == recov_message_a[i]){
	// 		k = 0;
	// 	}
	// }

	// if (k = 0 ){
	// 	printf("\nSuccessfully recovered message_a at node B");
	// }
	// else{
	// 	printf("\nNot recovered\n");
	// }



	free(pattern);
	// free(data_a);
	// free(data_b);
	// free(rand_data_a_x);
	free_bch(bch);


	//return data_a, message_a, message_b;
	//return message_a, message_b;
	return 0;
}