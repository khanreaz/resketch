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
// #include "2019-07-10_16-03-24_csi_log.c"
// #include "2019-07-10_15-12-45_csi_log.c"

#define MAX_ERRORS 2048
static int bitflip[MAX_ERRORS];

int main()
{
	time_t begin = time(NULL);
	int m = 7;
	int t = 9;
	int len = 0;
	int nrPkt = 128; // to produce 128 bits it needs to be set to 128
	int i, j, k, l, nerrors, niterations = 1; //c, nerrors_b,opt_cache_encode = 0, opt_decode = 0;
	uint8_t  *data_a, *data_a_block0, *data_a_block1, *data_b, *data_b_block0, *data_b_block1,*rand_data_a_x, *message_a, *message_b, *recov_message_a, *ecc_b, *ecc_a_block0, *ecc_a_block1;
	uint8_t *rand_data_a_k, *bit_mul_a, *bit_mul_b, *sketch_a_R, *sketch_a_r, *sketch_a_s, *sketch_b, *r1_b, *r2_b, *r3_b, *wclean_b; // for secresketch
	struct bch_control *bch;
	unsigned int *pattern = NULL, generator = 0, tmax;
	double *pe_data_a, *pe_data_b, *pe_sum;

	srand(time(NULL)); // Seed for the random number generator. Both node shall use the 2nd parameter: TOF after agreeing to a common TOF

	assert((m >= 5) && (m <= 15));
	if (len == 0) {
		len = 1 << (m-4);
	}
	assert((t > 0) && (t <= MAX_ERRORS));
	assert((len > 0) && (8*len+m*t <= (1 << m)-1));

	bch = init_bch(m, t, generator);
	assert(bch);

	message_a = malloc(nrPkt); // Alice: for quantized bits
	assert(message_a);
	message_b = malloc(nrPkt); // Bob: for quantized bits
	assert(message_b);

	recov_message_a = malloc(nrPkt);
	assert(recov_message_a);
	// org_message_b = malloc(bch->n);
	// assert(org_message_b);


	rand_data_a_x = malloc(nrPkt);
	assert(rand_data_a_x);
	rand_data_a_k = malloc(nrPkt);
	assert(rand_data_a_k);

	data_a = malloc(nrPkt);
	assert(data_a);
	data_a_block0 = malloc(bch->n - bch->ecc_bits);
	assert(data_a_block0);
	data_a_block1 = malloc(bch->n - bch->ecc_bits);
	assert(data_a_block1);

	data_b = malloc(nrPkt);
	assert(data_b);
	data_b_block0 = malloc(bch->n - bch->ecc_bits);
	assert(data_b_block0);
	data_b_block1 = malloc(bch->n - bch->ecc_bits);
	assert(data_b_block1);

	// ecc_a = malloc(bch->ecc_bits);
	// assert(ecc_a);
	ecc_a_block0 = malloc(bch->ecc_bits);
	assert(ecc_a_block0);
	ecc_a_block1 = malloc(bch->ecc_bits);
	assert(ecc_a_block1);

	ecc_b = malloc(bch->ecc_bits);
	assert(ecc_b);

	bit_mul_a = malloc(bch->n+bch->ecc_bits);
	assert(bit_mul_a);
	bit_mul_b = malloc(bch->n+bch->ecc_bits);
	assert(bit_mul_b);

	sketch_a_R = malloc(bch->n+bch->ecc_bits);
	assert(sketch_a_R);
	sketch_a_r = malloc(bch->n+bch->ecc_bits);
	assert(sketch_a_r);
	sketch_a_s = malloc(bch->n+bch->ecc_bits);
	assert(sketch_a_s);
	sketch_b = malloc(bch->n+bch->ecc_bits);
	assert(sketch_b);

	r1_b = malloc(bch->n+bch->ecc_bits);
	assert(r1_b);
	r2_b = malloc(bch->n+bch->ecc_bits);
	assert(r2_b);
	r3_b = malloc(bch->n+bch->ecc_bits);
	assert(r3_b);
	wclean_b = malloc(bch->n+bch->ecc_bits);
	assert(wclean_b);

	/* Quantization method (1): simple mean based thresholding */

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

	/* Quantization method (2): X point moving window  based */

	int wSize = 3; // Set window size: nrPkt > wSize > 1.
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
  	}  // EOQunatization(2)


	/* copy message_a to a new location recov_message_a */
	// memcpy(recov_message_a, message_a, bch->n);

	/* copy message_b to a new location org_message_b */
	// memcpy(org_message_b, message_b, bch->n);


	printf("\nquantized_a(message_a) 	= ");
    for (i = 0;	i < nrPkt; i++) {
        printf("%d", message_a[i]);
    }

	printf("\nquantized_b(message_b) 	= ");
    for (i = 0;	i < nrPkt; i++) {
        printf("%d", message_b[i]);
    }

	/* Check bit mismatch between A,B after quantization*/

	printf("\nBit mismatch message_a, message_b: \n");
	for (i = 0;	i < nrPkt; i++){
		if (message_a[i] != message_b[i]){
			printf("%d ", i);
		}
	}


			/* Create blocks from data_a*/
			/* 1st block length: bch->n - bch->ecc_bits */
			printf("\ndata_a_block0		= ");
			for (i =0; i < bch->n - bch->ecc_bits; i++) {
				data_a_block0[i] = message_a[i];
				printf("%d", data_a_block0[i]);
			}
			/* 2nd block length: nrPkt - (bch->n - bch->ecc_bits) */
			printf("\ndata_a_block1		= ");
			for (i =0; i < nrPkt - (bch->n - bch->ecc_bits); i++) {
				data_a_block1[i] = message_a[(bch->n - bch->ecc_bits)+i];
				printf("%d", data_a_block1[i]);
			}
				/* Padding for the 2nd block to make length: bch->n - bch->ecc_bits */
				int nrPadding =0;
				nrPadding = ((bch->n - bch->ecc_bits)  - (nrPkt - (bch->n - bch->ecc_bits)));
				// printf("\nnrPadding: %d", nrPadding);
				memset(data_a_block1 + nrPkt - (bch->n - bch->ecc_bits), 0, nrPadding);
				printf("\ndata_a_block1(padded)	= ");
				for (i =0; i < 71; i++) {
						printf("%d", data_a_block1[i]);
				}
				/* Generate Random  Vector as same length as message_a */
				generate_random_vector(bch->n - bch->ecc_bits , rand_data_a_x);
				printf("\nrand_data_a_x		= ");
				for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
					printf("%d", rand_data_a_x[i]);
				}
				// printf("\n");

					/* XOR "quantized bits with the  random_vector  */
					generate_xor_vector(bch->n - bch->ecc_bits, data_a_block0, rand_data_a_x, sketch_a_R);
					printf("\nsketch_a_R		= ");
					for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
						printf("%d", sketch_a_R[i]);
					}
					// printf("\n");




					/* */
						/* Generate additional Random  Vector  */
 						generate_random_vector(bch->n - bch->ecc_bits, rand_data_a_k);
 						memset(ecc_a_block0, 0, bch->ecc_bits);// memset(ecc_a_block0, 0, bch->ecc_bits); /* since parity bits add at the end of source bits, initially setting those to 0 */
						encodebits_bch(bch, rand_data_a_k, ecc_a_block0);

						printf("\nrand_data_a_k(r)	= ");
 						for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
 						    printf("%d", rand_data_a_k[i]);
 						}
 						printf("\n");


 						/* bitwise multiplication of x and r */
 						// memset(rand_data_a_x+bch->n, 0, bch->ecc_bits); /* making rand_data_a_x same length as rand_data_a_k  */
 						printf("r dot x         	= ");
 						for (i = 0; i < bch->n - bch->ecc_bits ; i++){
 							bit_mul_a[i] = rand_data_a_x[i] * rand_data_a_k[i];
 							printf("%d", bit_mul_a[i]);
 						}
 						printf("\n");
						/* test with encoding bit_mul_a instead of rand_data_a_k */
						memset(ecc_a_block1, 0, bch->ecc_bits);// memset(ecc_a_block0, 0, bch->ecc_bits); /* since parity bits add at the end of source bits, initially setting those to 0 */
						encodebits_bch(bch, bit_mul_a, ecc_a_block1);

 						/* XOR "message bits with x dot r  */
 						generate_xor_vector(bch->n - bch->ecc_bits, data_a_block0, bit_mul_a, sketch_a_s);

 						printf("\nsketch_a_s 		= ");
 						for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
 						    printf("%d", sketch_a_s[i]);
 						}
 						printf("\n");
					//EOEncoding at Node A



	/* Reconciliation at Node B */
			/* Create blocks from data_b*/
			/* 1st block length: bch->n - bch->ecc_bits */
			printf("\ndata_b_block0		= ");
			for (i =0; i < bch->n - bch->ecc_bits; i++) {
				data_b_block0[i] = message_b[i];
				printf("%d", data_b_block0[i]);
			}
			/* 2nd block length: nrPkt - (bch->n - bch->ecc_bits) */
			// printf("\ndata_b_block1		= ");
			// for (i =0; i < nrPkt - (bch->n - bch->ecc_bits); i++) {
			// 	data_b_block1[i] = message_b[(bch->n - bch->ecc_bits)+i];
			// 	printf("%d", data_b_block1[i]);
			// }

				/* Padding for the 2nd block to make length: bch->n - bch->ecc_bits */
				// int nrPadding_b =0;
				// nrPadding_b = ((bch->n - bch->ecc_bits)  - (nrPkt - (bch->n - bch->ecc_bits)));
				// printf("\nnrPadding_b: %d", nrPadding_b);
				// memset(data_b_block1 + nrPkt - (bch->n - bch->ecc_bits), 0, nrPadding_b);
				// printf("\ndata_b_block1(padded)	= ");
				// for (i =0; i < 71; i++) {
				// 		printf("%d", data_b_block1[i]);
				// }
		/*Generate r1 by XOR message_b and sketch_a_s (it is sent by Node A) */
		generate_xor_vector(bch->n - bch->ecc_bits, data_b_block0,sketch_a_s , r1_b);

		printf("\nr1			= ");
		for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
			printf("%d", r1_b[i]);
		}
		printf("\n");

		/* Check bit mismatch between mismatch scetch_a, r1_b */

		printf("\nBit mismatch scetch_a_s, r1_b: \n");
		for (i = 0;	i < bch->n - bch->ecc_bits; i++){
			if (r1_b[i] != sketch_a_s[i]){
				printf("%d ", i);
			}
		}
			/* BCH decode r1 */
			unsigned int errloc[t];
			memset(errloc, 0, t);
			nerrors = decodebits_bch(bch, r1_b, ecc_a_block1, errloc);
			printf("\nErrors detected: %d\nat bit-Position: ", nerrors);
			for(i = 0; i < nerrors; i++){
				printf(" %d ", errloc[i]);
			}
			printf("\n");


			/*correction*/
			corrupt_data(errloc, r1_b, nerrors);
			/*after BCH decoding r1 becomes r2*/
			printf("Decoded r1(=r2)		= ");
			for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
				printf("%d", r1_b[i]);
			}
			printf("\n");


				/* bch encode r2,  it becomes r3 */
				memset(ecc_b, 0, bch->ecc_bits);
				encodebits_bch(bch, r1_b, ecc_b);
				printf("\nr3			= ");
				for (j = 0;	j < bch->n - bch->ecc_bits; j++) {
					printf("%d", r1_b[j]);
				}
				printf("\n");


					/*Bitwise Multiplication of r3 and rand_data_a_x*/

					printf("r3 dot x         	= ");
					for (i = 0; i < bch->n - bch->ecc_bits ; i++){
						bit_mul_b[i] = r1_b[i] * rand_data_a_x[i];
						printf("%d", bit_mul_b[i]);

					}
					printf("\n");

						/* XOR sketch_a_s with (r3 dot x) to genrate wclean_b */
						generate_xor_vector(bch->n - bch->ecc_bits, sketch_a_s, bit_mul_b , wclean_b);
						printf("Wclean         		= ");
						for (i = 0; i < bch->n - bch->ecc_bits ; i++){
							printf("%d", wclean_b[i]);
						}
						printf("\n");

						/* Check bit mismatch between mismatch scetch_a, r1_b */

						printf("\nBit mismatch data_b_block0, wclean_b: \n");
						for (i = 0;	i < bch->n - bch->ecc_bits; i++){
							if (wclean_b[i] != data_b_block0[i]){
								printf("%d ", i);
							}
						}



							/* XOR wclean_b with rand_data_a_x  */
							generate_xor_vector(bch->n - bch->ecc_bits, wclean_b, rand_data_a_x , data_b);
							printf("\ndate_b         		= ");
							for (i = 0; i < bch->n - bch->ecc_bits ; i++){
								printf("%d", data_b[i]);
							}
							printf("\n");

							/* Check bit mismatch between mismatch scetch_a, r1_b */

							printf("\nBit mismatch data_b, sketch_a_R: \n");
							for (i = 0;	i < bch->n - bch->ecc_bits; i++){
								if (data_b[i] != sketch_a_R[i]){
									printf("%d ", i);
									k = 0;
								}
							}
					/* Decoding */
					// unsigned int errloc[t];
					// memset(errloc, 0, t);
					// nerrors = decodebits_bch(bch, data_b_block0, ecc_a_block0, errloc);

					// printf("\nNr. Errors in data_b_block0: %d\nat bit-Position: ", nerrors);
					// for(i = 0; i < nerrors; i++){
					// 	printf("%d ", errloc[i]);
					// }
					// printf("\n");

					/*  correcting errors  */
					// correctbits_bch(bch, data_b_block0, errloc, nerrors);

					/* Check if data_b has been successfully decoded (only for debugging)*/
					// printf("\nBit mismatch after decoding message: \n");
					// for (i = 0;	i < bch->n - bch->ecc_bits; i++){
					// 	if (data_a_block0[i] != data_b_block0[i]){
					// 		printf("%d ", i);
					// 		k =0;
					// 	}
					// }
					// if (k != 0 ){
					// 	printf("\nSuccessfully recovered data_a_block0 at node B\n");
					// }

					// unsigned int errloc[t];
					// memset(errloc, 0, t);
					// nerrors = decodebits_bch(bch, data_b_block1, ecc_a_block1, errloc);

					// printf("\nNr. Errors in message_b: %d\nat bit-Position: ", nerrors);
					// for(i = 0; i < nerrors; i++){
					// 	printf("%d ", errloc[i]);
					// }
					// printf("\n");

					/*  correcting errors  */
					// correctbits_bch(bch, data_b_block1, errloc, nerrors);

					/* Check if data_b has been successfully decoded (only for debugging)*/
					// printf("\nBit mismatch after decoding message: \n");
					// for (i = 0;	i < bch->n - bch->ecc_bits; i++){
					// 	if (data_a_block1[i] != data_b_block1[i]){
					// 		printf("%d ", i);
					// 		k =0;
					// 	}
					// }
					// if (k != 0 ){
					// 	printf("\nSuccessfully recovered data_a_block1 at node B\n");
					// }


						/* Reconstructing  */

						// for (i =0; i < bch->n - bch->ecc_bits; i++) {
						// 	data_b[i] = data_b_block0[i];
						// 	printf("%d", data_b[i]);
						// }
						// for (i =0; i < nrPkt - (bch->n - bch->ecc_bits); i++) {
						//  data_b[(bch->n - bch->ecc_bits)+i] = data_b_block1[i];
						// printf("%d", data_b[i]);
						// }
						/* Check if correctly reconstructed (DEBUG)*/
						// printf("\nBit mismatch after decoding message: \n");
						// for (i = 0;	i < nrPkt; i++){
						// 	if (data_a[i] != data_b[i]){
						// 		printf("%d ", i);
						// 		k =0;
						// 	}
						// }
						// if (k != 0 ){
						// 	printf("\nSuccessfully reconstructed \n");
						// }

						// 	/* Recovering message_a from data_b*/
						// 	generate_xor_vector(nrPkt, rand_data_a_x, data_b, recov_message_a);
						// 	printf("\nrecov_message_a			= ");
						// 	for (i = 0;	i < nrPkt; i++) {
						// 		printf("%d", recov_message_a[i]);
						// 	}
						// 	/* Check if correctly recovered (DEBUG)*/
						// 	printf("\nBit mismatch after decoding message: \n");
						// 	for (i = 0;	i < nrPkt; i++){
						// 		if (message_a[i] != recov_message_a[i]){
						// 			printf("%d ", i);
						// 			k =0;
						// 		}
						// 	}
						// 	if (k != 0 ){
						// 		printf("\nSuccessfully recovered message_a at node B \n");
						// 	}


	// function [sketch_a_s, rand_data_a_x, data_a] = secure_sketch_generate(message_a, k)
	// function  data_a = secure_sketch_reproduce(message_b, sketch_a_s, rand_data_a_x, k)


	free(pattern);
	free(data_a);
	free(data_a_block0);
	free(data_a_block1);
	free(data_b);
	free(data_b_block0);
	free(data_b_block1);
	free(rand_data_a_x);
	free_bch(bch);
	time_t end = time(NULL);
	printf("\nTime elapsed: %d sec.", (end - begin));
	printf("\n");

	//return data_a, data_b;
	//return message_a, message_b;
	return 0;
}