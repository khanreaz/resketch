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
#include "2019-07-10_16-03-24_csi_log.c"
// #include "2019-07-10_15-12-45_csi_log.c"

#define MAX_ERRORS 2048
static int bitflip[MAX_ERRORS];

int main()
{	/* Common */
	time_t begin = time(NULL);
	int m = 7; // BCH polynomial degree
	int t = 9; // erroc correction capability
	int len = 0;
	int nrPkt = 128; // to produce 128 bits
	int i, j, k=0, nerrors; //c, nerrors_b,opt_cache_encode = 0, opt_decode = 0;
	double *pe_sum;
	struct bch_control *bch;
	unsigned int *pattern = NULL, generator = 0, tmax, errloc[t];

	/* for Node A */
	uint8_t *message_a, *data_a_block0, *data_a_block1, *rand_data_a_x, *ecc_a_block0, *ecc_a_block1, *rand_data_a_k, *bit_mul_a, *sketch_a_R, *sketch_a_r, *sketch_a_s_block0, *sketch_a_s_block1;
	double *pe_data_a;

	/* for Node B*/
	uint8_t *message_b,*data_b_block0, *data_b_block1,  *ecc_b, *recov_sketch_R_block0, *recov_sketch_R_block1, *bit_mul_b, *r1_b, *r2_b, *r3_b, *wclean_b, *recov_block0, *recov_block1;
	double *pe_data_b;


	srand(time(NULL)); // Seed for the random number generator (needed at Node A)

	/* common */
	assert((m >= 5) && (m <= 15));
	if (len == 0) {
		len = 1 << (m-4);
	}
	assert((t > 0) && (t <= MAX_ERRORS));
	assert((len > 0) && (8*len+m*t <= (1 << m)-1));

	bch = init_bch(m, t, generator);
	assert(bch);

	/* for Node A */
	message_a = malloc(nrPkt); // Alice: for quantized bits
	assert(message_a);
	data_a_block0 = malloc(bch->n - bch->ecc_bits);
	assert(data_a_block0);
	data_a_block1 = malloc(bch->n - bch->ecc_bits);
	assert(data_a_block1);

	rand_data_a_x = malloc(nrPkt);
	assert(rand_data_a_x);
	rand_data_a_k = malloc(nrPkt);
	assert(rand_data_a_k);

	ecc_a_block0 = malloc(bch->ecc_bits);
	assert(ecc_a_block0);
	ecc_a_block1 = malloc(bch->ecc_bits);
	assert(ecc_a_block1);

	bit_mul_a = malloc(bch->n+bch->ecc_bits);
	assert(bit_mul_a);

	sketch_a_R = malloc(bch->n+bch->ecc_bits);
	assert(sketch_a_R);
	sketch_a_r = malloc(bch->n+bch->ecc_bits);
	assert(sketch_a_r);
	sketch_a_s_block0 = malloc(bch->n+bch->ecc_bits);
	assert(sketch_a_s_block0);
	sketch_a_s_block1 = malloc(bch->n+bch->ecc_bits);
	assert(sketch_a_s_block1);

	/* for Node B */
	message_b = malloc(nrPkt); // Bob: for quantized bits
	assert(message_b);
	data_b_block0 = malloc(bch->n - bch->ecc_bits);
	assert(data_b_block0);
	data_b_block1 = malloc(bch->n - bch->ecc_bits);
	assert(data_b_block1);

	recov_sketch_R_block0 = malloc(bch->n - bch->ecc_bits);
	assert(recov_sketch_R_block0);
	recov_sketch_R_block1 = malloc(bch->n - bch->ecc_bits);
	assert(recov_sketch_R_block1);

	recov_block0 = malloc(bch->n - bch->ecc_bits);
	assert(recov_block0);
	recov_block1 = malloc(bch->n - bch->ecc_bits);
	assert(recov_block1);

	ecc_b = malloc(bch->ecc_bits);
	assert(ecc_b);

	bit_mul_b = malloc(bch->n+bch->ecc_bits);
	assert(bit_mul_b);

	r1_b = malloc(bch->n+bch->ecc_bits);
	assert(r1_b);
	r2_b = malloc(bch->n+bch->ecc_bits);
	assert(r2_b);
	r3_b = malloc(bch->n+bch->ecc_bits);
	assert(r3_b);
	wclean_b = malloc(bch->n+bch->ecc_bits);
	assert(wclean_b);

	/* EOVariables*/

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
	/* EOQunatization (1) */

	/* Quantization method (2): "x" point moving window  based */
	printf("**** Quantization (Moving Window) ****");

	int wSize = 3; // Set window size, x: nrPkt > wSize > 1.
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


	printf("\nmessage_a 		= ");
    for (i = 0;	i < nrPkt; i++) {
        printf("%d", message_a[i]);
    }

	printf("\nmessage_b 		= ");
    for (i = 0;	i < nrPkt; i++) {
        printf("%d", message_b[i]);
    }

	/* Check bit mismatch between A,B after quantization (DEBUG)*/

	printf("\nBit mismatch message_a, message_b: ");
	for (i = 0;	i < nrPkt; i++){
		if (message_a[i] != message_b[i]){
			printf("%d ", i);
			k++;

		}
	}
	printf("\nNr. of errors		:%d", k-1);
	printf("\n");


			/* Create blocks from data_a*/
			/* 1st block length: bch->n - bch->ecc_bits */
			printf("\n**** Node A encodes data_a_block0 ****\n");
			printf("\ndata_a_block0		= ");
			for (i =0; i < bch->n - bch->ecc_bits; i++) {
				data_a_block0[i] = message_a[i];
				printf("%d", data_a_block0[i]);
			}

				/* Generate Random  Vector as same length as data block  */
				generate_random_vector(bch->n - bch->ecc_bits , rand_data_a_x);
				printf("\nrand_data_a_x		= ");
				for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
					printf("%d", rand_data_a_x[i]);
				}
				// printf("\n");

					/* XOR data block with the  random_vector  */
					generate_xor_vector(bch->n - bch->ecc_bits, data_a_block0, rand_data_a_x, sketch_a_R);
					printf("\nsketch_a_R		= ");
					for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
						printf("%d", sketch_a_R[i]);
					}
					// printf("\n");




					/* */
						/* Generate additional Random  Vector as same length as data block */
 						generate_random_vector(bch->n - bch->ecc_bits, rand_data_a_k);
 						// memset(ecc_a_block0, 0, bch->ecc_bits);// memset(ecc_a_block0, 0, bch->ecc_bits);
						// encodebits_bch(bch, rand_data_a_k, ecc_a_block0);

						printf("\nrand_data_a_k		= ");
 						for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
 						    printf("%d", rand_data_a_k[i]);
 						}
 						printf("\n");


 						/*  multiplication of x and r */
 						printf("r dot x         	= ");
 						for (i = 0; i < bch->n - bch->ecc_bits ; i++){
 							bit_mul_a[i] = rand_data_a_x[i] * rand_data_a_k[i];
 							printf("%d", bit_mul_a[i]);
 						}
 						printf("\n");
						/* encoding bit_mul_a  */
						memset(ecc_a_block0, 0, bch->ecc_bits);
						encodebits_bch(bch, bit_mul_a, ecc_a_block0);
 						printf("ecc_a_block0 		= ");
 						for (i = 0; i < bch->ecc_bits ; i++){
 							printf("%d", ecc_a_block0[i]);
 						}

 						/* XOR data block with encoded bit_mul_a to produce sketch_a_s_block0  */
 						generate_xor_vector(bch->n - bch->ecc_bits, data_a_block0, bit_mul_a, sketch_a_s_block0);
 						printf("\nsketch_a_s_block0	= ");
 						for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
 						    printf("%d", sketch_a_s_block0[i]);
 						}
 						printf("\n");
						/* EOEncoding at Node A */
						/* Node A sends sketch_a_s_block0, rand_data_a_x and ecc_a_block0 to Node B */


printf("\n/* Node A sends sketch_a_s_block0, rand_data_a_x and ecc_a_block0 to Node B */\n");

	/* Reconciliation at Node B */
			/* Create blocks from message_b*/
			/* 1st block length: bch->n - bch->ecc_bits */
			printf("\n**** Node B reconciles data_b_block0 ****\n");
			printf("\ndata_b_block0		= ");
			for (i =0; i < bch->n - bch->ecc_bits; i++) {
				data_b_block0[i] = message_b[i];
				printf("%d", data_b_block0[i]);
			}


		/* Reconciliation of data_block0*/
		/*Generate r1 by XOR data_b_block0 and sketch_a_s_block0 (sketch_a_s_block0 is sent by Node A) */
		generate_xor_vector(bch->n - bch->ecc_bits, data_b_block0,sketch_a_s_block0 , r1_b);

		// printf("\nr1			= ");
		for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
			r1_b[i];
			// printf("%d", r1_b[i]);
		}
		// printf("\n");

		/* Check bit mismatch between mismatch scetch_a, r1_b */

		// printf("Bit mismatch scetch_a_s, r1_b: ");
		// for (i = 0;	i < bch->n - bch->ecc_bits; i++){
			// if (r1_b[i] != sketch_a_s_block0[i]){
				// printf("%d ", i);
			// }
		// }
			/* BCH decode r1 with the recieved ecc_a_block0 */

			memset(errloc, 0, t);
			nerrors = decodebits_bch(bch, r1_b, ecc_a_block0, errloc);
			printf("\nErrors detected: %d, at bit-Position: ", nerrors);
			for(i = 0; i < nerrors; i++){
				printf(" %d ", errloc[i]);
			}
			printf("\n");


			/*correction*/
			corrupt_data(errloc, r1_b, nerrors);
			/*after BCH decoding r1 becomes r2, variable name DOES NOT change! */
			// printf("Decoded r1(=r2)		= ");
			// for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
			// 	printf("%d", r1_b[i]);
			// }
			// printf("\n");


				/* bch encode r2,  it becomes r3 , variable name DOES NOT change!*/
				memset(ecc_b, 0, bch->ecc_bits);
				encodebits_bch(bch, r1_b, ecc_b);
				// printf("r3			= ");
				// for (j = 0;	j < bch->n - bch->ecc_bits; j++) {
				// 	printf("%d", r1_b[j]);
				// }
				// printf("\n");


					/*Bitwise Multiplication of r3 and rand_data_a_x*/

					// printf("r3 dot x         	= ");
					for (i = 0; i < bch->n - bch->ecc_bits ; i++){
						bit_mul_b[i] = r1_b[i] * rand_data_a_x[i];
						// printf("%d", bit_mul_b[i]);

					}

						/* XOR sketch_a_s_block0 with (r3 dot x) to genrate wclean_b */
						generate_xor_vector(bch->n - bch->ecc_bits, sketch_a_s_block0, bit_mul_b , wclean_b);
						// printf("Wclean         		= ");
						// for (i = 0; i < bch->n - bch->ecc_bits ; i++){
						// 	printf("%d", wclean_b[i]);
						// }
						// printf("\n");

						/* Check bit mismatch between mismatch data_b_block0, wclean_b */

						// printf("Bit mismatch data_a_block0, wclean_b: ");
						// for (i = 0;	i < bch->n - bch->ecc_bits; i++){
						// 	if (wclean_b[i] != data_a_block0[i]){
						// 		printf("%d ", i);
						// 	}
						// }


							/* XOR wclean_b with rand_data_a_x  */
							generate_xor_vector(bch->n - bch->ecc_bits, wclean_b, rand_data_a_x , recov_sketch_R_block0);
							printf("\nrecov_sketch_R_block0	= ");
							for (i = 0; i < bch->n - bch->ecc_bits ; i++){
								printf("%d", recov_sketch_R_block0[i]);
							}
							printf("\n");

							/* Check bit mismatch between mismatch recov_sketch_R_block0, sketch_a_R */

							printf("Bit mismatch recov_sketch_R_block0, sketch_a_R: \n");
							for (i = 0;	i < bch->n - bch->ecc_bits; i++){
								if (recov_sketch_R_block0[i] != sketch_a_R[i]){
									printf("%d ", i);
								}
							}

							/* Recovering data_a_block0 */
							generate_xor_vector(bch->n - bch->ecc_bits, rand_data_a_x, recov_sketch_R_block0, recov_block0);
							printf("\nrecov_block0		= ");
							for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
								printf("%d", recov_block0[i]);
							}
							printf("\ndata_a_block0		= ");
							for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
								printf("%d", data_a_block0[i]);
							}
							/* Check if correctly recovered (DEBUG)*/
							printf("\nBit mismatch after recovering : \n");
							for (i = 0;	i < bch->n - bch->ecc_bits; i++){
								if (data_a_block0[i] != recov_block0[i]){
									printf("%d ", i);
									k =0;
								}
							}
							if (k != 0 ){
								printf("\nSuccessfully recovered data_a_block0 at node B \n");
							}
							// EOReconciliation of data_block0

	/*  **** Node A encodes data_a_block1 **** */

	printf("\n**** Node A encodes data_a_block1 ****\n");

			/* 2nd block length: nrPkt - (bch->n - bch->ecc_bits) */
			// printf("\ndata_a_block1		= ");
			for (i =0; i < nrPkt - (bch->n - bch->ecc_bits); i++) {
				data_a_block1[i] = message_a[(bch->n - bch->ecc_bits)+i];
				// printf("%d", data_a_block1[i]);
			}
				/* Padding for the 2nd block to make length: bch->n - bch->ecc_bits */
				int nrPadding =0;
				nrPadding = ((bch->n - bch->ecc_bits)  - (nrPkt - (bch->n - bch->ecc_bits)));
				// printf("\nnrPadding: %d", nrPadding);
				memset(data_a_block1 + nrPkt - (bch->n - bch->ecc_bits), 0, nrPadding);
				printf("\ndata_a_block1(padded)	= ");
				for (i =0; i < bch->n - bch->ecc_bits; i++) { // TODO: remove the constant!
					data_a_block1[i];
					printf("%d", data_a_block1[i]);
				}
				/* Generate Random  Vector as same length as data block  */
				generate_random_vector(bch->n - bch->ecc_bits , rand_data_a_x);
				printf("\nrand_data_a_x		= ");
				for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
					printf("%d", rand_data_a_x[i]);
				}
				// printf("\n");

					/* XOR "data_a_block1 bits with the  random_vector  */
					generate_xor_vector(bch->n - bch->ecc_bits, data_a_block1, rand_data_a_x, sketch_a_R);
					printf("\nsketch_a_R		= ");
					for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
						printf("%d", sketch_a_R[i]);
					}
					// printf("\n");

						/* Generate additional Random  Vector  */
 						generate_random_vector(bch->n - bch->ecc_bits, rand_data_a_k);

						printf("\nrand_data_a_k		= ");
 						for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
 						    printf("%d", rand_data_a_k[i]);
 						}
 						printf("\n");


 						/* bitwise multiplication of x and r */
 						printf("r dot x         	= ");
 						for (i = 0; i < bch->n - bch->ecc_bits ; i++){
 							bit_mul_a[i] = rand_data_a_x[i] * rand_data_a_k[i];
 							printf("%d", bit_mul_a[i]);
 						}

						/* encoding bit_mul_a */
						memset(ecc_a_block1, 0, bch->ecc_bits);
						encodebits_bch(bch, bit_mul_a, ecc_a_block1);
						printf("\necc_a_block1 		= ");
 						for (i = 0; i < bch->ecc_bits ; i++){
 							printf("%d", ecc_a_block1[i]);
 						}


 						/* XOR data_a_block1 bits with bit_mul_a  */
 						generate_xor_vector(bch->n - bch->ecc_bits, data_a_block1, bit_mul_a, sketch_a_s_block1);

 						printf("\nsketch_a_s_block1 	= ");
 						for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
 						    printf("%d", sketch_a_s_block1[i]);
 						}
						 printf("\n");
						//EOEncoding at Node A for data_block1

						printf("\n/* Node A sends sketch_a_s_block1, rand_data_a_x(it has changed for block1!) and ecc_a_block1 to Node B */\n");

/* At Node B*/
printf("\n**** Node B reconciles data_b_block1 ****\n");

			/* 2nd block length: nrPkt - (bch->n - bch->ecc_bits) */
			// printf("\ndata_b_block1		= ");
			for (i =0; i < nrPkt - (bch->n - bch->ecc_bits); i++) {
				data_b_block1[i] = message_b[(bch->n - bch->ecc_bits)+i];
				// printf("%d", data_b_block1[i]);
			}

				/* Padding for the 2nd block to make length: bch->n - bch->ecc_bits */
				int nrPadding_b =0;
				nrPadding_b = ((bch->n - bch->ecc_bits)  - (nrPkt - (bch->n - bch->ecc_bits)));
				// printf("\nnrPadding_b: %d", nrPadding_b);
				memset(data_b_block1 + nrPkt - (bch->n - bch->ecc_bits), 0, nrPadding_b);
				printf("\ndata_b_block1(padded)	= ");
				for (i =0; i < bch->n - bch->ecc_bits; i++) {
					data_b_block1[i];
					printf("%d", data_b_block1[i]);
				}

				/*Generate r1 by XOR data_b_block1 and sketch_a_s_block1 (sketch_a_s_block1 is sent by Node A) */
				generate_xor_vector(bch->n - bch->ecc_bits, data_b_block1,sketch_a_s_block1 , r1_b);

				// printf("\nr1			= ");
				for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
					r1_b[i];
					// printf("%d", r1_b[i]);
				}
				printf("\n");

		/* Check bit mismatch between mismatch sketch_a_s_block1, r1_b */

		// printf("Bit mismatch scetch_a_s, r1_b: ");
		// for (i = 0;	i < bch->n - bch->ecc_bits; i++){
		// 	if (r1_b[i] != sketch_a_s_block1[i]){
		// 		// printf("%d ", i);
		// 	}
		// }
			/* BCH decode r1 with ecc_a_block1 (sent by A)  */
			// unsigned int errloc[t];
			memset(errloc, 0, t);
			nerrors = decodebits_bch(bch, r1_b, ecc_a_block1, errloc);
			printf("\nErrors detected: %d, at bit-Position: ", nerrors);
			for(i = 0; i < nerrors; i++){
				printf(" %d ", errloc[i]);
			}

			/*correction*/
			corrupt_data(errloc, r1_b, nerrors);
			/*after BCH decoding r1 becomes r2*/
			// printf("Decoded r1(=r2)		= ");
			for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
				r1_b[i];
				// printf("%d", r1_b[i]);
			}
			// printf("\n");


				/* bch encode r2,  it becomes r3 */
				memset(ecc_b, 0, bch->ecc_bits);
				encodebits_bch(bch, r1_b, ecc_b);
				// printf("r3			= ");
				for (j = 0;	j < bch->n - bch->ecc_bits; j++) {
					r1_b[i];
					// printf("%d", r1_b[j]);
				}
				// printf("\n");


					/* Multiplication of r3 and rand_data_a_x*/

					// printf("r3 dot x         	= ");
					for (i = 0; i < bch->n - bch->ecc_bits ; i++){
						bit_mul_b[i] = r1_b[i] * rand_data_a_x[i];
						// printf("%d", bit_mul_b[i]);

					}
					// printf("\n");

						/* XOR sketch_a_s_block1 with (r3 dot x) to genrate wclean_b */
						generate_xor_vector(bch->n - bch->ecc_bits, sketch_a_s_block1, bit_mul_b , wclean_b);
						// printf("\nwclean_b         	= ");
						// for (i = 0; i < bch->n - bch->ecc_bits ; i++){
						// 	printf("%d", wclean_b[i]);
						// }

						/* Check bit mismatch between mismatch data_a_block1, wclean_b */

						// printf("Bit mismatch data_a_block1, wclean_b: ");
						// for (i = 0;	i < bch->n - bch->ecc_bits; i++){
						// 	if (wclean_b[i] != data_b_block1[i]){
						// 		printf("%d ", i);
						// 	}
						// }

							/* XOR wclean_b with rand_data_a_x  */
							generate_xor_vector(bch->n - bch->ecc_bits, wclean_b, rand_data_a_x , recov_sketch_R_block1);
							printf("\nrecov_sketch_R_block1	= ");
							for (i = 0; i < bch->n - bch->ecc_bits ; i++){
								printf("%d", recov_sketch_R_block1[i]);
							}
							printf("\n");

							/* Check bit mismatch between mismatch recov_sketch_R_block1, sketch_a_R */

							printf("Bit mismatch recov_sketch_R_block1, sketch_a_R: \n");
							for (i = 0;	i < bch->n - bch->ecc_bits; i++){
								if (recov_sketch_R_block1[i] != sketch_a_R[i]){
									printf("%d ", i);
									k = 0;
								}
							}
							/* Recovering data_a_block1 */
							generate_xor_vector(bch->n - bch->ecc_bits, rand_data_a_x, recov_sketch_R_block1, recov_block1);
							printf("\nrecov_block1		= ");
							for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
								printf("%d", recov_block1[i]);
							}
							printf("\ndata_a_block1		= ");
							for (i = 0;	i < bch->n - bch->ecc_bits; i++) {
								printf("%d", data_a_block1[i]);
							}
							/* Check if correctly recovered (DEBUG)*/
							printf("\nBit mismatch after recvering data_a_block1: \n");
							for (i = 0;	i < bch->n - bch->ecc_bits; i++){
								if (data_a_block1[i] != recov_block1[i]){
									printf("%d ", i);
									k =0;
								}
							}
							if (k != 0 ){
								printf("\nSuccessfully recovered data_a_block1 at node B \n");
							}

						/* Reconstruction  */

						for (i =0; i < bch->n - bch->ecc_bits; i++) {
							message_b[i] = recov_block0[i];
							// printf("%d", message_b[i]);
						}
						for (i =0; i < nrPkt - (bch->n - bch->ecc_bits); i++) {
							message_b[(bch->n - bch->ecc_bits)+i] = recov_block1[i];
						// printf("%d", message_b[i]);
						}
						/* Check if correctly reconstructed */
						printf("\nBit mismatch after reconstruction: \n");
						for (i = 0;	i < nrPkt; i++){
							if (message_a[i] != message_b[i]){
								printf("%d ", i);
								k =0;
							}
						}
						if (k != 0 ){
							printf("\nSuccessfully reconstructed from data blocks\n");
						}

							printf("\nNow message_b is same as message_a");
							printf("\nmessage_a 		= ");
							for (i = 0;	i < nrPkt; i++) {
								printf("%d", message_a[i]);
							}

							printf("\nmessage_b 		= ");
							for (i = 0;	i < nrPkt; i++) {
								printf("%d", message_b[i]);
							}
							printf("\nNr. of Bits		= %d", i);


	/* Calculate Metric Entropy */
	printf("\n**** Metric entropy of message_a ****\n");
	 double entropy_H;
	// H(X) = [-[(freq_0)log2(freq_0))+((freq_1)log2(freq_1))]] / number of bits ;


		/* Calculate frequency of each alphabet */
		double freq_0;
		double c;

		for (int i = 0; i < nrPkt; i++) {
			if (message_a[i] == 0)
				++c;
		}
		freq_0 = c / nrPkt;
		printf("\nfreq_0: %f", freq_0);


		double freq_1;
		c =0.0;

		for (int i = 0; i < nrPkt; i++) {
			if (message_a[i] == 1)
				++c;
		}
		freq_1 = c / nrPkt;
		printf("\nfreq_1: %f", freq_1);

		entropy_H = (- ((freq_0 * log2(freq_0)) + (freq_1 * log2(freq_1)))) / nrPkt;
		printf("\nEntropy of messsage_a: %f", entropy_H);

		/*  Compare entropy of message_a with PRNG*/
		printf("\n**** Entropy of entropy_PRNG ****\n");
		double entropy_PRNG;

		generate_random_vector(nrPkt , rand_data_a_x);

		/* Calculate frequency of each alphabet */
		// double freq_0;
		// double c;

		for (int i = 0; i < nrPkt; i++) {
			if (rand_data_a_x[i] == 0)
				++c;
		}
		freq_0 = c / nrPkt;
		printf("\nfreq_0: %f", freq_0);


		// double freq_1;
		c =0.0;

		for (int i = 0; i < nrPkt; i++) {
			if (rand_data_a_x[i] == 1)
				++c;
		}
		freq_1 = c / nrPkt;
		printf("\nfreq_1: %f", freq_1);

		entropy_PRNG = (- ((freq_0 * log2(freq_0)) + (freq_1 * log2(freq_1)))) / nrPkt;
		printf("\nMetic entropy of PRNG: %f\n", entropy_PRNG);







	/* Common */
	free(pattern);
	free_bch(bch);

	/* Node A */
	free(data_a_block0);
	free(data_a_block1);
	free(rand_data_a_k);
	free(rand_data_a_x);

	/* Node B */
	free(data_b_block0);
	free(data_b_block1);
	// free(rand_data_a_x); // Node B also has it so needs to be freed!


	time_t end = time(NULL);
	printf("\nTotal time elapsed: %d sec.", (end - begin));
	printf("\n");


	//return message_a;
	//return message_b;
	return 0;
}