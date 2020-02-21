#ifndef __SECURE_SKETCH_H_
#define __SECURE_SKETCH_H_

void generate_random_vector(int len, uint8_t *data);
void generate_xor_vector(int len, uint8_t *message, uint8_t *sketch_data, uint8_t *data);
void corrupt_data(unsigned int *bitflip, uint8_t *data, unsigned int ncorrupt);
void generate_mul_vector(int len, uint8_t *data1, uint8_t *data2, uint8_t *result);

#endif
