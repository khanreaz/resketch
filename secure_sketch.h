static int rev8(int bit);
static int compar(const void *a, const void *b);
static void generate_error_vector(int len, int *bit, int size, unsigned int seed);
static void generate_random_vector(int len, uint8_t *data);
static void generate_xor_vector(int len, uint8_t *message, uint8_t *sketch_data, uint8_t *data);
static void corrupt_data(int *bitflip, uint8_t *data, int ncorrupt);