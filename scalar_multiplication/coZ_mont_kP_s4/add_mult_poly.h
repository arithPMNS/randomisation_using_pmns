#ifndef POLY_MULT_ADD
#define POLY_MULT_ADD


void sub_poly(int64_t *rop, int64_t *pa, int64_t *pb);
void add_poly(int64_t *rop, int64_t *pa, int64_t *pb);
void double_add_poly(int64_t *rop, int64_t *pa, int64_t *pb);
void neg_poly(int64_t *rop, int64_t *op);
void scalar_mult_poly(int64_t *rop, int64_t *op, int64_t scalar);
void double_poly_coeffs(int64_t *rop, int64_t *op);
void lshift_poly_coeffs(int64_t *rop, int64_t *op, int nb_pos);

void compute_rand_zero(int64_t *rand_zero, int64_t *rand_pol);

void gen_rand_pol(int64_t *rand_pol);

void mult_mod_poly(int64_t *rop, int64_t *pa, int64_t *pb);

void square_mod_poly(int64_t *rop, int64_t *pa);

void internal_reduction(int64_t *rop, int128 *op);

void exact_coeffs_reduction(int64_t *rop, int64_t *op);

#endif

