#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>

// Optimized modular reduction for Mersenne numbers
void mersenne_mod(mpz_t s, mpz_t m, int p) {
    mpz_t q, r;
    mpz_inits(q, r, NULL);
    mpz_tdiv_q_2exp(q, s, p); // q = s >> p
    mpz_tdiv_r_2exp(r, s, p); // r = s & ((1 << p) - 1)
    mpz_add(s, r, q); // s = r + q
    if (mpz_cmp(s, m) >= 0) {
        mpz_sub(s, s, m);
    }
    mpz_clears(q, r, NULL);
}

int lucas_lehmer(int p) {
    if (p == 2) return 1;

    mpz_t s, m;
    mpz_init_set_ui(s, 4);
    mpz_init(m);
    mpz_setbit(m, p); // m = 2^p
    mpz_sub_ui(m, m, 1); // m = 2^p - 1

    for (int i = 1; i < p - 1; i++) {
        mpz_mul(s, s, s); // s = s^2
        mpz_sub_ui(s, s, 2); // s = s - 2
        mersenne_mod(s, m, p); // Optimized modulus for Mersenne numbers
    }

    int result = mpz_cmp_ui(s, 0) == 0;
    mpz_clears(s, m, NULL);
    return result;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <exponent>\n", argv[0]);
        exit(1);
    }

    int p = atoi(argv[1]);
    if (p < 2) {
        fprintf(stderr, "Exponent must be greater than 1.\n");
        exit(1);
    }

    int result = lucas_lehmer(p);
    if (result) {
        printf("2^%d - 1 is a Mersenne prime.\n", p);
    } else {
        printf("2^%d - 1 is not a Mersenne prime.\n", p);
    }

    return 0;
}
