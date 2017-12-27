#include <gmp.h>
#include <lela/solutions/echelon-form.h>
#include <lela/ring/integers.h>
#include <lela/vector/sparse.h>
#include <lela/ring/gf2.h>
#include <iostream>
#include <random>

using namespace std;
using namespace LELA;

template <class Ring>
void MyAlgorithm (const Ring &R, size_t n)
{
    std::vector<typename Ring::Element> v (n); // A dense vector in R^n

    // Initialise some entries
    R.init (v[0], 1);
    R.init (v[1], -1);
    R.init (v[2], -3);
}

using namespace LELA;

template <class Ring, class Matrix>
void ComputeRowEchelonForm (const Ring &R, Matrix &A)
{
    Context<Ring> ctx (R);
    EchelonForm<Ring> EF (ctx);
    EchelonForm<Ring, Matrix>::echelonize (A, true, EchelonForm<Ring, Matrix>::METHOD_FAUGERE_LACHARTRE); // A now replaced by its reduced row echelon form
}


/**
 *
 * @param d
 * @param n
 * @return
 */
mpz_t* getRandomNumber(int d, int n) {

    random_device rd;                                                       ///random device used for initializing uniform distribution
    unsigned int random_num;                                                ///random number for seed
    unsigned int lower_bound = 1;                                           ///lower bound for seed
    unsigned int upper_bound = USHRT_MAX;                                   ///upper bound for seed
    uniform_int_distribution<unsigned int> unif(lower_bound, upper_bound);  ///declaration of uniform distribution
    int bit = n * (2 - d);                                                  ///number of bits

    mp_bitcnt_t bits; bits = bit;   ///number of bits for gmp
    mpz_t seed;                     ///represents the seed used in each iteration
    mpz_t* random_int;               ///represents the random integer retrieved from gmp
    gmp_randstate_t grt;            ///represents the state of the gmp random machine

    mpz_init_set_ui(seed, random_num);      ///set the number to the seed
    gmp_randinit_default(grt);              ///initialize state
    gmp_randseed(grt, seed);                ///initialize random machine
    mpz_init(*random_int);                   ///initialize the random integer
    mpz_urandomb(*random_int, grt, bits);    ///get a random num of n bits

    return random_int;
}

int main() {

    Integers Z;
    mpz_t MPZ_T;
    GF2 F;

    size_t n = 5;

    MyAlgorithm<mpz_t>(MPZ_T, n);
    DenseMatrix<mpz_t> denseMatrix;
    denseMatrix.resize(10, 10);

    for(int i=0; i<10; ++i) {
        mpz_t* trueRandom = getRandomNumber(1, 80);
        gmp_printf("The shit is now random ... %d", *trueRandom);
        ///@G-Eazy - Random
        denseMatrix[i] = trueRandom;
    }
    ComputeRowEchelonForm<mpz_t, DenseMatrix>(MPZ_T, denseMatrix);

    return 0;
}