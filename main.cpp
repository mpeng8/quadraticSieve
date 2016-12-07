#include <iostream>
#include <vector>
#include <stack>
#include <ctime>
#include <math.h>
#include <limits>
#include <stdlib.h>
#include <stdint.h>
#include <gmpxx.h>
#include "qs.h"

const static uint32_t MINIMAL_BOUND = 300;
const static uint32_t INTERVAL_LENGTH = 65536;
const static uint32_t TRIAL_THRESHOLD = 1000000000;
// Maximum input size we can handle (bits).
const static uint32_t MAX_DIGITS = 100;
extern gmp_randclass rng;

int main(int arg, char *argv[]) {
    // Find some primes for trial division.
    std::vector<uint32_t> primes;
    uint32_t max = ceil(sqrt(TRIAL_THRESHOLD)) + 1;
    std::vector<bool> sieve(max, false);
    for (uint32_t p = 2; p < max; ++p) {
        if (sieve[p])
            continue;
        primes.push_back(p);
        for (uint32_t i = p; i < max; i += p)
            sieve[i] = true;
    }
    //std::cout<<std::numeric_limits<unsigned long long>::max()<<"\n";
    mpz_class N = 612209628037453;

    if (mpz_probab_prime_p(N.get_mpz_t(), 10)) {
        // N is prime.
        std::cout << N << std::endl << std::endl;
    } else{
	std::cout<<"Not a prime!\n";
        mpz_class result = quadraticSieve(N);
        std::cout <<result<<"\n";
    }
    return 0;
}
