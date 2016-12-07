#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <stdint.h>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
using namespace std;
float getLogBig(mpz_class num) {

    mpfr_t result, numPP;
    mpfr_init(result);
    mpfr_init(numPP);
    mpfr_set_z(numPP, num.get_mpz_t (), (mpfr_rnd_t) 1);
    mpfr_log2 (result,numPP, (mpfr_rnd_t) 1);
    float pp =  mpfr_get_flt(result, (mpfr_rnd_t)1);
    return pp;

}

uint64_t modularPow(uint64_t b, uint64_t e, uint64_t m) {
    uint64_t result = 1;
    while (e > 0) {
        if (e & 1) {
            result = (result * b) % m;
	}
        e >>= 1;
        b = (b * b) % m;
    }
    return result;
}


int32_t legendreSymbol(uint32_t a, uint32_t p) {
    uint64_t result = modularPow(a, (p - 1) / 2, p);
    return result > 1 ? -1 : result;
}

pair<uint32_t, uint32_t> tonelliShanks(uint32_t n, uint32_t p) {
    if (p == 2) {
        return std::make_pair(n,n);
    } else {
	uint64_t temp = p - 1;
	uint64_t power = 0;
        while (temp%2 != 1) {
	    temp/= 2;
	    power += 1;
	}
	uint64_t z = 2;
	while(legendreSymbol(z,p) != -1) {
	    z+= 1;
	}
        uint64_t c = modularPow(z, temp, p);
	uint64_t t = modularPow(n, temp, p);
        uint64_t R = modularPow(n, (temp+1)/2, p);
	uint64_t M = power;

	while(t %p != 1) {
	    uint64_t i = 1;
	    while (modularPow(t, pow(2, i), p) != 1) {
                i+= 1;
	    }
	    uint64_t b = modularPow(c, pow(2, M-i-1), p);
	    R = (R*b)%p;
	    t = (t*b*b)%p;
	    c = (b * b)%p;
	    M = i;
	}
	return make_pair(R, p-R);
    }
}

vector<uint32_t> generateFactorBase(const mpz_class& N, uint32_t B) {
    vector<uint32_t> result;
    vector<bool> sieve(B + 1, false);
    for(uint32_t i = 2; i <= B; i++) {
	if(sieve[i]) {
            //not a prime or visited
	    continue;
	}
        if(mpz_legendre(N.get_mpz_t(), mpz_class(i).get_mpz_t()) == 1) {
	    result.push_back(i);
	}
        for(uint32_t temp = i; temp <= B; temp+=i) {
	    sieve[temp] = true;
        }
    }
    return result;
}
