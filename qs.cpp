#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdint.h>
#include <gmp.h>
#include <gmpxx.h>
#include <mpfr.h>

#define ENABLE_TIMER
#include "timer.h"
#include "mathHelper.h"
#include "matrix.h"


const static uint32_t MINIMAL_BOUND = 300;
const static uint32_t INTERVAL_LENGTH = 65536;


//  METHOD: Main method to perform quadratic sieve
//
//  Parameters:
//    Input, const mpz_class& N
//    Output, mpz_class num, a nontrivial factor
mpz_class quadraticSieve(const mpz_class& N) {

    const float logN = mpz_sizeinbase(N.get_mpz_t(), 2) * std::log(2);
    const float loglogN = std::log(logN);
    const mpz_class sqrtN = sqrt(N);
    const uint32_t B = MINIMAL_BOUND + std::ceil(std::exp(0.55*std::sqrt(logN * loglogN)));
    //Data Collection Stage
    std::cout <<"Stage One: Data Collection\n";
    //step 1: generate factor base
    std::cout <<"Step 1.1: generating factor base\n";
    START();
    std::vector<uint32_t> factorBase = generateFactorBase(N, B);
    STOP("Generated factor base");
    std::cout <<"Done\n\n";

    std::cout <<"Step 1.2: restore B-smooth numbers\n";
    START();
    std::pair<std::vector<uint32_t>, std::vector<uint32_t> > startIndex(
        std::vector<uint32_t>(factorBase.size()), // Vector of first start index.
        std::vector<uint32_t>(factorBase.size())  // Vector of second start index.
    );
    //Initialization of startIndex vector
    //store x such that (x + sqrt(N))^2 = N (mod p) -> x^2 + 2x*sqrt(N) = 0 (mod p)
    for (uint32_t i = 0; i < factorBase.size(); ++i) {
        uint32_t p = factorBase[i];
        uint32_t N_mod_p = mpz_class(N % p).get_ui();
	//solve for (temp)^2 = N mod(p)
        std::pair<uint32_t, uint32_t> temp = tonelliShanks(N_mod_p, p);
	//store the value of x = temp - sqrt(N)(mod p), keep them positive
        startIndex.first[i] = mpz_class((((temp.first - sqrtN) % p) + p) % p).get_ui();
        startIndex.second[i] = mpz_class((((temp.second - sqrtN) % p) + p) % p).get_ui();
    }
    STOP("Calculated indices");

    //Sieve Stage
    std::cout <<"Stage Two: Sieve\n";

    //for this implementation, we will take a subinterval, 
    uint32_t intervalStart = 0;
    uint32_t intervalEnd = INTERVAL_LENGTH;

    std::vector<uint32_t> smooth;                      // B-smooth numbers.
    std::vector<std::vector<uint32_t> > smoothFactors; // Factorization of each B-smooth number.
    std::vector<float> logApprox(INTERVAL_LENGTH, 0);

    START();
    while (smooth.size() < factorBase.size() + 20) {
        uint32_t x = intervalStart + 1;
        std::cout<<intervalStart<<","<<intervalEnd<<"\n";
	//compute log(Q(x)) for each x in the interval, where Q(x) = (x+sqrt(N))^2 - N
        for (uint32_t i = 1; i < INTERVAL_LENGTH; ++i) {
            const mpz_class Q = (x + sqrtN) * (x + sqrtN) - N;
            logApprox[i] = getLogBig(Q);
	    ++x;
        }
        //std::cout<<"\n";
        //std::cout<<nextLogQfunc<<",\n";
        std::cout<<"\nlog approx done\n";
        //Now loop over all primes in the factor base and for each x in the interval where p|Q(x) we subtract log(p) from log(Q(x))
        for (uint32_t i = 0; i < factorBase.size(); ++i) {
            const uint32_t p = factorBase[i];
            const float log_p = std::log(p) / std::log(2);

            // Sieve first sequence.
            while (startIndex.first[i] < intervalEnd) {
                logApprox[startIndex.first[i] - intervalStart] -= log_p;
                startIndex.first[i] += p;
            }

            if (p == 2) continue; // x^2 = N (mod 2) only has one root.

            // Sieve second sequence.
            while (startIndex.second[i] < intervalEnd) {
                logApprox[startIndex.second[i] - intervalStart] -= log_p;
                startIndex.second[i] += p;
            }
        }
        //A value of Q(x) that factors completely over the factor base should theoretically have it's 'log(Q(x))' reduced to zero by this procedure
        const float threshold = std::log(factorBase.back()) / std::log(2);
        //threshold def
        x = intervalStart;
        for (uint32_t i = 0; i < INTERVAL_LENGTH; ++i) {
            //only consider values of Q(x) where our initial log(Q(x)) have been reduced below that threshold
            if (std::fabs(logApprox[i]) < threshold) {
                mpz_class Q = (x + sqrtN) * (x + sqrtN) - N;
                std::vector<uint32_t> factors;

                // For each factor p in the factor base.
                for (uint32_t j = 0; j < factorBase.size(); ++j) {
                    // Repeatedly divide Q by p until it's not possible anymore.
                    //const uint32_t p = factorBase[j];
                    while (mpz_divisible_ui_p(Q.get_mpz_t(), factorBase[j])) {
                        mpz_divexact_ui(Q.get_mpz_t(), Q.get_mpz_t(), factorBase[j]);
                        factors.push_back(j); // The j:th factor base number was a factor.
                    }
                }
                if (Q == 1) {
                    // Q really was B-smooth, so save its factors and the corresponding x.
                    smoothFactors.push_back(factors);
                    smooth.push_back(x);
                }
                if (smooth.size() >= factorBase.size() + 20)
                    break; // We have enough smooth numbers, so stop factoring.
            }
            ++x;
        }

        // Move on to next interval.
        intervalStart += INTERVAL_LENGTH;
        intervalEnd += INTERVAL_LENGTH;
    }
    std::cout<<intervalEnd<<"\n";
    STOP("sieveing end");
    //for(uint32_t i = 0; i < smooth.size(); i++) {
      //  std::cout<<smooth[i]<<",";
    //}
    std::cout<<"size:" <<smooth.size()<<"\n";

    //Gaussian Stage
    std::cout <<"Stage Three: Gaussian Elminiation\n";
    //Construct a binary matrix M with M_ij = the parity of the i:th prime factor from the factor base in the factorization of the j:th B-smooth number.

    Matrix M(factorBase.size(), smoothFactors.size() + 1);
    for (uint32_t i = 0; i < smoothFactors.size(); ++i) {
        for (uint32_t j = 0; j < smoothFactors[i].size(); ++j) {
            M(smoothFactors[i][j], i).flip();
        }
    }

    //Reduce the matrix to row echelon form and solve it repeatedly until a factor is found.

    M.reduce();
    mpz_class a = 1;
    mpz_class b = 1;
    std::vector<uint32_t> x = M.solve();

    //b = prod(smooth[i]+sqrt(N))
    std::vector<uint32_t> decomp(factorBase.size(), 0);
    for (uint32_t i = 0; i < smoothFactors.size(); ++i) {
        if (x[i] == 1) {
            for(uint32_t p = 0; p < smoothFactors[i].size(); ++p) {
                    ++decomp[smoothFactors[i][p]];
	    }
            b *= (smooth[i] + sqrtN);
        }
    }
    //a = sqrt(product(factorBase[p])).
    for(uint32_t p = 0; p < factorBase.size(); ++p) {
        for(uint32_t i = 0; i < (decomp[p] / 2); ++i) {
            a *= factorBase[p];
	}
    }

    while (a % N == b % N || a % N == (- b) % N + N) {
        x = M.solve();
        a = 1;
        b = 1;
        std::vector<uint32_t> decomp(factorBase.size(), 0);
        for (uint32_t i = 0; i < smoothFactors.size(); ++i) {
            if (x[i] == 1) {
                for(uint32_t p = 0; p < smoothFactors[i].size(); ++p) {
                    ++decomp[smoothFactors[i][p]];
		}
                b *= (smooth[i] + sqrtN);
            }
        }
        for(uint32_t p = 0; p < factorBase.size(); ++p) {
            for(uint32_t i = 0; i < (decomp[p] / 2); ++i) {
                a *= factorBase[p];
	    }
        }
    }
    
    std::cout <<"Found factor, returning...\n";
    mpz_class factor;
    mpz_gcd(factor.get_mpz_t(), mpz_class(b - a).get_mpz_t(), N.get_mpz_t());

    return factor;
    
}
