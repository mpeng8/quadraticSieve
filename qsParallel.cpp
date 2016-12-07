#include <vector>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <stdint.h>
#include <gmp.h>
#include <gmpxx.h>
#include <mpfr.h>
#define ENABLE_TIMER
#include "timer.h"
 #include <set> 
#include <iterator>
#include "mathHelper.h"


// Minimal smoothness bound.
const static uint32_t MINIMAL_BOUND = 300;

// Sieving interval length.
const static uint32_t INTERVAL_LENGTH = 100;

// Input size threshold below which we resort to trial division.
const static uint32_t TRIAL_THRESHOLD = 1000000000;

// Maximum input size we can handle (bits).
const static uint32_t MAX_DIGITS = 100;

// GMP random number generator.
extern gmp_randclass rng;

int main(int argc, char* argv[]) {
    mpz_class N  = 63787;
    const float logN = mpz_sizeinbase(N.get_mpz_t(), 2) * std::log(2);
    const float loglogN = std::log(logN);
    const mpz_class sqrtN = sqrt(N);

    const uint32_t B = MINIMAL_BOUND + std::ceil(std::exp(0.55*std::sqrt(logN * loglogN)));
    std::vector<uint32_t> factorBase = generateFactorBase(N, B);


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
    //Find the optimal interval to sieve 
    uint32_t intervalStart = 0;
    uint32_t intervalEnd = 131072;

    uint32_t intervalLength;

    std::vector<uint32_t> masterVec;
    
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    MPI_Status status;

    intervalLength = (intervalEnd - intervalStart+1)/(size);
    uint32_t* intervals;

    START();
    if (rank == 0) {
        //printf("rank 0 -> master node\n");
        intervals = new uint32_t[size];
        printf("Sieving Intervals:");
        for(int i = 0; i < size; i++) {
	    intervals[i] = intervalStart + i*intervalLength;
 	    //printf("%d->%d,", intervals[i], intervals[i]+intervalLength-1);
        }
	//std::cout<<"\n";
    }
    if (rank != 0){
	intervals = new uint32_t[size];
    }
    /* 
    The master process broadcasts the computed initial values 
    to all the other processes.
    */
    MPI_Bcast(intervals, size, MPI_UINT32_T, 0, comm);
    std::cout<<"Processor "<<rank<<" with "<<intervals[rank]<<"->"<<intervals[rank]+intervalLength-1<<"\n";
    std::vector<uint32_t> result;

    std::vector<float> logApprox(intervalLength, 0);  // Approx. 2-logarithms of a^2 - N.
    //Generate log approximations of Q = (a + sqrt(N))^2 - N in the current interval.
    uint32_t x = intervals[rank]+1;

    //compute log(Q(x)) for each x in the interval, where Q(x) = (x+sqrt(N))^2 - N
    for (uint32_t i = 1; i < intervalLength; ++i) {
        const mpz_class Q = (x + sqrtN) * (x + sqrtN) - N;
        logApprox[i] = getLogBig(Q);
        ++x;
    }
    //std::cout<<"rank "<<rank<<": "<<logApprox[1]<<"\n";
    //std::cout<<rank<<" done log approx\n";

    for(uint32_t i = 0; i < factorBase.size(); i++) {
        const uint32_t p = factorBase[i];
        const float log_p = std::log(p) / std::log(2);
        int counter = 0;
        while(startIndex.first[i] < intervals[rank]) {
            startIndex.first[i]+= p;
            counter+=1;
        }
        logApprox[startIndex.first[i]-intervals[rank]] -= log_p*counter;
        while(startIndex.first[i] < intervals[rank] + intervalLength) {
            logApprox[startIndex.first[i] - intervals[rank]] -= log_p;
	    startIndex.first[i]+= p;
        }

	//std::cout<<"ddd";
        if(p == 2) continue;
        counter = 0;
        while(startIndex.second[i] < intervals[rank]) {
            startIndex.second[i]+= p;
            counter += 1;
        }
        logApprox[startIndex.second[i] - intervals[rank]] -= log_p*counter;
        while(startIndex.second[i] < intervalLength + intervals[rank]) {
            logApprox[startIndex.second[i] - intervals[rank]] -= log_p;
 	    startIndex.second[i] +=p;
        }

    }
    //std::cout<<"rank "<<rank<<": "<<logApprox[1]<<"\n";
    //std::cout<<"\n";
    //A value of Q(x) that factors completely over the factor base should theoretically have it's 'log(Q(x))' reduced to zero by this procedure
    const float threshold = std::log(factorBase.back()) / std::log(2);
    //threshold def
    x = intervals[rank];
    for (uint32_t i = 0; i < intervalLength; ++i) {
            //only consider values of Q(x) where our initial log(Q(x)) have been reduced below that threshold
            if (std::fabs(logApprox[i]) < threshold) {
                mpz_class Q = (x + sqrtN) * (x + sqrtN) - N;
                std::vector<uint32_t> factors;
                for (uint32_t j = 0; j < factorBase.size(); ++j) {
                    while (mpz_divisible_ui_p(Q.get_mpz_t(), factorBase[j])) {
                        mpz_divexact_ui(Q.get_mpz_t(), Q.get_mpz_t(), factorBase[j]);
                        factors.push_back(j);
                    }
                }
                if (Q == 1) {
		    if((result.size()!= 0 && result.back() < x)|| result.size() == 0) {
                    result.push_back(x);}
                }
            }
            ++x;
    }

    /* 
    Each worker process sends its sum back to the master process.
    */

    if(rank != 0) {
	MPI_Send(&result[0],result.size(), MPI_UINT32_T, 0, 1, comm);
    } else {
       masterVec = result;
       for(int i = 1; i < size; i++) {
	    MPI_Recv(&result[0], result.size(), MPI_UINT32_T, MPI_ANY_SOURCE, 1, comm, &status);
            masterVec.insert(masterVec.end(), result.begin(), result.end());
       }
    }


    if (rank == 0) {
	std::set<uint32_t> rrrr(masterVec.begin(), masterVec.end());
        STOP("sieve parallel done");
	std::set<uint32_t>::iterator iter;
	for(iter = rrrr.begin(); iter!= rrrr.end(); ++iter){
	    std::cout<<	(*iter)<<",";
	}
        printf("\n");
	result.clear();
    }

  
    MPI_Finalize();

    return 0;

}
