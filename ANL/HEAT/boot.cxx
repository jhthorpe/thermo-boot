/*
 * boot.cxx
 * 	JHT, March 27, 2024
 *
 * Given a set of data contained in N files called data_x.csv, 
 * where x ranges from 0 -> N-1, generate bootstrapped statistics for with 
 * the whole dataset sampled in the ratios provided by the input. For example, 
 * calling this exetuable with :
 *
 * boot.exe 10 5 7
 * 
 * Would contstruct hypothetical samples of 10 points from data_0.csv, 
 * 5 points from data_1.csv, and 7 points from data_2.csv.
 * 
 */

#include <omp.h>

#include <stdio.h>
#include <string>
#include <vector>
#include <cstdlib>

#include "randutils.hpp"
#include "bootstrap/c++/bootstrap.hpp"

/*
 * Types
 */
using data_t        = float;
using idx_t         = bootstrap::types::idx_t;
using data_buffer_t = bootstrap::types::data_buffer_t<data_t>; 
using idx_buffer_t  = bootstrap::types::idx_buffer_t; 
using data_strata_t = std::vector<data_buffer_t>; 
using idx_strata_t  = std::vector<idx_buffer_t>; 
using rng_t         = bootstrap::types::rng_t; 

//Number of times we resample the population
//Unsure how many iterations are needed to converge this
#define N_RESAMPLE 1

//Number of trials to generate statistics per population. Testing
//confirmed this converges by 4096 
#define N_TRIAL 10000

//#define DEBUG true

/*
 * breakpoint print statment
 */
void breakpoint(const std::string& str)
{
    #if defined DEBUG
    printf("BREAKPOINT:%s\n", str.c_str());
    bool bad = false;
    char line[256];
   
    if (nullptr == fgets(line, sizeof(line), stdin)) exit(EXIT_FAILURE); 
    #endif
}

/*
 * Count number of lines (\n characters) in a file
 */
#define BUF_SIZE 2048
int count_lines(FILE* file)
{
    char buf[BUF_SIZE];
    int counter = 0;
    for(;;)
    {
        size_t res = fread(buf, 1, BUF_SIZE, file);
        if (ferror(file))
            return -1;

        int i;
        for(i = 0; i < res; i++)
            if (buf[i] == '\n')
                counter++;

        if (feof(file))
            break;
    }

    return counter;
}

/*
 * Reads a vector of samples from a given file. This assumes that the first
 * entry of the file contains the number of samples
 */
auto get_samples_from_file(const std::string& fname)
{
    FILE* fptr = fopen(fname.c_str(), "r");
    if (nullptr == fptr) exit(1); 

    int nsamples = count_lines(fptr);
    printf("There will be %d samples\n", nsamples);

    if (nsamples <= 0) exit(1);

    data_buffer_t samples(nsamples);

    rewind(fptr); 

    for (auto i = 0; i < samples.size(); i++)
    {
        fscanf(fptr, "%e", &samples[i]);
    }

    printf("Samples read in were...\n");
    for (auto i = 0; i < samples.size(); i++)
    {
        printf("%d %e\n", i, samples[i]);
    }
    printf("\n");

    fclose(fptr);

    return samples;
}

/*
 * returns a vector of the indices of a unique (no repetitions) 
 * subsamples from a set of samples 
 */
auto get_unique_subsamples_idx(idx_t n, 
                               const data_buffer_t& s,
                               rng_t& rng)
{
    idx_buffer_t idx(s.size()); 
    idx_buffer_t smp(n); 
    for (auto i = 0; i < idx.size(); i++) idx[i] = i;

    for (auto i = 0; i < n-1; i++)
    {
        int ii = rng.uniform(0, (int) idx.size()-1);
        smp[i] = idx[ii]; 
        idx.erase(idx.begin()+ii);
    }

    int ii = rng.uniform(0, (int) idx.size()-1);
    smp[n-1] = idx[ii];

    return smp;
}

/*
 * fills a vector of indices with a unique (no repetitions) 
 * subsamples from a set of samples 
 * Parameters:
 * 	n 	: number of subsamples to take
 *      s	: size of original sample space
 * 	indices : idx_buffer to fill
 *      rng	: random number generator
 */
void fill_unique_subsamples_idx(const idx_t         n, 
                                const idx_t         s, 
                                      idx_buffer_t& indices,
                                      rng_t&        rng)
{
    idx_buffer_t idx(s); 
    for (auto i = 0; i < s; i++) idx[i] = i;

    idx_t end = s - 1;
    for (auto i = 0; i < n-1; i++)
    {
        const auto ii = rng.uniform((idx_t) 0, (idx_t) end); 
        indices[i] = idx[ii]; 
        std::swap(idx[ii], idx[end]);        
        end--;
    }

    const auto ii = rng.uniform((idx_t) 0, (idx_t) end);
    indices[n-1] = idx[ii];

}

/*
 * Given set of samples, copy random subsamples into a buffer
 * n_strata		: number of strata
 * count		: counts of samples to draw from each strata
 * samples		: stratified samples  
 * data			: data buffer to store results in
 * rng			: random number generator
 */
void unique_subsamples(const idx_t          n_strata, 
                       const idx_buffer_t&  counts, 
                       const data_strata_t& samples,
                             data_buffer_t& data,
                             rng_t&         rng) 
{
    //Determine indices of elements to copy per strata
    idx_strata_t indices(n_strata); 
    for (int i = 0; i < n_strata; i++)
    {
        indices[i].resize(counts[i]);
        fill_unique_subsamples_idx(counts[i], samples[i].size(), indices[i], rng);
    }

    //Copy elements according to indices selected 
    idx_t idx = 0;
    for (int i = 0; i < n_strata; i++)
    for (int j = 0; j < counts[i]; j++)
    {
        data[idx] = samples[i][ indices[i][j] ]; 
        idx++;
    }
}

/*
 * Put data buffer in formatted file
 */
void to_file(const std::string& fname, data_buffer_t& data)
{
    FILE* fptr = fopen(fname.c_str(), "w");
    if (nullptr == fptr) exit(1);
    for (auto i = 0; i < data.size(); i++)
        fprintf(fptr, "%f\n", (float) data[i]); 
    fclose(fptr);
}

/*
 * Process program input
 */ 
auto proc_input(int argc, char* argv[])
{
    int N = argc - 1;

    if (N <= 0) 
    {
        printf("Number of samples <=  0\n");
        exit(1);
    }

    idx_buffer_t N_vec(N); 

    for (auto i = 1; i <= N; i++)
    {
        int nn = atoi(argv[i]);
        if (nn <= 0)
        {
            printf("Input %d has sample size <= 0 : %d \n", i, nn);
            exit(1);
        } 

        N_vec[i-1] = (idx_t) nn;
  
    } 
 
    return N_vec;
}

auto trad_conf = [](const idx_t n, const data_buffer_t& vec)
{
    return bootstrap::mean<data_t>(n, vec); 
//    return 2.0 * bootstrap::stddev_p<data_t>(n, vec); 
};


int main(int argc, char* argv[])
{
    //RNG
    rng_t rng;

    //Process input
    auto N_vec = proc_input(argc, argv);
    auto N_strata = N_vec.size();

    //Total number of species across all strata
    int N_total = 0;
    for (auto N_sub : N_vec)
        N_total += N_sub;
    
    printf("N_total %d\n", N_total);

    //Load stratified species from file
    data_strata_t samples(N_strata);
    for (int i = 0; i < N_strata; i++)
        samples[i] = get_samples_from_file("data_" + std::to_string(i) + ".csv"); 

   
    /*
     * Buffers for bootstrap across molecular samples
     */
    data_buffer_t molecule_trials_original(N_RESAMPLE);
    data_buffer_t molecule_trials_nrm_lo(N_RESAMPLE);
    data_buffer_t molecule_trials_nrm_hi(N_RESAMPLE);
    data_buffer_t molecule_trials_pct_lo(N_RESAMPLE);
    data_buffer_t molecule_trials_pct_hi(N_RESAMPLE);
    data_buffer_t molecule_trials_bca_lo(N_RESAMPLE);
    data_buffer_t molecule_trials_bca_hi(N_RESAMPLE);
    data_buffer_t molecule_set(N_total);

    /*
     * BS class is sufficient for the inner loops 
     */
    bootstrap::bs_sim<data_t> bs;

    // Loop over resamples
    for (idx_t resample = 0; resample < N_RESAMPLE; resample++)
    {
        for (idx_t i = 0; i < N_total; i++) molecule_set[i] = (data_t) i;

        //Bootstrap the 95% confidence interval
        bs.sample(N_TRIAL, molecule_set.begin(), molecule_set.begin() + N_total, trad_conf, rng);
        bs.BCa_analysis(trad_conf, rng);

        //Load a random subset of the molecules    
        unique_subsamples(N_strata, N_vec, samples, molecule_set, rng); 

        //Bootstrap the 95% confidence interval
        bs.sample(N_TRIAL, molecule_set.begin(), molecule_set.begin() + N_total, trad_conf, rng);
        bs.BCa_analysis(trad_conf, rng);
        
        molecule_trials_original[resample] = bs.original_average();
        molecule_trials_nrm_lo[resample]   = bs.nrm_bias_corrected_lo_ci(0.95);
        molecule_trials_nrm_hi[resample]   = bs.nrm_bias_corrected_hi_ci(0.95);
        molecule_trials_pct_lo[resample]   = bs.percent_lo_ci(0.95);
        molecule_trials_pct_hi[resample]   = bs.percent_hi_ci(0.95);
        molecule_trials_bca_lo[resample]   = bs.bca_lo_ci(0.95);
        molecule_trials_bca_hi[resample]   = bs.bca_hi_ci(0.95);

    }

    std::sort(molecule_trials_original.begin(), molecule_trials_original.end());
    std::sort(molecule_trials_nrm_lo.begin(), molecule_trials_nrm_lo.end());
    std::sort(molecule_trials_nrm_hi.begin(), molecule_trials_nrm_hi.end());
    std::sort(molecule_trials_pct_lo.begin(), molecule_trials_pct_lo.end());
    std::sort(molecule_trials_pct_hi.begin(), molecule_trials_pct_hi.end());
    std::sort(molecule_trials_bca_lo.begin(), molecule_trials_bca_lo.end());
    std::sort(molecule_trials_bca_hi.begin(), molecule_trials_bca_hi.end());

    auto original_average = bootstrap::mean<data_t>(N_RESAMPLE, molecule_trials_original);
    auto nrm_lo_average = bootstrap::mean<data_t>(N_RESAMPLE, molecule_trials_nrm_lo);
    auto nrm_hi_average = bootstrap::mean<data_t>(N_RESAMPLE, molecule_trials_nrm_hi);
    auto pct_lo_average = bootstrap::mean<data_t>(N_RESAMPLE, molecule_trials_pct_lo);
    auto pct_hi_average = bootstrap::mean<data_t>(N_RESAMPLE, molecule_trials_pct_hi);
    auto bca_lo_average = bootstrap::mean<data_t>(N_RESAMPLE, molecule_trials_bca_lo);
    auto bca_hi_average = bootstrap::mean<data_t>(N_RESAMPLE, molecule_trials_bca_hi);

    printf("Accross all samples average estimates of confidence intervals are:\n");
    printf("ORG[avg=%f               ]\n", original_average);
    printf("NRM[        lo=%f,  hi=%f]\n", nrm_lo_average, nrm_hi_average);
    printf("PCT[        lo=%f,  hi=%f]\n", pct_lo_average, pct_hi_average);
    printf("BCA[        lo=%f,  hi=%f]\n", bca_lo_average, bca_hi_average);

    printf("");
    printf("\n%d, %f, %f, %f, %f, %f, %f, %f\n", (int) N_total, 
                                             original_average,
                                             nrm_lo_average, nrm_hi_average,
                                             pct_lo_average, pct_hi_average,
                                             bca_lo_average, bca_hi_average);

    return EXIT_SUCCESS;    
}
    

