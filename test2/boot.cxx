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

#include <stdio.h>
#include <string>
#include <vector>

#include "randutils.hpp"
#include "bootstrap.hpp"

//Number of times we resample the population
//Unsure how many iterations are needed to converge this
#define N_RESAMPLE 100

//Number of trials to generate statistics per population. Testing
//confirmed this converges by 4096 
#define N_TRIAL 4096

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
std::vector<float> get_samples_from_file(const std::string& fname)
{
    FILE* fptr = fopen(fname.c_str(), "r");
    if (nullptr == fptr) exit(1); 

    int nsamples = count_lines(fptr);
    printf("There will be %d samples\n", nsamples);

    if (nsamples <= 0) exit(1);

    std::vector<float> samples(nsamples);

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
std::vector<int> get_unique_subsamples_idx(int n, 
                                           const std::vector<float>& s,
                                           randutils::mt19937_rng& rng)
{

    std::vector<int> idx(s.size()); 
    std::vector<int> smp(n); 
    for (auto i = 0; i < idx.size(); i++) idx[i] = i;

    for (auto i = 0; i < n-1; i++)
    {
        int ii = rng.uniform(0, (int) idx.size()-1);
//        printf("itr %d, we will select position %d, size is %d\n", i, ii, (int) idx.size());
        smp[i] = idx[ii]; 
        idx.erase(idx.begin()+ii);
    }

    int ii = rng.uniform(0, (int) idx.size()-1);
//    printf("Last itr, we will select position %d, size is %d\n", ii, (int) idx.size());
    smp[n-1] = idx[ii];

    return smp;
}

void to_file(const std::string& fname, const std::vector<float>& data)
{
    FILE* fptr = fopen(fname.c_str(), "w");
    if (nullptr == fptr) exit(1);
    for (auto i = 0; i < data.size(); i++)
        fprintf(fptr, "%f\n", data[i]); 
    fclose(fptr);
}

/*
 * Process program input
 */ 
std::vector<int> proc_input(int argc, char* argv[])
{
    int N = argc - 1;

    if (N <= 0) 
    {
        printf("Number of samples <=  0\n");
        exit(1);
    }

    std::vector<int> N_vec(N); 

    for (auto i = 1; i <= N; i++)
    {
        int nn = atoi(argv[i]);
        if (nn <= 0)
        {
            printf("Input %d has sample size <= 0 : %d \n", i, nn);
            exit(1);
        } 

        N_vec[i-1] = nn;
  
    } 
 
    return N_vec;
}


int main(int argc, char* argv[])
{
    //Process input
    auto N_vec = proc_input(argc, argv);
    auto N_samples = N_vec.size();

    //Total number of species
    int N_total = 0;
    for (auto N_sub : N_vec)
        N_total += N_sub;

    //Buffer to contain values for hypothetical sample data
    std::vector<float> buf(N_total);
 
    //Buffer of offsets for random subsamples
    std::vector<size_t> off(N_samples);
    off[0] = 0;
    for (size_t i = 1; i < N_samples; i++)
        off[i] = off[i-1] + N_vec[i-1]; 

    //Buffer that holds indices of random subsamples
    std::vector<std::vector<int>> subsamples_idx; 
    subsamples_idx.resize(N_samples);
    for (size_t i = 0; i < N_samples; i++)
       subsamples_idx[i].resize(N_vec[i]);

    //Initialize random number generator
    randutils::mt19937_rng rng;

    //Buffers for the trials
    std::vector<float> mse(N_TRIAL);
    std::vector<float> mae(N_TRIAL);
    std::vector<float> dev(N_TRIAL);

    //Buffers for tracking across resamples 
    std::vector<float> mse_avg_samples(N_RESAMPLE);
    std::vector<float> mse_dev_samples(N_RESAMPLE);
    std::vector<float> mae_avg_samples(N_RESAMPLE);
    std::vector<float> mae_dev_samples(N_RESAMPLE);
    std::vector<float> dev_avg_samples(N_RESAMPLE);
    std::vector<float> dev_dev_samples(N_RESAMPLE);


    //Read in samples from file
    std::vector<std::vector<float>> samples;
    for (int sub_idx = 0; sub_idx < N_samples; sub_idx++) 
    {
        std::string fname = "data_" + std::to_string(sub_idx) + ".csv";
        samples.push_back(get_samples_from_file(fname));
    }

    //Iterate over the number of resamples for this set
    for (size_t resample_itr = 0; resample_itr < N_RESAMPLE; resample_itr++)
    {
        //Generate the random subsamples for this itr
        for (size_t i = 0; i < N_samples; i++)
            subsamples_idx[i] = get_unique_subsamples_idx(N_vec[i], samples[i], rng);

        //Copy data into buffer
        for (size_t sub = 0; sub < N_samples; sub++)
        for (size_t idx = 0; idx < N_vec[sub]; idx++)
           buf[off[sub] + idx] = samples[sub][subsamples_idx[sub][idx]]; 

        //Run statistics
        bootstrap::simple_samples_avg_stddev<float>(N_TRIAL, 
                                                    buf.begin(),
                                                    buf.end(), 
                                                    mse, 
                                                    mae, 
                                                    dev,
                                                    rng);

        //Generate statistics for this subsample
        float itr_mse = 0.;
        float itr_mae = 0.;
        float itr_dev = 0.;

        for (const auto& elm : mse)
            itr_mse += elm;
        itr_mse /= mse.size();

        for (const auto& elm : mae)
            itr_mae += elm;
        itr_mae /= mae.size();

        for (const auto& elm : dev)
            itr_dev += elm;
        itr_dev /= dev.size();

        mse_avg_samples[resample_itr] = itr_mse;
        mae_avg_samples[resample_itr] = itr_mae;
        dev_avg_samples[resample_itr] = itr_dev;
        
    }

    //From these, generate Bootstrapped MAE, MSE, and DEV of these ratios of subsamples
    float total_mse = 0.;
    float total_mae = 0.;
    float total_dev = 0.;

    for (const auto& elm : mse_avg_samples)
        total_mse += elm;
    total_mse /= mse_avg_samples.size();

    for (const auto& elm : mae_avg_samples)
        total_mae += elm;
    total_mae /= mae_avg_samples.size();

    for (const auto& elm : dev_avg_samples)
        total_dev += elm;
    total_dev /= dev_avg_samples.size();

    printf("Total input size: %zu\n", N_samples);
    printf("Subsample sizes [");
    for (auto& sub : N_vec)
        printf("%d, ", sub);
    printf("]\n");
    printf("MSE  MAE   DEV\n");
    printf("%f %f %f\n", total_mse, total_mae, total_dev);

}
