/*
 * boot.cxx
 * 	JHT, March 27, 2024
 *
 * Generates bootstrap statistics for a given file. 
 * File is expected to be in .csv format, with the first line containing
 * the number of lines of the file and the remaining lines contain the full
 * dataset
 * 
 */

#include <stdio.h>
#include <string>
#include <vector>

#include "randutils.hpp"
#include "bootstrap.hpp"

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


int main()
{
    randutils::mt19937_rng rng;

    int n_sub = 150; //number of subsamples
    int n_resample = 10; //number of times we regenerate subsamples
    int n_averages = 10; //number of times we test a number of trials
    std::vector<int> n_trials_vec; //vector of bootstrap trials per subsamples

    //Generate vector of boostrap trials per subsample we want to run. This will
    //eventually be an integer, but we need to see how this converges
    printf("numbers of trials to test\n");
//    for (auto i = 2; i <= 2048; i*=2)
    for (auto i = 2; i <= 128; i*=2)
    {
        printf("%d,", i);
        n_trials_vec.push_back(i);
    }
    printf("\n\n");

    printf("n_trials_vec: [");
    for (auto elm : n_trials_vec)
        printf("%d,", elm);
    printf("]\n");

    //Matrix of results : average
    std::vector<std::vector<float>> avg_resample_trials(n_resample);
    for (auto& samp : avg_resample_trials)
    for (auto trial : n_trials_vec)
        samp.push_back(0.);

    //matrix of results : stddev
    std::vector<std::vector<float>> stddev_resample_trials(n_resample);
    for (auto& samp : stddev_resample_trials)
    for (auto trial : n_trials_vec)
        samp.push_back(0.);

    //Read samples from file
    auto samples = get_samples_from_file("all.csv");

    //Loop over random subsamples 
    for (auto resample_itr = 0; resample_itr < n_resample; resample_itr++)
    {
        //Generate a random subsample 
        auto sub_samples_idx = get_unique_subsamples_idx(n_sub, samples, rng);

        std::vector<float> sub_samples(n_sub);
        for (auto i = 0; i < n_sub; i++)
            sub_samples[i] = samples[sub_samples_idx[i]];
      
        printf("Random subsample %d\n", resample_itr);
        printf("Species idx in this random subsample: \n[");
        for (auto idx : sub_samples_idx) 
            printf("%d,", idx);
        printf("]\n");
        
        //loop over the number of bootstrap trials
        for (auto trials_itr = 0; trials_itr < n_trials_vec.size(); trials_itr++)
        { 
            auto ntrials = n_trials_vec[trials_itr];

            printf("Num Trials : %d ,", ntrials);
            std::vector<float> trials_avg(ntrials);

            //Buffer for bootstrapping
            std::vector<float> bs_data(ntrials);
 
            //For this random subsample, run bootstrap a couple times
            for (auto averages_itr = 0; averages_itr < n_averages; averages_itr++)
            {
                bootstrap::simple_samples<float>(ntrials, 
                                                 sub_samples.begin(),
                                                 sub_samples.end(), 
                                                 bs_data, 
                                                 rng);
 
                float avg = 0.;
                for (const auto& elm : bs_data)
                   avg += elm;
                avg /= ntrials;

                trials_avg[averages_itr] = avg;
            }

            //Calculate average of averages
            float avg = 0.;
            for (const auto& elm : trials_avg)
                avg += elm;
            avg /= n_averages;
           
            avg_resample_trials[resample_itr][trials_itr] = avg;

            //Calculate std.dev. of averages
            float stddev = 0.;
            for (const auto& elm : trials_avg)
                stddev += (elm - avg) * (elm - avg);
            stddev = sqrt(stddev / n_averages);
                
            stddev_resample_trials[resample_itr][trials_itr] = stddev;

            printf("After averaging %d trials, avg: %f with stddev : %f\n", ntrials, avg, stddev);

        }

        printf("\n");

    }//end loop over random subsamples
    printf("\n\n");

    //Print results
    printf("Average results matrix:\n");
    for (auto& samp : avg_resample_trials)
    {
        for (auto avg : samp)
            printf("%f,", avg);
        printf("\n");
    }

    printf("\n\n");
    printf("Std.Dev results matrix:\n");
    for (auto& samp : stddev_resample_trials)
    {
        for (auto avg : samp)
            printf("%f,", avg);
        printf("\n");
    }


/*

    printf("with %d n_run, avergage was %f\n", n_run, avg);
*/
    

}
