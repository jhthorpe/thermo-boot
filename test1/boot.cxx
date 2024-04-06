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

#define N_SUB 15

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


int main()
{
    randutils::mt19937_rng rng;

    int n_sub = N_SUB; //number of subsamples
    int n_resample = 15; //number of times we regenerate subsamples
    int n_averages = 1; //number of times we test a number of trials
    std::vector<int> n_trials_vec; //vector of bootstrap trials per subsamples

    //Generate vector of boostrap trials per subsample we want to run. This will
    //eventually be an integer, but we need to see how this converges
    printf("numbers of trials to test\n");
    for (auto i = 2; i <= 4096; i*=2)
    {
        printf("%d,", i);
        n_trials_vec.push_back(i);
    }
    printf("\n\n");

    printf("Number of subsamples to test : %d \n", n_sub);

    printf("n_trials_vec: [");
    for (auto elm : n_trials_vec)
        printf("%d,", elm);
    printf("]\n");

    //Matrix of results : average of MSE
    std::vector<std::vector<float>> avg_mse_resample_trials(n_resample);
    for (auto& samp : avg_mse_resample_trials)
    for (auto trial : n_trials_vec)
        samp.push_back(0.);

    //Matrix of results : average of MAE
    std::vector<std::vector<float>> avg_mae_resample_trials(n_resample);
    for (auto& samp : avg_mae_resample_trials)
    for (auto trial : n_trials_vec)
        samp.push_back(0.);

    //matrix of results : average of stddev
    std::vector<std::vector<float>> avg_stddev_resample_trials(n_resample);
    for (auto& samp : avg_stddev_resample_trials)
    for (auto trial : n_trials_vec)
        samp.push_back(0.);

    //non-bootstrapped results
    std::vector<float> nobs_mse_resample(n_resample);
    std::vector<float> nobs_mae_resample(n_resample);
    std::vector<float> nobs_dev_resample(n_resample);

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

        //NON-BS MSE
        float no_mse = 0.;
        for (const auto& elm : sub_samples)
            no_mse += elm;
        no_mse /= n_sub;
        nobs_mse_resample[resample_itr] = no_mse;

        //NON-BS MAE
        float no_mae = 0.;
        for (const auto& elm : sub_samples)
            no_mae += fabs(elm);
        no_mae /= n_sub;
        nobs_mae_resample[resample_itr] = no_mae;

        //NON-BS DEV
        float no_dev = 0.;
        for (const auto& elm : sub_samples)
            no_dev += (elm - no_mse) * (elm - no_mse);
        no_dev /= n_sub;
        nobs_dev_resample[resample_itr] = sqrt(no_dev);

        
        //loop over the number of bootstrap trials
        for (auto trials_itr = 0; trials_itr < n_trials_vec.size(); trials_itr++)
        { 
            auto ntrials = n_trials_vec[trials_itr];

            printf("Num Trials : %d ,", ntrials);
            std::vector<float> trials_mse(n_averages);
            std::vector<float> trials_mae(n_averages);
            std::vector<float> trials_stddev(n_averages);

            //Buffer for bootstrapping
            std::vector<float> bs_mse(ntrials);
            std::vector<float> bs_mae(ntrials);
            std::vector<float> bs_dev(ntrials);
 
            //For this random subsample, run bootstrap a couple times
            for (auto averages_itr = 0; averages_itr < n_averages; averages_itr++)
            {
                #if 1
                bootstrap::simple_samples_avg_stddev<float>(ntrials, 
                                                            sub_samples.begin(),
                                                            sub_samples.end(), 
                                                            bs_mse, 
                                                            bs_mae, 
                                                            bs_dev,
                                                            rng);
                #else
                bootstrap::simple_samples<float>(ntrials, 
                                                 sub_samples.begin(),
                                                 sub_samples.end(), 
                                                 bs_avg, 
                                                 rng);

                #endif 

               
   
                //get average of the mse 
                float mse = 0.;
                for (const auto& elm : bs_mse)
                   mse += elm;
                mse /= ntrials;
                trials_mse[averages_itr] = mse;

                //get average of mae 
                float mae = 0.;
                for (const auto& elm : bs_mae)
                   mae += elm;
                mae /= ntrials;
                trials_mae[averages_itr] = mae;
                
                //get average of dev
                float dev = 0.;
                for (const auto& elm : bs_dev)
                    dev += elm;
                dev /= ntrials;
                trials_stddev[averages_itr] = dev;

                 //Save the best trials
                 if (trials_itr == n_trials_vec.size() - 1)
                 {
                     to_file("mse" + std::to_string(resample_itr) + ".txt", bs_mse);
                     to_file("mae" + std::to_string(resample_itr) + ".txt", bs_mae);
                     to_file("dev" + std::to_string(resample_itr) + ".txt", bs_dev);
                 }
                
            }

            //Calculate average of MSE
            float mse = 0.;
            for (const auto& elm : trials_mse)
                mse += elm;
            mse /= n_averages;
            avg_mse_resample_trials[resample_itr][trials_itr] = mse;

            //Calculate average of MAE
            float mae = 0.;
            for (const auto& elm : trials_mae)
                mae += elm;
            mae /= n_averages;
            avg_mae_resample_trials[resample_itr][trials_itr] = mae;

            //Calculate average of stddevs
            float dev = 0.;
            for (const auto& elm : trials_stddev)
                dev += elm;
            dev /= n_averages;
            avg_stddev_resample_trials[resample_itr][trials_itr] = dev;

        }

        printf("\n");

    }//end loop over random subsamples
    printf("\n\n");

    //Print MSE results
    printf("Non-bootstrapped MSE\n");
    for (auto& mse : nobs_mse_resample)
        printf("%f\n", mse);
    printf("\n");
    printf("Average of MSE results matrix:\n");
    for (auto& samp : avg_mse_resample_trials)
    {
        for (auto avg : samp)
            printf("%f,", avg);
        printf("\n");
    }

    //Print MAE results
    printf("Non-bootstrapped MAE\n");
    for (auto& mae : nobs_mae_resample)
        printf("%f\n", mae);
    printf("\n");
    printf("Average of MAE results matrix:\n");
    for (auto& samp : avg_mae_resample_trials)
    {
        for (auto avg : samp)
            printf("%f,", avg);
        printf("\n");
    }

    //STDDEV results
    printf("\nNon-bootstrapped stdev\n");
    for (auto& dev : nobs_dev_resample)
        printf("%f\n", dev);
    printf("\n");
    printf("Avg. of Std.Dev results matrix:\n");
    for (auto& samp : avg_stddev_resample_trials)
    {
        for (auto avg : samp)
            printf("%f,", avg);
        printf("\n");
    }


}
