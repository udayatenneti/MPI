#include <iostream>
#include <argparse.h>
#include <threads.h>
#include <io.h>
#include <chrono>
#include <cstring>
#include <pthread.h>
#include <math.h>
#include <atomic>

using namespace std;



int main(int argc, char **argv)
{
    // Parse args
    struct options_t opts;
    get_opts(argc, argv, &opts);

    //printf("spin: %d\n", opts.spin);
    bool sequential = false;
    if (opts.n_threads == 0) {
        opts.n_threads = 1;
        sequential = true;
    }

    // Setup threads
    pthread_t *threads = sequential ? NULL : alloc_threads(opts.n_threads);


    // Setup args & read input data
    prefix_sum_args_t *ps_args = alloc_args(opts.n_threads);
    int n_vals;
    int *input_vals, *output_vals;
    bool pad_with_zeroes;

    read_file(&opts, &n_vals, &input_vals, &output_vals, &pad_with_zeroes, &sequential);
    int next_power = n_vals;
    if(pad_with_zeroes){
        //printf("came here wrong 3!\n");
        next_power = next_power_of_two_or_equal(n_vals);
    }
    //printf("next_power %d n_vals %d\n", next_power, n_vals);

    //"op" is the operator you have to use, but you can use "add" to test
    int (*scan_operator)(int, int, int);
    scan_operator = op;
    //scan_operator = add;

    pthread_barrier_t   barrier; // the barrier synchronization object
    pthread_barrier_init(&barrier, NULL, opts.n_threads);  
    int level = log2 (next_power); //padded with zeroes potentially
    std::atomic<int>* go_arr = alloc_barrier_arr(opts.n_threads);
    std::atomic<int> global_counter = 0;
    //TODO: STEP2: Once work at level is complete, replace level-arg with level-1, work-at-level-1 

    fill_args(ps_args, opts.n_threads, n_vals, input_vals, output_vals,
        opts.spin, scan_operator, opts.n_loops, level, &barrier, &global_counter, go_arr, pad_with_zeroes, next_power);

    // Start timer
    auto start = std::chrono::high_resolution_clock::now();

    if (sequential)  {
        //sequential prefix scan
        output_vals[0] = input_vals[0];
        for (int i = 1; i < n_vals; ++i) {
            //y_i = y_{i-1}  <op>  x_i
            output_vals[i] = scan_operator(output_vals[i-1], input_vals[i], ps_args->n_loops);
        }
    }
    else {
        int last_value = input_vals[n_vals-1];
        start_threads(threads, opts.n_threads, ps_args, compute_prefix_sum);

        // Wait for threads to finish
        join_threads(threads, opts.n_threads);
        for (int i = 0; i < n_vals - 1; ++i) {
            //y_i = y_{i-1}  <op>  x_i
            output_vals[i] = input_vals[i+1];
        }
        output_vals[n_vals -1] = output_vals[n_vals - 2] + last_value;
    }

    //End timer and print out elapsed
    auto end = std::chrono::high_resolution_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    //printf("no of threads: %d\n", opts.n_threads);
    std::cout << "time: " << diff.count() << std::endl;

    // Write output data
    write_file(&opts, &(ps_args[0]));

    // Free other buffers
    free(threads);
    free(ps_args);
    pthread_barrier_destroy(&barrier);
}
