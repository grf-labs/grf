#ifndef GRF_PROGRESSBAR_H
#define GRF_PROGRESSBAR_H
#include <atomic>
#include <mutex>
#include <iostream>
#include <utility>
#include "commons/utility.h"
#include "commons/globals.h"
#include <Rcpp.h>


using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::seconds;


class ProgressBar {

    public:
        ProgressBar (std::string, size_t, bool);

        void increment_progress(size_t value, std::mutex& mutex, std::condition_variable& condition_var) {
            std::unique_lock<std::mutex> lock{mutex};
            progress_ += value;
            condition_var.notify_one();
        }

        void init_time_vars() {
            start_time_ = steady_clock::now();
            last_time_ = steady_clock::now();
        }

        void showProgress(std::mutex& mutex, std::condition_variable& condition_var){
            std::unique_lock<std::mutex> lock{mutex};
            // Wait for message from threads and show output if enough time elapsed
            while (progress_ < max_progress) {
              condition_var.wait(lock);
              elapsed_time_ = duration_cast<seconds>(steady_clock::now() - last_time_);
              if (progress_ > 0 && elapsed_time_.count() > grf::STATUS_INTERVAL){
                double relative_progress = (double) progress_ / (double) max_progress;
                seconds time_from_start = duration_cast<seconds>(steady_clock::now() - start_time_);
                uint remaining_time = (1 / relative_progress - 1) * time_from_start.count();
                if (verbose) {
                    Rcpp::Rcout << operation << " Progress: " << round(100 * relative_progress) << "%. Estimated remaining time: "
                        << grf::beautifyTime(remaining_time) << "." << std::endl;
                }
                last_time_ = steady_clock::now();
              }
            }
        }



    private:
        std::string operation;
        size_t max_progress;
        bool verbose;
        steady_clock::time_point start_time_;
        steady_clock::time_point last_time_;
        seconds elapsed_time_{};
        size_t progress_{0};
};




#endif //GRF_PROGRESSBAR_H
