#ifndef GRF_PROGRESSBAR_H
#define GRF_PROGRESSBAR_H
#include <atomic>
#include <mutex>
#include <iostream>
#include <utility>
#include "commons/utility.h"
#include "commons/globals.h"


using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::seconds;


class ProgressBar {

    public:
        ProgressBar (std::string, size_t, bool);

        void set_progress(size_t value) {
            std::unique_lock<std::mutex> lock{mutex_};
            progress_ += value;
        }

        void set_initial_times() {
            start_time_ = steady_clock::now();
            last_time_ = steady_clock::now();
        }

        void update(size_t value, std::ostream &os = std::cout) {
            set_progress(value);
            write_time_estimate(os);
        }

        void write_time_estimate(std::ostream &os = std::cout) {
            std::unique_lock<std::mutex> lock{mutex_};
            elapsed_time_ = duration_cast<seconds>(steady_clock::now() - last_time_);
            if (progress_ > 0 && elapsed_time_.count() > grf::STATUS_INTERVAL){
                double relative_progress = (double) progress_ / (double) max_progress;
                seconds time_from_start = duration_cast<seconds>(steady_clock::now() - start_time_);
                uint remaining_time = (1 / relative_progress - 1) * time_from_start.count();
                if (verbose) {
                    os << operation << " Progress: " << round(100 * relative_progress) << "%. Estimated remaining time: "
                        << grf::beautifyTime(remaining_time) << "." << std::endl;
                }
                last_time_ = steady_clock::now();
            }
        }


    private:
        std::string operation;
        size_t max_progress;
        bool verbose;
        steady_clock::time_point start_time_;
        steady_clock::time_point last_time_;
        seconds elapsed_time_{};
        std::mutex mutex_;
        size_t progress_{0};
};




#endif //GRF_PROGRESSBAR_H
