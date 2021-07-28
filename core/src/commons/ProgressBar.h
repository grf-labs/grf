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

        void set_progress(float value_pct, size_t value_int) {
            std::unique_lock<std::mutex> lock{mutex_};
            progress_pct_ = value_pct;
            progress_ += value_int;
        }

        void set_initial_times() {
            start_time_ = steady_clock::now();
            last_time_ = steady_clock::now();
        }

        void set_bar_width(size_t width) {
            std::unique_lock<std::mutex> lock{mutex_};
            bar_width_ = width;
        }

        void fill_bar_progress_with(const std::string &chars) {
            std::unique_lock<std::mutex> lock{mutex_};
            fill_ = chars;
        }

        void fill_bar_remainder_with(const std::string &chars) {
            std::unique_lock<std::mutex> lock{mutex_};
            remainder_ = chars;
        }

        void set_status_text(const std::string &status) {
            std::unique_lock<std::mutex> lock{mutex_};
            status_text_ = status;
        }

        void update(float value_pct, size_t value_int, std::ostream &os = std::cout) {
            set_progress(value_pct, value_int);
//            write_progress(os);
            write_time_estimate(os);
        }

        void write_time_estimate(std::ostream &os = Rcpp::Rcout) {
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

        void write_progress(std::ostream &os = std::cout) {
            std::unique_lock<std::mutex> lock{mutex_};

            // No need to write once progress is 100%
            if (progress_pct_ > 100.0f) return;

            // Move cursor to the first position on the same line and flush
            os << "\r" << std::flush;

            // Start bar
            os << "[";

            const auto completed = static_cast<size_t>(progress_pct_ * static_cast<float>(bar_width_) / 100.0);
            for (size_t i = 0; i < bar_width_; ++i) {
                if (i <= completed)
                    os << fill_;
                else
                    os << remainder_;
            }

            // End bar
            os << "]";

            // Write progress percentage
            os << " " << std::min(static_cast<size_t>(progress_pct_), size_t(100)) << "%";

            // Write status text
            os << " " << status_text_;
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
        float progress_pct_{0.0f};
        size_t bar_width_{60};
        std::string fill_{"#"}, remainder_{" "}, status_text_{""};
};




#endif //GRF_PROGRESSBAR_H
