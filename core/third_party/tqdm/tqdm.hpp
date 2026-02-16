#pragma once

/*
 *Copyright (c) 2018-2019 <Miguel Raggi> <mraggi@gmail.com>
 *
 *Permission is hereby granted, free of charge, to any person
 *obtaining a copy of this software and associated documentation
 *files (the "Software"), to deal in the Software without
 *restriction, including without limitation the rights to use,
 *copy, modify, merge, publish, distribute, sublicense, and/or sell
 *copies of the Software, and to permit persons to whom the
 *Software is furnished to do so, subject to the following
 *conditions:
 *
 *The above copyright notice and this permission notice shall be
 *included in all copies or substantial portions of the Software.
 *
 *THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *OTHER DEALINGS IN THE SOFTWARE.
 */

#include <chrono>
#include <cmath>
#include <codecvt>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <string>
#include <sstream>
#include <type_traits>

// -------------------- chrono stuff --------------------
namespace tq
{
using index = std::ptrdiff_t; // maybe std::size_t, but I hate unsigned types.
using time_point_t = std::chrono::time_point<std::chrono::steady_clock>;

inline double elapsed_seconds(time_point_t from, time_point_t to)
{
    using seconds = std::chrono::duration<double>;
    return std::chrono::duration_cast<seconds>(to - from).count();
}

class Chronometer
{
public:
    Chronometer() : start_(std::chrono::steady_clock::now()) {}

    double reset()
    {
        auto previous = start_;
        start_ = std::chrono::steady_clock::now();

        return elapsed_seconds(previous, start_);
    }

    double peek() const
    {
        auto now = std::chrono::steady_clock::now();

        return elapsed_seconds(start_, now);
    }

    time_point_t start_;
};

// -------------------- progress_bar --------------------
inline void clamp(double& x, double a, double b)
{
    if (x < a) x = a;
    if (x > b) x = b;
}

class progress_bar
{
public:
    ~progress_bar() {
        if (display_) (*os_) << std::endl;
    }

    void restart()
    {
        chronometer_.reset();
        refresh_.reset();
    }

    void update(
        int iters_done,
        int num_iters
    )
    {
        if (!display_) return;
        if (time_since_refresh() > min_time_per_update_ || iters_done==0 || iters_done==num_iters)
        {
            reset_refresh_timer();
            display(iters_done, num_iters);
        }
        if (iters_done < num_iters) {
            suffix_.str("");
        }
    }

    void set_ostream(std::ostream& os) { os_ = &os; }
    void set_prefix(std::string s) { prefix_ = std::move(s); }
    void set_bar_size(int size) { bar_size_ = size; }
    void set_min_update_time(double time) { min_time_per_update_ = time; }
    void set_display(bool b) { display_ = b; }
    void set_bar_symbol(std::string s) { bar_symbol_ = std::move(s); }

    template <class T>
    progress_bar& operator<<(const T& t)
    {
        suffix_ << t;
        return *this;
    }

    double elapsed_time() const { return chronometer_.peek(); }

private:
    void display(int iters_done, int num_iters)
    {
        double progress = iters_done / (num_iters + 1e-9);
        clamp(progress,0.0,1.0);

        auto flags = os_->flags();

        double t = chronometer_.peek();
        double eta = std::max(t/progress - t, 0.0);

        std::stringstream line_ss;

        const auto get_ss_length = [](auto& ss) {
            auto pos = ss.tellg();
            ss.seekg(0, ss.end);
            const auto size = ss.tellg();
            ss.seekg(pos);
            return size;
        };

        line_ss << '\r' << prefix_
            << std::fixed << std::setprecision(0) << std::setw(3)
            << 100*progress
            << '%'
            ;

        print_bar(line_ss, progress);

        const auto line_ss_initial_bytesize = get_ss_length(line_ss);

        line_ss << ' ' << iters_done << '/' << num_iters << ' ';

        const auto print_time = [&](auto s) {
            const int s_hours = s / 3600.;
            const int s_minutes = (s - s_hours * 3600) / 60.;
            const int s_seconds = s - s_hours * 3600 - s_minutes * 60;
            line_ss << std::setfill('0') << std::setw(2) << s_hours << ':'
                << std::setfill('0') << std::setw(2) << s_minutes << ':'
                << std::setfill('0') << std::setw(2) << s_seconds
                ;
        };

        line_ss << '[';
        print_time(t);
        line_ss << '<';

        if (std::isinf(eta) || std::isnan(eta)) {
            line_ss << '?';
        } else {
            print_time(eta);
        }

        line_ss << ", ";

        if (iters_done == 0) {
            line_ss << '?';
        } else {
            line_ss << std::setprecision(2)
                    << iters_done / t;
        }

        line_ss << "it/s]" << suffix_.str();

        const auto line_size = get_ss_length(line_ss) - line_ss_initial_bytesize;
        max_chars_ = std::max<index>(max_chars_, line_size);

        // NOTE: MUST append the spaces instead of doing << twice with two separate strings.
        // No idea why but if we do not do this, then the buffer isn't overwritten with spaces at the end.
        line_ss << std::string(max_chars_-line_size, ' ');
        (*os_) << line_ss.str() << std::flush;

        os_->flags(flags);
    }

    void print_bar(std::stringstream& ss, double filled) const
    {
        auto num_filled = static_cast<index>(std::round(filled*bar_size_));
        ss << '|';
        for (int i = 0; i < num_filled; ++i) {
            ss << bar_symbol_;
        }
        ss << std::string(bar_size_ - num_filled, ' ') << '|';
    }

    double time_since_refresh() const { return refresh_.peek(); }
    void reset_refresh_timer() { refresh_.reset(); }

    Chronometer chronometer_{};
    Chronometer refresh_{};
    double min_time_per_update_{1e-1}; // found experimentally
    bool display_{true};

    std::ostream* os_ = nullptr;

    index bar_size_{10};
    index max_chars_{0};

    std::string prefix_{};
    std::stringstream suffix_{};
    std::string bar_symbol_ = "\033[1;32m\u2588\033[0m";
};


// -------------------- iter_wrapper --------------------

template <class ForwardIter, class Parent>
class iter_wrapper
{
public:
    using iterator_category = typename ForwardIter::iterator_category;
    using value_type = typename ForwardIter::value_type;
    using difference_type = typename ForwardIter::difference_type;
    using pointer = typename ForwardIter::pointer;
    using reference = typename ForwardIter::reference;

    iter_wrapper(ForwardIter it, Parent* parent)
        : current_(it), parent_(parent)
    {}

    auto operator*() { return *current_; }

    void operator++() { ++current_; }

    template <class Other>
    bool operator!=(const Other& other) const
    {
        parent_->update(); // here and not in ++ because I need to run update before first advancement!
        return current_ != other;
    }

    bool operator!=(const iter_wrapper<ForwardIter, Parent>& other) const
    {
        parent_->update(); // here and not in ++ because I need to run update before first advancement!
        return current_ != other.current_;
    }

    const ForwardIter& get() const { return current_; }

private:
    friend Parent;
    ForwardIter current_;
    Parent* parent_;
};

// -------------------- tqdm_for_lvalues --------------------

template <class ForwardIter, class EndIter=ForwardIter>
class tqdm_for_lvalues
{
public:
    using this_t = tqdm_for_lvalues<ForwardIter,EndIter>;
    using iterator = iter_wrapper<ForwardIter, this_t>;
    using value_type = typename ForwardIter::value_type;
    using size_type = index;
    using difference_type = index;

    tqdm_for_lvalues(ForwardIter begin, EndIter end)
        : first_(begin, this)
        , last_(end)
        , num_iters_(std::distance(begin, end))
    {}

    tqdm_for_lvalues(ForwardIter begin, EndIter end, index total)
        : first_(begin, this), last_(end), num_iters_(total)
    {}

    template <class Container>
    explicit tqdm_for_lvalues(Container& C)
        : first_(C.begin(), this), last_(C.end()), num_iters_(C.size())
    {}

    template <class Container>
    explicit tqdm_for_lvalues(const Container& C)
        : first_(C.begin(), this), last_(C.end()), num_iters_(C.size())
    {}

    tqdm_for_lvalues(const tqdm_for_lvalues&) = delete;
    tqdm_for_lvalues(tqdm_for_lvalues&&) = delete;
    tqdm_for_lvalues& operator=(tqdm_for_lvalues&&) = delete;
    tqdm_for_lvalues& operator=(const tqdm_for_lvalues&) = delete;

    template <class Container>
    tqdm_for_lvalues(Container&&) = delete; // prevent misuse!

    ~tqdm_for_lvalues() {
        bar_.set_min_update_time(0); // force printing in next line
        bar_.update(iters_done_-1, num_iters_);
    }

    iterator begin()
    {
        bar_.restart();
        iters_done_ = 0;
        return first_;
    }

    EndIter end() const { return last_; }

    void update()
    {
        bar_.update(iters_done_, num_iters_);
        ++iters_done_;
    }

    // NOTE: possibly invalidates other members.
    // Safe usage is to call .begin(), .end() afterwards to reset everything.
    void set_range(ForwardIter begin, EndIter end)
    {
        first_ = iterator(begin, this);
        last_ = end;
        num_iters_ = std::distance(begin, end);
    }
    void set_ostream(std::ostream& os) { bar_.set_ostream(os); }
    void set_prefix(std::string s) { bar_.set_prefix(std::move(s)); }
    void set_bar_size(int size) { bar_.set_bar_size(size); }
    void set_min_update_time(double time) { bar_.set_min_update_time(time); }
    void set_display(bool b) { bar_.set_display(b); }
    void set_symbol(std::string s) { bar_.set_bar_symbol(std::move(s)); }

    template <class T>
    tqdm_for_lvalues& operator<<(const T& t)
    {
        bar_ << t;
        return *this;
    }

    void manually_set_progress(int to)
    {
        iters_done_ = to;
    }

private:
    iterator first_;
    EndIter last_;
    index num_iters_;
    index iters_done_{0};
    progress_bar bar_;
};

template <class Container>
tqdm_for_lvalues(Container&)->tqdm_for_lvalues<typename Container::iterator>;

template <class Container>
tqdm_for_lvalues(const Container&)
  ->tqdm_for_lvalues<typename Container::const_iterator>;

// -------------------- tqdm_for_rvalues --------------------

template <class Container>
class tqdm_for_rvalues
{
public:
    using iterator = typename Container::iterator;
    using const_iterator = typename Container::const_iterator;
    using value_type = typename Container::value_type;

    explicit tqdm_for_rvalues(Container&& C)
        : C_(std::forward<Container>(C)), tqdm_(C_)
    {}

    auto begin() { return tqdm_.begin(); }

    auto end() { return tqdm_.end(); }

    void update() { return tqdm_.update(); }

    void set_range(iterator begin, iterator end) { tqdm_.set_range(begin, end); }
    void set_ostream(std::ostream& os) { tqdm_.set_ostream(os); }
    void set_prefix(std::string s) { tqdm_.set_prefix(std::move(s)); }
    void set_bar_size(int size) { tqdm_.set_bar_size(size); }
    void set_min_update_time(double time) { tqdm_.set_min_update_time(time); }
    void set_display(bool b) { tqdm_.set_display(b) ; }

    template <class T>
    auto& operator<<(const T& t)
    {
        return tqdm_ << t;
    }

    void advance(index amount) { tqdm_.advance(amount); }

    void manually_set_progress(double to)
    {
        tqdm_.manually_set_progress(to);
    }

private:
    Container C_;
    tqdm_for_lvalues<iterator> tqdm_;
};

template <class Container>
tqdm_for_rvalues(Container &&)->tqdm_for_rvalues<Container>;

// -------------------- tqdm --------------------
template <class ForwardIter>
inline auto tqdm(const ForwardIter& first, const ForwardIter& last)
{
    return tqdm_for_lvalues(first, last);
}

template <class ForwardIter>
inline auto tqdm(const ForwardIter& first, const ForwardIter& last, index total)
{
    return tqdm_for_lvalues(first, last, total);
}

template <class Container>
inline auto tqdm(const Container& C)
{
    return tqdm_for_lvalues(C);
}

template <class Container>
inline auto tqdm(Container& C)
{
    return tqdm_for_lvalues(C);
}

template <class Container>
inline auto tqdm(Container&& C)
{
    return tqdm_for_rvalues(std::forward<Container>(C));
}

// -------------------- int_iterator --------------------

template <class IntType>
class int_iterator
{
public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = IntType;
    using difference_type = IntType;
    using pointer = IntType*;
    using reference = IntType&;

    explicit int_iterator(IntType val) : value_(val) {}

    IntType& operator*() { return value_; }

    int_iterator& operator++()
    {
        ++value_;
        return *this;
    }
    int_iterator& operator--()
    {
        --value_;
        return *this;
    }

    int_iterator& operator+=(difference_type d)
    {
        value_ += d;
        return *this;
    }

    difference_type operator-(const int_iterator& other) const
    {
        return value_ - other.value_;
    }

    bool operator!=(const int_iterator& other) const
    {
        return value_ != other.value_;
    }

private:
    IntType value_;
};

// -------------------- range --------------------
template <class IntType>
class range
{
public:
    using iterator = int_iterator<IntType>;
    using const_iterator = iterator;
    using value_type = IntType;

    range(IntType first, IntType last) : first_(first), last_(last) {}
    explicit range(IntType last) : first_(0), last_(last) {}

    iterator begin() const { return first_; }
    iterator end() const { return last_; }
    index size() const { return last_ - first_; }

private:
    iterator first_;
    iterator last_;
};

template <class IntType>
auto trange(IntType first, IntType last)
{
    return tqdm(range(first, last));
}

template <class IntType>
auto trange(IntType last)
{
    return tqdm(range(last));
}

// -------------------- timing_iterator --------------------

struct timing_iterator_end_sentinel
{
public:
    explicit timing_iterator_end_sentinel(double num_seconds) : num_seconds_(num_seconds) {}

    double num_seconds_;
};

class timing_iterator
{
public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = double;
    using difference_type = double;
    using pointer = double*;
    using reference = double&;

    double& operator*() const
    {
        workaround_ = chrono_.peek();
        return workaround_;
    }

    timing_iterator& operator++()
    {
        return *this;
    }

    bool operator!=(const timing_iterator_end_sentinel& other) const
    {
        return chrono_.peek() < other.num_seconds_;
    }

private:
    tq::Chronometer chrono_;
    mutable double workaround_;
};

// -------------------- timer -------------------
struct timer
{
public:
    using iterator = timing_iterator;
    using end_iterator = timing_iterator_end_sentinel;
    using const_iterator = iterator;
    using value_type = double;

    explicit timer(double num_seconds) : num_seconds(num_seconds)
    {
    }

    iterator begin() const { return iterator(); }
    end_iterator end() const { return end_iterator(num_seconds); }

    double num_seconds;
};

class tqdm_timer
{
public:
    using iterator = iter_wrapper<timing_iterator,tqdm_timer>;
    using end_iterator = timer::end_iterator;
    using value_type = typename timing_iterator::value_type;
    using size_type = index;
    using difference_type = index;

    explicit tqdm_timer(double num_seconds)
        : num_seconds_(num_seconds)
    {}

    tqdm_timer(const tqdm_timer&) = delete;
    tqdm_timer(tqdm_timer&&) = delete;
    tqdm_timer& operator=(tqdm_timer&&) = delete;
    tqdm_timer& operator=(const tqdm_timer&) = delete;

    template <class Container>
    tqdm_timer(Container&&) = delete; // prevent misuse!

    iterator begin()
    {
        bar_.restart();
        return iterator(timing_iterator(), this);
    }

    end_iterator end() const { return end_iterator(num_seconds_); }

    void update()
    {
        double t = bar_.elapsed_time();

        bar_.update(t, num_seconds_);
    }

    void set_ostream(std::ostream& os) { bar_.set_ostream(os); }
    void set_prefix(std::string s) { bar_.set_prefix(std::move(s)); }
    void set_bar_size(int size) { bar_.set_bar_size(size); }
    void set_min_update_time(double time) { bar_.set_min_update_time(time); }

    template <class T>
    tqdm_timer& operator<<(const T& t)
    {
        bar_ << t;
        return *this;
    }


private:
    double num_seconds_;
    progress_bar bar_;
};

inline auto tqdm(timer t)
{
    return tqdm_timer(t.num_seconds);
}

template <class IntType>
using progress_bar_int_t = std::decay_t<decltype(trange(std::declval<IntType>()))>;
using progress_bar_t = progress_bar_int_t<int>;

} // namespace tq
