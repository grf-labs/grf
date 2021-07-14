#include <limits>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <vector>
#include <initializer_list>
#include <numeric>
#include <iostream>
#include <random>
#include <type_traits>

// MSVC compiler directives from:
// https://github.com/llvm-mirror/libcxx/blob/master/include/__config#L180
// https://github.com/llvm-mirror/libcxx/blob/master/include/support/win32/limits_msvc_win32.h#L25
#ifdef _MSC_VER
#define _LIBCPP_COMPILER_MSVC
#define __CHAR_BIT__       CHAR_BIT
#endif

// https://github.com/llvm-mirror/libcxx/blob/master/include/__config#L278
#if defined(_WIN32)
#  if (defined(_M_AMD64) || defined(__x86_64__)) || (defined(_M_ARM) || defined(__arm__))
#    define _LIBCPP_HAS_BITSCAN64
#  endif
#  endif

#ifndef _LIBGRF_RANDOM
#define _LIBGRF_RANDOM

namespace nonstd {


// Precondition:  __x != 0
    inline
    unsigned __clz(unsigned __x) {
#ifndef _LIBCPP_COMPILER_MSVC
        return static_cast<unsigned>(__builtin_clz(__x));
#else
        static_assert(sizeof(unsigned) == sizeof(unsigned long), "");
  static_assert(sizeof(unsigned long) == 4, "");
  unsigned long where;
  // Search from LSB to MSB for first set bit.
  // Returns zero if no set bit is found.
  if (_BitScanReverse(&where, __x))
    return 31 - where;
  return 32; // Undefined Behavior.
#endif
    }

    inline
    unsigned long __clz(unsigned long __x) {
#ifndef _LIBCPP_COMPILER_MSVC
        return static_cast<unsigned long>(__builtin_clzl(__x));
#else
        static_assert(sizeof(unsigned) == sizeof(unsigned long), "");
    return __clz(static_cast<unsigned>(__x));
#endif
    }

    inline
    unsigned long long __clz(unsigned long long __x) {
#ifndef _LIBCPP_COMPILER_MSVC
        return static_cast<unsigned long long>(__builtin_clzll(__x));
#else
        unsigned long where;
// BitScanReverse scans from MSB to LSB for first set bit.
// Returns 0 if no set bit is found.
#if defined(_LIBCPP_HAS_BITSCAN64)
  if (_BitScanReverse64(&where, __x))
    return static_cast<int>(63 - where);
#else
  // Scan the high 32 bits.
  if (_BitScanReverse(&where, static_cast<unsigned long>(__x >> 32)))
    return 63 - (where + 32); // Create a bit offset from the MSB.
  // Scan the low 32 bits.
  if (_BitScanReverse(&where, static_cast<unsigned long>(__x)))
    return 63 - where;
#endif
  return 64; // Undefined Behavior.
#endif // _LIBCPP_COMPILER_MSVC
    }


// __independent_bits_engine

    template<unsigned long long _Xp, size_t _Rp>
    struct __log2_imp {
        static const size_t value = _Xp & ((unsigned long long) (1) << _Rp) ? _Rp
                                                                            : __log2_imp<_Xp, _Rp - 1>::value;
    };

    template<unsigned long long _Xp>
    struct __log2_imp<_Xp, 0> {
        static const size_t value = 0;
    };

    template<size_t _Rp>
    struct __log2_imp<0, _Rp> {
        static const size_t value = _Rp + 1;
    };

    template<class _UIntType, _UIntType _Xp>
    struct __log2 {
        static const size_t value = __log2_imp<_Xp,
                sizeof(_UIntType) * __CHAR_BIT__ - 1>::value;
    };

    template<class _Engine, class _UIntType>
    class __independent_bits_engine {
    public:
        // types
        typedef _UIntType result_type;

    private:
        typedef typename _Engine::result_type _Engine_result_type;
        typedef typename std::conditional
                <
                        sizeof(_Engine_result_type) <= sizeof(result_type),
                        result_type,
                        _Engine_result_type
                >::type _Working_result_type;

        _Engine &__e_;
        size_t __w_;
        size_t __w0_;
        size_t __n_;
        size_t __n0_;
        _Working_result_type __y0_;
        _Working_result_type __y1_;
        _Engine_result_type __mask0_;
        _Engine_result_type __mask1_;

#ifdef _LIBCPP_CXX03_LANG
        static const _Working_result_type _Rp = _Engine::_Max - _Engine::_Min
                                          + _Working_result_type(1);
#else
        static const _Working_result_type _Rp = _Engine::max() - _Engine::min()
                                                                  + _Working_result_type(1);
#endif
        static const size_t __m = __log2<_Working_result_type, _Rp>::value;
        static const size_t _WDt = std::numeric_limits<_Working_result_type>::digits;
        static const size_t _EDt = std::numeric_limits<_Engine_result_type>::digits;

    public:
        // constructors and seeding functions
        __independent_bits_engine(_Engine &__e, size_t __w);

        // generating functions
        result_type operator()() { return __eval(std::integral_constant<bool, _Rp != 0>()); }

    private:
        result_type __eval(std::false_type);

        result_type __eval(std::true_type);
    };

    template<class _Engine, class _UIntType>
    __independent_bits_engine<_Engine, _UIntType>
    ::__independent_bits_engine(_Engine &__e, size_t __w)
            : __e_(__e),
              __w_(__w) {
        __n_ = __w_ / __m + (__w_ % __m != 0);
        __w0_ = __w_ / __n_;
        if (_Rp == 0)
            __y0_ = _Rp;
        else if (__w0_ < _WDt)
            __y0_ = (_Rp >> __w0_) << __w0_;
        else
            __y0_ = 0;
        if (_Rp - __y0_ > __y0_ / __n_) {
            ++__n_;
            __w0_ = __w_ / __n_;
            if (__w0_ < _WDt)
                __y0_ = (_Rp >> __w0_) << __w0_;
            else
                __y0_ = 0;
        }
        __n0_ = __n_ - __w_ % __n_;
        if (__w0_ < _WDt - 1)
            __y1_ = (_Rp >> (__w0_ + 1)) << (__w0_ + 1);
        else
            __y1_ = 0;
        __mask0_ = __w0_ > 0 ? _Engine_result_type(~0) >> (_EDt - __w0_) :
                   _Engine_result_type(0);
        __mask1_ = __w0_ < _EDt - 1 ?
                   _Engine_result_type(~0) >> (_EDt - (__w0_ + 1)) :
                   _Engine_result_type(~0);
    }

    template<class _Engine, class _UIntType>
    inline
    _UIntType
    __independent_bits_engine<_Engine, _UIntType>::__eval(std::false_type) {
        return static_cast<result_type>(__e_() & __mask0_);
    }

    template<class _Engine, class _UIntType>
    _UIntType
    __independent_bits_engine<_Engine, _UIntType>::__eval(std::true_type) {
        const size_t _WRt = std::numeric_limits<result_type>::digits;
        result_type _Sp = 0;
        for (size_t __k = 0; __k < __n0_; ++__k) {
            _Engine_result_type __u;
            do {
                __u = __e_() - _Engine::min();
            } while (__u >= __y0_);
            if (__w0_ < _WRt)
                _Sp <<= __w0_;
            else
                _Sp = 0;
            _Sp += __u & __mask0_;
        }
        for (size_t __k = __n0_; __k < __n_; ++__k) {
            _Engine_result_type __u;
            do {
                __u = __e_() - _Engine::min();
            } while (__u >= __y1_);
            if (__w0_ < _WRt - 1)
                _Sp <<= __w0_ + 1;
            else
                _Sp = 0;
            _Sp += __u & __mask1_;
        }
        return _Sp;
    }


// uniform_int_distribution

    template<class _IntType = int>
    class uniform_int_distribution {
    public:
        // types
        typedef _IntType result_type;

        class param_type {
            result_type __a_;
            result_type __b_;
        public:
            typedef uniform_int_distribution distribution_type;

            explicit param_type(result_type __a = 0,
                                result_type __b = std::numeric_limits<result_type>::max())
                    : __a_(__a), __b_(__b) {}

            result_type a() const { return __a_; }

            result_type b() const { return __b_; }

            friend bool operator==(const param_type &__x, const param_type &__y) {
                return __x.__a_ == __y.__a_ && __x.__b_ == __y.__b_;
            }

            friend bool operator!=(const param_type &__x, const param_type &__y) { return !(__x == __y); }
        };

    private:
        param_type __p_;

    public:
        // constructors and reset functions
        explicit uniform_int_distribution(result_type __a = 0,
                                          result_type __b = std::numeric_limits<result_type>::max())
                : __p_(param_type(__a, __b)) {}

        explicit uniform_int_distribution(const param_type &__p) : __p_(__p) {}

        void reset() {}

        // generating functions
        template<class _URNG>
        result_type operator()(_URNG &__g) { return (*this)(__g, __p_); }

        template<class _URNG>
        result_type operator()(_URNG &__g, const param_type &__p);

        // property functions
        result_type a() const { return __p_.a(); }

        result_type b() const { return __p_.b(); }

        param_type param() const { return __p_; }

        void param(const param_type &__p) { __p_ = __p; }

        result_type min() const { return a(); }

        result_type max() const { return b(); }

        friend bool operator==(const uniform_int_distribution &__x,
                               const uniform_int_distribution &__y) { return __x.__p_ == __y.__p_; }

        friend bool operator!=(const uniform_int_distribution &__x,
                               const uniform_int_distribution &__y) { return !(__x == __y); }
    };

    template<class _IntType>
    template<class _URNG>
    typename uniform_int_distribution<_IntType>::result_type
    uniform_int_distribution<_IntType>::operator()(_URNG &__g, const param_type &__p) {
        typedef typename std::conditional<sizeof(result_type) <= sizeof(uint32_t),
                uint32_t, uint64_t>::type _UIntType;
        const _UIntType _Rp = __p.b() - __p.a() + _UIntType(1);
        if (_Rp == 1)
            return __p.a();
        const size_t _Dt = std::numeric_limits<_UIntType>::digits;
        typedef __independent_bits_engine<_URNG, _UIntType> _Eng;
        if (_Rp == 0)
            return static_cast<result_type>(_Eng(__g, _Dt)());
        size_t __w = _Dt - __clz(_Rp) - 1;
        if ((_Rp & (std::numeric_limits<_UIntType>::max() >> (_Dt - __w))) != 0)
            ++__w;
        _Eng __e(__g, __w);
        _UIntType __u;
        do {
            __u = __e();
        } while (__u >= _Rp);
        return static_cast<result_type>(__u + __p.a());
    }

#if _LIBCPP_STD_VER <= 14 || defined(_LIBCPP_ENABLE_CXX17_REMOVED_RANDOM_SHUFFLE) \
 || defined(_LIBCPP_BUILDING_LIBRARY)

    class  __rs_default;

     __rs_default __rs_get();

    class  __rs_default {
        static unsigned __c_;

        __rs_default();

    public:
        typedef uint_fast32_t result_type;

        static const result_type _Min = 0;
        static const result_type _Max = 0xFFFFFFFF;

        __rs_default(const __rs_default &);

        ~__rs_default();

        result_type operator()();

        static result_type min() { return _Min; }

        static result_type max() { return _Max; }

        friend  __rs_default __rs_get();
    };

     __rs_default __rs_get();

    template<class _RandomAccessIterator>
    void
    random_shuffle(_RandomAccessIterator __first, _RandomAccessIterator __last) {
        typedef typename std::iterator_traits<_RandomAccessIterator>::difference_type difference_type;
        typedef uniform_int_distribution<ptrdiff_t> _Dp;
        typedef typename _Dp::param_type _Pp;
        difference_type __d = __last - __first;
        if (__d > 1) {
            _Dp __uid;
            __rs_default __g = __rs_get();
            for (--__last, --__d; __first < __last; ++__first, --__d) {
                difference_type __i = __uid(__g, _Pp(0, __d));
                if (__i != difference_type(0))
                    swap(*__first, *(__first + __i));
            }
        }
    }

    template<class _RandomAccessIterator, class _RandomNumberGenerator>
    void
    random_shuffle(_RandomAccessIterator __first, _RandomAccessIterator __last,
#ifndef _LIBCPP_CXX03_LANG
                   _RandomNumberGenerator &&__rand)
#else
    _RandomNumberGenerator& __rand)
#endif
    {
        typedef typename std::iterator_traits<_RandomAccessIterator>::difference_type difference_type;
        difference_type __d = __last - __first;
        if (__d > 1) {
            for (--__last; __first < __last; ++__first, --__d) {
                difference_type __i = __rand(__d);
                swap(*__first, *(__first + __i));
            }
        }
    }

#endif

    template<class _PopulationIterator, class _SampleIterator, class _Distance,
            class _UniformRandomNumberGenerator>

    _SampleIterator __sample(_PopulationIterator __first,
                             _PopulationIterator __last, _SampleIterator __output_iter,
                             _Distance __n,
                             _UniformRandomNumberGenerator &__g,
                             std::input_iterator_tag) {

        _Distance __k = 0;
        for (; __first != __last && __k < __n; ++__first, (void) ++__k)
            __output_iter[__k] = *__first;
        _Distance __sz = __k;
        for (; __first != __last; ++__first, (void) ++__k) {
            _Distance __r = uniform_int_distribution<_Distance>(0, __k)(__g);
            if (__r < __sz)
                __output_iter[__r] = *__first;
        }
        return __output_iter + std::min(__n, __k);
    }



// generate_canonical

    template<class _RealType, size_t __bits, class _URNG>
    _RealType
    __generate_canonical(_URNG &__g) {
        const size_t _Dt = std::numeric_limits<_RealType>::digits;
        const size_t __b = _Dt < __bits ? _Dt : __bits;
#ifdef _LIBCPP_CXX03_LANG
        const size_t __logR = __log2<uint64_t, _URNG::_Max - _URNG::_Min + uint64_t(1)>::value;
#else
        const size_t __logR = __log2<uint64_t, _URNG::max() - _URNG::min() + uint64_t(1)>::value;
#endif
        const size_t __k = __b / __logR + (__b % __logR != 0) + (__b == 0);
        const _RealType _Rp = _URNG::max() - _URNG::min() + _RealType(1);
        _RealType __base = _Rp;
        _RealType _Sp = __g() - _URNG::min();
        for (size_t __i = 1; __i < __k; ++__i, __base *= _Rp)
            _Sp += (__g() - _URNG::min()) * __base;
        return _Sp / __base;
    }




// poisson distribution

    template<class _IntType = int>
    class  poisson_distribution {
    public:
        // types
        typedef _IntType result_type;

        class  param_type {
            double __mean_;
            double __s_;
            double __d_;
            double __l_;
            double __omega_;
            double __c0_;
            double __c1_;
            double __c2_;
            double __c3_;
            double __c_;

        public:
            typedef poisson_distribution distribution_type;

            explicit param_type(double __mean = 1.0);


            double mean() const { return __mean_; }

            friend
            bool operator==(const param_type &__x, const param_type &__y) { return __x.__mean_ == __y.__mean_; }

            friend
            bool operator!=(const param_type &__x, const param_type &__y) { return !(__x == __y); }

            friend class poisson_distribution;
        };

    private:
        param_type __p_;

    public:
        // constructors and reset functions

        explicit poisson_distribution(double __mean = 1.0) : __p_(__mean) {}


        explicit poisson_distribution(const param_type &__p) : __p_(__p) {}


        void reset() {}

        // generating functions
        template<class _URNG>

        result_type operator()(_URNG &__g) { return (*this)(__g, __p_); }

        template<class _URNG>
        result_type operator()(_URNG &__g, const param_type &__p);

        // property functions

        double mean() const { return __p_.mean(); }


        param_type param() const { return __p_; }


        void param(const param_type &__p) { __p_ = __p; }


        result_type min() const { return 0; }


        result_type max() const { return std::numeric_limits<result_type>::max(); }

        friend
        bool operator==(const poisson_distribution &__x,
                        const poisson_distribution &__y) { return __x.__p_ == __y.__p_; }

        friend
        bool operator!=(const poisson_distribution &__x,
                        const poisson_distribution &__y) { return !(__x == __y); }
    };

//#define std std::_LIBCPP_NAMESPACE

    template<class _IntType>
    poisson_distribution<_IntType>::param_type::param_type(double __mean)
            : __mean_(__mean) {
        if (__mean_ < 10) {
            __s_ = 0;
            __d_ = 0;
            __l_ = std::exp(-__mean_);
            __omega_ = 0;
            __c3_ = 0;
            __c2_ = 0;
            __c1_ = 0;
            __c0_ = 0;
            __c_ = 0;
        } else {
            __s_ = std::sqrt(__mean_);
            __d_ = 6 * __mean_ * __mean_;
            __l_ = static_cast<result_type>(__mean_ - 1.1484);
            __omega_ = .3989423 / __s_;
            double __b1_ = .4166667E-1 / __mean_;
            double __b2_ = .3 * __b1_ * __b1_;
            __c3_ = .1428571 * __b1_ * __b2_;
            __c2_ = __b2_ - 15. * __c3_;
            __c1_ = __b1_ - 6. * __b2_ + 45. * __c3_;
            __c0_ = 1. - __b1_ + 3. * __b2_ - 15. * __c3_;
            __c_ = .1069 / __mean_;
        }
    }


// uniform_real_distribution

    template<class _RealType = double>
    class  uniform_real_distribution {
    public:
        // types
        typedef _RealType result_type;

        class  param_type {
            result_type __a_;
            result_type __b_;
        public:
            typedef uniform_real_distribution distribution_type;


            explicit param_type(result_type __a = 0,
                                result_type __b = 1)
                    : __a_(__a), __b_(__b) {}


            result_type a() const { return __a_; }


            result_type b() const { return __b_; }

            friend
            bool operator==(const param_type &__x, const param_type &__y) {
                return __x.__a_ == __y.__a_ && __x.__b_ == __y.__b_;
            }

            friend
            bool operator!=(const param_type &__x, const param_type &__y) { return !(__x == __y); }
        };

    private:
        param_type __p_;

    public:
        // constructors and reset functions

        explicit uniform_real_distribution(result_type __a = 0, result_type __b = 1)
                : __p_(param_type(__a, __b)) {}


        explicit uniform_real_distribution(const param_type &__p) : __p_(__p) {}


        void reset() {}

        // generating functions
        template<class _URNG>

        result_type operator()(_URNG &__g) { return (*this)(__g, __p_); }

        template<class _URNG>
         result_type operator()(_URNG &__g, const param_type &__p);

        // property functions

        result_type a() const { return __p_.a(); }


        result_type b() const { return __p_.b(); }


        param_type param() const { return __p_; }


        void param(const param_type &__p) { __p_ = __p; }


        result_type min() const { return a(); }


        result_type max() const { return b(); }

        friend
        bool operator==(const uniform_real_distribution &__x,
                        const uniform_real_distribution &__y) { return __x.__p_ == __y.__p_; }

        friend
        bool operator!=(const uniform_real_distribution &__x,
                        const uniform_real_distribution &__y) { return !(__x == __y); }
    };

    template<class _RealType>
    template<class _URNG>
    inline
    typename uniform_real_distribution<_RealType>::result_type
    uniform_real_distribution<_RealType>::operator()(_URNG &__g, const param_type &__p) {
        return (__p.b() - __p.a())
               * __generate_canonical<_RealType, std::numeric_limits<_RealType>::digits>(__g)
               + __p.a();
    }


// normal_distribution

    template<class _RealType = double>
    class  normal_distribution {
    public:
        // types
        typedef _RealType result_type;

        class  param_type {
            result_type __mean_;
            result_type __stddev_;
        public:
            typedef normal_distribution distribution_type;


            explicit param_type(result_type __mean = 0, result_type __stddev = 1)
                    : __mean_(__mean), __stddev_(__stddev) {}


            result_type mean() const { return __mean_; }


            result_type stddev() const { return __stddev_; }

            friend
            bool operator==(const param_type &__x, const param_type &__y) {
                return __x.__mean_ == __y.__mean_ && __x.__stddev_ == __y.__stddev_;
            }

            friend
            bool operator!=(const param_type &__x, const param_type &__y) { return !(__x == __y); }
        };

    private:
        param_type __p_;
        result_type _V_;
        bool _V_hot_;

    public:
        // constructors and reset functions

        explicit normal_distribution(result_type __mean = 0, result_type __stddev = 1)
                : __p_(param_type(__mean, __stddev)), _V_hot_(false) {}


        explicit normal_distribution(const param_type &__p)
                : __p_(__p), _V_hot_(false) {}


        void reset() { _V_hot_ = false; }

        // generating functions
        template<class _URNG>

        result_type operator()(_URNG &__g) { return (*this)(__g, __p_); }

        template<class _URNG>
        result_type operator()(_URNG &__g, const param_type &__p);

        // property functions

        result_type mean() const { return __p_.mean(); }


        result_type stddev() const { return __p_.stddev(); }


        param_type param() const { return __p_; }


        void param(const param_type &__p) { __p_ = __p; }


        result_type min() const { return -std::numeric_limits<result_type>::infinity(); }


        result_type max() const { return std::numeric_limits<result_type>::infinity(); }

        friend
        bool operator==(const normal_distribution &__x,
                        const normal_distribution &__y) {
            return __x.__p_ == __y.__p_ && __x._V_hot_ == __y._V_hot_ &&
                   (!__x._V_hot_ || __x._V_ == __y._V_);
        }

        friend
        bool operator!=(const normal_distribution &__x,
                        const normal_distribution &__y) { return !(__x == __y); }

    };

    template<class _RealType>
    template<class _URNG>
    _RealType
    normal_distribution<_RealType>::operator()(_URNG &__g, const param_type &__p) {
        result_type _Up;
        if (_V_hot_) {
            _V_hot_ = false;
            _Up = _V_;
        } else {
            uniform_real_distribution<result_type> _Uni(-1, 1);
            result_type __u;
            result_type __v;
            result_type __s;
            do {
                __u = _Uni(__g);
                __v = _Uni(__g);
                __s = __u * __u + __v * __v;
            } while (__s > 1 || __s == 0);
            result_type _Fp = std::sqrt(-2 * std::log(__s) / __s);
            _V_ = __v * _Fp;
            _V_hot_ = true;
            _Up = __u * _Fp;
        }
        return _Up * __p.stddev() + __p.mean();
    }




// exponential_distribution

    template<class _RealType = double>
    class  exponential_distribution {
    public:
        // types
        typedef _RealType result_type;

        class  param_type {
            result_type __lambda_;
        public:
            typedef exponential_distribution distribution_type;


            explicit param_type(result_type __lambda = 1) : __lambda_(__lambda) {}


            result_type lambda() const { return __lambda_; }

            friend
            bool operator==(const param_type &__x, const param_type &__y) { return __x.__lambda_ == __y.__lambda_; }

            friend
            bool operator!=(const param_type &__x, const param_type &__y) { return !(__x == __y); }
        };

    private:
        param_type __p_;

    public:
        // constructors and reset functions

        explicit exponential_distribution(result_type __lambda = 1)
                : __p_(param_type(__lambda)) {}


        explicit exponential_distribution(const param_type &__p) : __p_(__p) {}


        void reset() {}

        // generating functions
        template<class _URNG>

        result_type operator()(_URNG &__g) { return (*this)(__g, __p_); }

        template<class _URNG>
        result_type operator()(_URNG &__g, const param_type &__p);

        // property functions

        result_type lambda() const { return __p_.lambda(); }


        param_type param() const { return __p_; }


        void param(const param_type &__p) { __p_ = __p; }


        result_type min() const { return 0; }


        result_type max() const { return std::numeric_limits<result_type>::infinity(); }

        friend
        bool operator==(const exponential_distribution &__x,
                        const exponential_distribution &__y) { return __x.__p_ == __y.__p_; }

        friend
        bool operator!=(const exponential_distribution &__x,
                        const exponential_distribution &__y) { return !(__x == __y); }
    };

    template<class _RealType>
    template<class _URNG>
    _RealType
    exponential_distribution<_RealType>::operator()(_URNG &__g, const param_type &__p) {
        return -std::log
                (
                        result_type(1) -
                        __generate_canonical<result_type,
                                std::numeric_limits<result_type>::digits>(__g)
                )
               / __p.lambda();
    }


// poisson distribution

    template<class _IntType>
    template<class _URNG>
    _IntType
    poisson_distribution<_IntType>::operator()(_URNG &__urng, const param_type &__pr) {
        result_type __x;
        uniform_real_distribution<double> __urd;
        if (__pr.__mean_ < 10) {
            __x = 0;
            for (double __p = __urd(__urng); __p > __pr.__l_; ++__x)
                __p *= __urd(__urng);
        } else {
            double __difmuk;
            double __g = __pr.__mean_ + __pr.__s_ * normal_distribution<double>()(__urng);
            double __u;
            if (__g > 0) {
                __x = static_cast<result_type>(__g);
                if (__x >= __pr.__l_)
                    return __x;
                __difmuk = __pr.__mean_ - __x;
                __u = __urd(__urng);
                if (__pr.__d_ * __u >= __difmuk * __difmuk * __difmuk)
                    return __x;
            }
            exponential_distribution<double> __edist;
            for (bool __using_exp_dist = false; true; __using_exp_dist = true) {
                double __e;
                if (__using_exp_dist || __g < 0) {
                    double __t;
                    do {
                        __e = __edist(__urng);
                        __u = __urd(__urng);
                        __u += __u - 1;
                        __t = 1.8 + (__u < 0 ? -__e : __e);
                    } while (__t <= -.6744);
                    __x = __pr.__mean_ + __pr.__s_ * __t;
                    __difmuk = __pr.__mean_ - __x;
                    __using_exp_dist = true;
                }
                double __px;
                double __py;
                if (__x < 10) {
                    const double __fac[] = {1, 1, 2, 6, 24, 120, 720, 5040,
                                            40320, 362880};
                    __px = -__pr.__mean_;
                    __py = std::pow(__pr.__mean_, (double) __x) / __fac[__x];
                } else {
                    double __del = .8333333E-1 / __x;
                    __del -= 4.8 * __del * __del * __del;
                    double __v = __difmuk / __x;
                    if (std::abs(__v) > 0.25)
                        __px = __x * std::log(1 + __v) - __difmuk - __del;
                    else
                        __px = __x * __v * __v * (((((((.1250060 * __v + -.1384794) *
                                                       __v + .1421878) * __v + -.1661269) * __v + .2000118) *
                                                    __v + -.2500068) * __v + .3333333) * __v + -.5) - __del;
                    __py = .3989423 / std::sqrt(__x);
                }
                double __r = (0.5 - __difmuk) / __pr.__s_;
                double __r2 = __r * __r;
                double __fx = -0.5 * __r2;
                double __fy = __pr.__omega_ * (((__pr.__c3_ * __r2 + __pr.__c2_) *
                                                __r2 + __pr.__c1_) * __r2 + __pr.__c0_);
                if (__using_exp_dist) {
                    if (__pr.__c_ * std::abs(__u) <= __py * std::exp(__px + __e) -
                                                     __fy * std::exp(__fx + __e))
                        break;
                } else {
                    if (__fy - __u * __fy <= __py * std::exp(__px - __fx))
                        break;
                }
            }
        }
        return __x;
    }


// discrete_distribution

    template<class _IntType = int>
    class  discrete_distribution {
    public:
        // types
        typedef _IntType result_type;

        class  param_type {
            std::vector<double> __p_;
        public:
            typedef discrete_distribution distribution_type;


            param_type() {}

            template<class _InputIterator>

            param_type(_InputIterator __f, _InputIterator __l)
                    : __p_(__f, __l) { __init(); }

#ifndef _LIBCPP_CXX03_LANG


            param_type(std::initializer_list<double> __wl)
                    : __p_(__wl.begin(), __wl.end()) { __init(); }

#endif  // _LIBCPP_CXX03_LANG

            template<class _UnaryOperation>
            param_type(size_t __nw, double __xmin, double __xmax,
                       _UnaryOperation __fw);

            std::vector<double> probabilities() const;

            friend
            bool operator==(const param_type &__x, const param_type &__y) { return __x.__p_ == __y.__p_; }

            friend
            bool operator!=(const param_type &__x, const param_type &__y) { return !(__x == __y); }

        private:
            void __init();

            friend class discrete_distribution;

        };

    private:
        param_type __p_;

    public:
        // constructor and reset functions

        discrete_distribution() {}

        template<class _InputIterator>

        discrete_distribution(_InputIterator __f, _InputIterator __l)
                : __p_(__f, __l) {}

#ifndef _LIBCPP_CXX03_LANG


        discrete_distribution(std::initializer_list<double> __wl)
                : __p_(__wl) {}

#endif  // _LIBCPP_CXX03_LANG

        template<class _UnaryOperation>

        discrete_distribution(size_t __nw, double __xmin, double __xmax,
                              _UnaryOperation __fw)
                : __p_(__nw, __xmin, __xmax, __fw) {}


        explicit discrete_distribution(const param_type &__p)
                : __p_(__p) {}


        void reset() {}

        // generating functions
        template<class _URNG>

        result_type operator()(_URNG &__g) { return (*this)(__g, __p_); }

        template<class _URNG>
        result_type operator()(_URNG &__g, const param_type &__p);

        // property functions

        std::vector<double> probabilities() const { return __p_.probabilities(); }


        param_type param() const { return __p_; }


        void param(const param_type &__p) { __p_ = __p; }


        result_type min() const { return 0; }


        result_type max() const { return __p_.__p_.size(); }

        friend
        bool operator==(const discrete_distribution &__x,
                        const discrete_distribution &__y) { return __x.__p_ == __y.__p_; }

        friend
        bool operator!=(const discrete_distribution &__x,
                        const discrete_distribution &__y) { return !(__x == __y); }

    };

    template<class _IntType>
    template<class _UnaryOperation>
    discrete_distribution<_IntType>::param_type::param_type(size_t __nw,
                                                            double __xmin,
                                                            double __xmax,
                                                            _UnaryOperation __fw) {
        if (__nw > 1) {
            __p_.reserve(__nw - 1);
            double __d = (__xmax - __xmin) / __nw;
            double __d2 = __d / 2;
            for (size_t __k = 0; __k < __nw; ++__k)
                __p_.push_back(__fw(__xmin + __k * __d + __d2));
            __init();
        }
    }

    template<class _IntType>
    void
    discrete_distribution<_IntType>::param_type::__init() {
        if (!__p_.empty()) {
            if (__p_.size() > 1) {
                double __s = std::accumulate(__p_.begin(), __p_.end(), 0.0);
                for (std::vector<double>::iterator __i = __p_.begin(), __e = __p_.end();
                     __i < __e; ++__i)
                    *__i /= __s;
                std::vector<double> __t(__p_.size() - 1);
                std::partial_sum(__p_.begin(), __p_.end() - 1, __t.begin());
                swap(__p_, __t);
            } else {
                __p_.clear();
                __p_.shrink_to_fit();
            }
        }
    }

    template<class _IntType>
    std::vector<double>
    discrete_distribution<_IntType>::param_type::probabilities() const {
        size_t __n = __p_.size();
        std::vector<double> __p(__n + 1);
        std::adjacent_difference(__p_.begin(), __p_.end(), __p.begin());
        if (__n > 0)
            __p[__n] = 1 - __p_[__n - 1];
        else
            __p[0] = 1;
        return __p;
    }

    template<class _IntType>
    template<class _URNG>
    _IntType
    discrete_distribution<_IntType>::operator()(_URNG &__g, const param_type &__p) {
        uniform_real_distribution<double> __gen;
        return static_cast<_IntType>(
                std::upper_bound(__p.__p_.begin(), __p.__p_.end(), __gen(__g)) -
                __p.__p_.begin());
    }


} // namespace random
# endif
