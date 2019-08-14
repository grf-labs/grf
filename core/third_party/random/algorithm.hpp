#include <limits>
#include <cmath>
#include <algorithm>
#include <iterator>
#include "random/random.hpp"

#ifndef _GRFSTD_ALGORITHM
#define _GRFSTD_ALGORITHM

namespace nonstd {

        template<class _RandomAccessIterator, class _UniformRandomNumberGenerator>
        void shuffle(_RandomAccessIterator __first, _RandomAccessIterator __last,
#ifndef _LIBCPP_CXX03_LANG
                     _UniformRandomNumberGenerator &&__g)
#else
        _UniformRandomNumberGenerator& __g)
#endif
        {
            typedef typename std::iterator_traits<_RandomAccessIterator>::difference_type difference_type;
            typedef nonstd::uniform_int_distribution<ptrdiff_t> _Dp;
            typedef typename _Dp::param_type _Pp;
            difference_type __d = __last - __first;
            if (__d > 1) {
                _Dp __uid;
                for (--__last, --__d; __first < __last; ++__first, --__d) {
                    difference_type __i = __uid(__g, _Pp(0, __d));
                    if (__i != difference_type(0))
                        std::swap(*__first, *(__first + __i));
                }
            }
        }
}

#endif  // _GRFSTD_ALGORITHM
