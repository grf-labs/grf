#ifndef GRF_GLOBALS_H_
#define GRF_GLOBALS_H_

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
    TypeName(const TypeName&);             \
    void operator=(const TypeName&)

namespace grf {

typedef unsigned int uint;

static const uint DEFAULT_NUM_THREADS = 0;

} // namespace grf
#endif /* GRF_GLOBALS_H_ */
