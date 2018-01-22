#ifndef GRF_GLOBALS_H_
#define GRF_GLOBALS_H_

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
    TypeName(const TypeName&);             \
    void operator=(const TypeName&)

typedef unsigned int uint;

static const uint DEFAULT_NUM_THREADS = 0;
static const double Q_THRESHOLD = 0.02;

#endif /* GRF_GLOBALS_H_ */
