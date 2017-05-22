#ifndef GRF_GLOBALS_H_
#define GRF_GLOBALS_H_

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
    TypeName(const TypeName&);             \
    void operator=(const TypeName&)

// Pi
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef unsigned int uint;

enum TreeType {
  TREE_QUANTILE = 11,
  TREE_INSTRUMENTAL = 15,
};

// Mask and Offset to store 2 bit values in bytes
static const int mask[4] = {192,48,12,3};
static const int offset[4] = {6,4,2,0};


// Default values
const uint DEFAULT_NUM_TREE = 500;
const uint DEFAULT_NUM_THREADS = 0;

const uint DEFAULT_MIN_NODE_SIZE_CLASSIFICATION = 1;
const uint DEFAULT_MIN_NODE_SIZE_REGRESSION = 5;
const uint DEFAULT_MIN_NODE_SIZE_PROBABILITY = 10;

// Interval to print progress in seconds
const double STATUS_INTERVAL = 30.0;

// Threshold for q value split method switch
const double Q_THRESHOLD = 0.02;

#endif /* GRF_GLOBALS_H_ */
