#!/bin/bash

PURPLE='\033[1;35m'
NC='\033[0m'

echo -e "${PURPLE}Testing quantile forests...${NC}"
./build/ranger --verbose --file ../smoke_test/quantile.dat --depvarname Y --treetype 11 --ntree 4 --nthreads 4 --write
./build/ranger --verbose --file ../smoke_test/quantile.dat --depvarname Y --treetype 11 --ntree 4 --nthreads 4 --predict ranger.forest

echo ""
echo -e "${PURPLE}Testing instrumental forests...${NC}"
./build/ranger --verbose --file ../smoke_test/causal.dat --depvarname Y --statusvarname W --instrumentvarname W --treetype 15 --ntree 4 --nthreads 4 --write
./build/ranger --verbose --file ../smoke_test/causal.dat --depvarname Y --statusvarname W --instrumentvarname W --treetype 15 --ntree 4 --nthreads 4 --predict ranger.forest