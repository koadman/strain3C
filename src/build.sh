#!/bin/bash
#
# compiles the stan program and applies the patch to the C++ to implement
# stochastic gradient descent
#
STAN_PREFIX=/home/koadman/software/cmdstan-2.21.0


$STAN_PREFIX/bin/stanc  --o=src/genotypes3C.hpp src/genotypes3C.stan && \
src/stochastic_patch.py < src/genotypes3C.hpp > src/genotypes3C_stochastic.hpp && \
g++ -std=c++1y -pthread -D_REENTRANT -Wno-sign-compare     -I $STAN_PREFIX/stan/lib/stan_math/lib/tbb_2019_U8/include -O3 -I $STAN_PREFIX/src -I $STAN_PREFIX/stan/src -I $STAN_PREFIX/lib/rapidjson_1.1.0/ -I $STAN_PREFIX/stan/lib/stan_math/ -I $STAN_PREFIX/stan/lib/stan_math/lib/eigen_3.3.3 -I $STAN_PREFIX/stan/lib/stan_math/lib/boost_1.69.0 -I $STAN_PREFIX/stan/lib/stan_math/lib/sundials_4.1.0/include    -DBOOST_DISABLE_ASSERTS      -c  -x c++ -o src/genotypes3C.o src/genotypes3C_stochastic.hpp && \
g++ -std=c++1y -pthread -D_REENTRANT -Wno-sign-compare     -I $STAN_PREFIX/stan/lib/stan_math/lib/tbb_2019_U8/include -O3 -I $STAN_PREFIX/src -I $STAN_PREFIX/stan/src -I $STAN_PREFIX/lib/rapidjson_1.1.0/ -I $STAN_PREFIX/stan/lib/stan_math/ -I $STAN_PREFIX/stan/lib/stan_math/lib/eigen_3.3.3 -I $STAN_PREFIX/stan/lib/stan_math/lib/boost_1.69.0 -I $STAN_PREFIX/stan/lib/stan_math/lib/sundials_4.1.0/include    -DBOOST_DISABLE_ASSERTS            -Wl,-L,"$STAN_PREFIX/stan/lib/stan_math/lib/tbb" -Wl,-rpath,"$STAN_PREFIX/stan/lib/stan_math/lib/tbb"  src/genotypes3C.o $STAN_PREFIX/src/cmdstan/main.o         $STAN_PREFIX/stan/lib/stan_math/lib/sundials_4.1.0/lib/libsundials_nvecserial.a $STAN_PREFIX/stan/lib/stan_math/lib/sundials_4.1.0/lib/libsundials_cvodes.a $STAN_PREFIX/stan/lib/stan_math/lib/sundials_4.1.0/lib/libsundials_idas.a $STAN_PREFIX/stan/lib/stan_math/lib/sundials_4.1.0/lib/libsundials_kinsol.a  $STAN_PREFIX/stan/lib/stan_math/lib/tbb/libtbb.so.2 -o bin/genotypes3C


$STAN_PREFIX/bin/stanc  --o=src/hierarchical_genotypes.hpp src/hierarchical_genotypes.stan && \
src/stochastic_patch.py < src/hierarchical_genotypes.hpp > src/hierarchical_genotypes_stochastic.hpp && \
g++ -std=c++1y -pthread -D_REENTRANT -Wno-sign-compare     -I $STAN_PREFIX/stan/lib/stan_math/lib/tbb_2019_U8/include -O3 -I $STAN_PREFIX/src -I $STAN_PREFIX/stan/src -I $STAN_PREFIX/lib/rapidjson_1.1.0/ -I $STAN_PREFIX/stan/lib/stan_math/ -I $STAN_PREFIX/stan/lib/stan_math/lib/eigen_3.3.3 -I $STAN_PREFIX/stan/lib/stan_math/lib/boost_1.69.0 -I $STAN_PREFIX/stan/lib/stan_math/lib/sundials_4.1.0/include    -DBOOST_DISABLE_ASSERTS      -c  -x c++ -o src/hierarchical_genotypes.o src/hierarchical_genotypes_stochastic.hpp && \
g++ -std=c++1y -pthread -D_REENTRANT -Wno-sign-compare     -I $STAN_PREFIX/stan/lib/stan_math/lib/tbb_2019_U8/include -O3 -I $STAN_PREFIX/src -I $STAN_PREFIX/stan/src -I $STAN_PREFIX/lib/rapidjson_1.1.0/ -I $STAN_PREFIX/stan/lib/stan_math/ -I $STAN_PREFIX/stan/lib/stan_math/lib/eigen_3.3.3 -I $STAN_PREFIX/stan/lib/stan_math/lib/boost_1.69.0 -I $STAN_PREFIX/stan/lib/stan_math/lib/sundials_4.1.0/include    -DBOOST_DISABLE_ASSERTS            -Wl,-L,"$STAN_PREFIX/stan/lib/stan_math/lib/tbb" -Wl,-rpath,"$STAN_PREFIX/stan/lib/stan_math/lib/tbb"  src/hierarchical_genotypes.o $STAN_PREFIX/src/cmdstan/main.o         $STAN_PREFIX/stan/lib/stan_math/lib/sundials_4.1.0/lib/libsundials_nvecserial.a $STAN_PREFIX/stan/lib/stan_math/lib/sundials_4.1.0/lib/libsundials_cvodes.a $STAN_PREFIX/stan/lib/stan_math/lib/sundials_4.1.0/lib/libsundials_idas.a $STAN_PREFIX/stan/lib/stan_math/lib/sundials_4.1.0/lib/libsundials_kinsol.a  $STAN_PREFIX/stan/lib/stan_math/lib/tbb/libtbb.so.2 -o bin/hierarchical_genotypes
