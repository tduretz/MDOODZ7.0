git clone https://github.com/jlblancoc/suitesparse-metis-for-windows lib/SuiteSparse
cmake -S lib/SuiteSparse -B lib/SuiteSparse/_build
cd lib/SuiteSparse/_build
make amd umfpack cholmod ccolamd btf lapack cxsparse cblas