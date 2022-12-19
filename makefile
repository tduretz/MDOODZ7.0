run-tests:
	cd cmake-build && ctest --extra-verbose --output-on-failure

run:
	cd cmake-exec/$(SET) && ./$(SET)

run-vis:
	cd visualtests-out && ./visualtests

build-dev:
	cmake -G "MinGW Makefiles" -B ./cmake-build -DOPT=$(OPT) -DOMP=$(OMP) -DVIS=$(VIS) -DSET=$(SET) -DCMAKE_C_COMPILER=c:/MinGW/bin/gcc -DCMAKE_CXX_COMPILER=g++ && cmake --build ./cmake-build

build:
	cmake -DOPT=ON -DOMP=ON -B ./cmake-build -DSET=$(SET) && cmake --build ./cmake-build

clean:
	rm -rf *build*/ && rm -rf *exec*/

deps:
	rm -rf deps && git clone https://github.com/kulakovri/MDOODZ-dependencies deps && cd deps && make install-hdf5 && make install-suitesparse

install-suitesparse:
	rm -rf deps && git clone https://github.com/kulakovri/MDOODZ-dependencies deps && cd deps && make install-suitesparse