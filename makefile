run-tests:
	cd cmake-build && ctest --extra-verbose --output-on-failure

run:
	echo "Copying $(SET).txt from SETS" && cp -r SETS/$(SET).txt ./cmake-exec/$(SET)/ && cd cmake-exec/$(SET) && ./$(SET) $(TXT)

run-vis:
	cd visualtests-out && ./visualtests

build-dev:
	cmake -B ./cmake-build -DOPT=$(OPT) -DOMP=$(OMP) -DVIS=$(VIS) -DSET=$(SET) -DTXT=$(TXT) -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ && cmake --build ./cmake-build

build:
	cmake -DOPT=ON -DOMP=ON -B ./cmake-build -DSET=$(SET) -DTXT=$(TXT) && cmake --build ./cmake-build

clean:
	rm -rf *build*/ && rm -rf *exec*/

deps:
	rm -rf deps && git clone https://github.com/kulakovri/MDOODZ-dependencies deps && cd deps && make install-hdf5 && make install-suitesparse

install-suitesparse:
	rm -rf deps && git clone https://github.com/kulakovri/MDOODZ-dependencies deps && cd deps && make install-suitesparse