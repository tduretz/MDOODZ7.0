run-tests:
	cd cmake-build && ctest --extra-verbose

run:
	cd cmake-exec/$(SET) && ./$(SET)

build-dev:
	cmake -B ./cmake-build -DOPT=$(OPT) -DOMP=$(OMP) && cmake --build ./cmake-build

build:
	cmake -DOPT=ON -DOMP=ON -B ./cmake-build && cmake --build ./cmake-build

clean:
	rm -rf *build*/

deps:
	rm -rf lib && git clone https://github.com/kulakovri/MDOODZ-dependencies deps && cd deps && make install-hdf5 && make install-suitesparse