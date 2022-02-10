run-tests:
	cd cmake-build && ./ShearTemplate_test

run:
	cd cmake-build && ./MDOODZ

build:
	cmake -G "Unix Makefiles" -DMODEL=$(MODEL) -B ./cmake-build && cmake --build ./cmake-build -- -j 6

clean:
	rm -rf *build*/

install-deps:
	rm -rf lib && git clone https://github.com/kulakovri/MDOODZ-dependencies lib && cd lib && make install-hdf5 && make install-suitesparse