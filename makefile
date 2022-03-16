run-tests:
	cd cmake-build && ./shearPwl_test && ./ShearTemplate_test

run:
	cd cmake-build && ./MDOODZ

build-dev:
	cmake -G "Unix Makefiles" -DMODEL=$(MODEL) -B ./cmake-build && cmake --build ./cmake-build -- -j 6

build:
	cmake -G "Unix Makefiles" -DMODEL=$(MODEL) -DOPT=ON -DOMP=ON -B ./cmake-build && cmake --build ./cmake-build -- -j 6

clean:
	rm -rf *build*/

deps:
	rm -rf lib && git clone https://github.com/kulakovri/MDOODZ-dependencies lib && cd lib && make install-hdf5 && make install-suitesparse