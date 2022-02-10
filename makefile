run-tests:
	cd cmake-build && ./ShearTemplate_test

run:
	cd cmake-build && ./MDOODZ

build:
	cmake -G "Unix Makefiles" -DMODEL=$(MODEL) -B ./cmake-build && cmake --build ./cmake-build -- -j 6

clean:
	rm -rf *build*/

install-hdf5:
	git clone https://github.com/HDFGroup/hdf5 lib/hdf5 && cd lib/hdf5 && mkdir install && cmake -S . -B _build -DCMAKE_INSTALL_PREFIX=install && cd _build && cmake --build . -j --target install

install-suitesparse:
	git clone https://github.com/jlblancoc/suitesparse-metis-for-windows/ lib/suitesparse && cd lib/suitesparse && mkdir install && cmake BUILD_METIS=0 -S . -B _build -DCMAKE_INSTALL_PREFIX=install && cd _build && cmake --build . -j --target install