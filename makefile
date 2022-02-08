run-tests:
	cd cmake-build && ./ShearTemplate_test

run:
	cd cmake-build && ./MDOODZ

build:
	cmake -G "Unix Makefiles" -DMODEL=$(MODEL) -B ./cmake-build && cmake --build ./cmake-build -- -j 6

clean:
	rm -rf *build*/