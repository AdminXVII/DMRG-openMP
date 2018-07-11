KOKKOS_SETTINGS=KOKKOS_SRC_PATH=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos CXX=g++ KOKKOS_DEVICES=OpenMP KOKKOS_ARCH=Power8 GTEST_PATH=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos/tpls/gtest KOKKOS_OPTIONS=compiler_warnings,disable_profiling KOKKOS_PATH=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos

kokkoslib:
	cd core; \
	$(MAKE) -f /lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos/core/src/Makefile KOKKOS_SRC_PATH=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos CXX=g++ KOKKOS_DEVICES=OpenMP KOKKOS_ARCH=Power8 GTEST_PATH=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos/tpls/gtest KOKKOS_OPTIONS=compiler_warnings KOKKOS_PATH=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos PREFIX=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos/build/gcc_6.3/openMP/install build-lib

install: kokkoslib
	cd core; \
	$(MAKE) -f /lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos/core/src/Makefile KOKKOS_SRC_PATH=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos CXX=g++ KOKKOS_DEVICES=OpenMP KOKKOS_ARCH=Power8 GTEST_PATH=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos/tpls/gtest KOKKOS_OPTIONS=compiler_warnings KOKKOS_PATH=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos PREFIX=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos/build/gcc_6.3/openMP/install install

kokkoslib-test:
	cd core; \
	$(MAKE) -f /lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos/core/src/Makefile KOKKOS_SRC_PATH=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos CXX=g++ KOKKOS_DEVICES=OpenMP KOKKOS_ARCH=Power8 GTEST_PATH=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos/tpls/gtest KOKKOS_OPTIONS=compiler_warnings KOKKOS_PATH=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos PREFIX=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos/build/gcc_6.3/openMP/install build-lib

install-test: kokkoslib-test
	cd core; \
	$(MAKE) -f /lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos/core/src/Makefile KOKKOS_SRC_PATH=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos CXX=g++ KOKKOS_DEVICES=OpenMP KOKKOS_ARCH=Power8 GTEST_PATH=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos/tpls/gtest KOKKOS_OPTIONS=compiler_warnings KOKKOS_PATH=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos PREFIX=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos/build/gcc_6.3/openMP/install install

build-test: install-test
	$(MAKE) -C core/unit_test
	$(MAKE) -C core/perf_test
	$(MAKE) -C containers/unit_tests
	$(MAKE) -C containers/performance_tests
	$(MAKE) -C algorithms/unit_tests
	$(MAKE) -C example/fixture
	$(MAKE) -C example/feint
	$(MAKE) -C example/fenl
	$(MAKE) -C example/tutorial build

test: build-test
	$(MAKE) -C core/unit_test test
	$(MAKE) -C core/perf_test test
	$(MAKE) -C containers/unit_tests test
	$(MAKE) -C containers/performance_tests test
	$(MAKE) -C algorithms/unit_tests test
	$(MAKE) -C example/fixture test
	$(MAKE) -C example/feint test
	$(MAKE) -C example/fenl test
	$(MAKE) -C example/tutorial test

unit-tests-only:
	$(MAKE) -C core/unit_test test
	$(MAKE) -C containers/unit_tests test
	$(MAKE) -C algorithms/unit_tests test

clean:
	$(MAKE) -C core/unit_test clean
	$(MAKE) -C core/perf_test clean
	$(MAKE) -C containers/unit_tests clean
	$(MAKE) -C containers/performance_tests clean
	$(MAKE) -C algorithms/unit_tests clean
	$(MAKE) -C example/fixture clean
	$(MAKE) -C example/feint clean
	$(MAKE) -C example/fenl clean
	$(MAKE) -C example/tutorial clean
	cd core; \
	$(MAKE) -f /lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos/core/src/Makefile KOKKOS_SRC_PATH=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos CXX=g++ KOKKOS_DEVICES=OpenMP KOKKOS_ARCH=Power8 GTEST_PATH=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos/tpls/gtest KOKKOS_OPTIONS=compiler_warnings KOKKOS_PATH=/lustre/atlas/scratch/arghyac/csc205/summit_miniapp/miniapp_ronnie/kokkos/kokkos clean
