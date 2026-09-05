all: help

help:
	@echo "make install   -- install dependencies and PyMultiNest"
	@echo "make test      -- run demo tests"

install:
	git clone https://github.com/JohannesBuchner/MultiNest || true
	mkdir -p MultiNest/build
	cd MultiNest/build && cmake .. && make
	test -e MultiNest/lib/libmultinest.so
	python -m pip install .

test:
	LD_LIBRARY_PATH=MultiNest/lib/:$(LD_LIBRARY_PATH) python pymultinest_demo_minimal.py
	LD_LIBRARY_PATH=MultiNest/lib/:$(LD_LIBRARY_PATH) python pymultinest_demo.py
	LD_LIBRARY_PATH=MultiNest/lib/:$(LD_LIBRARY_PATH) python tests/pymultinest_long_name_test.py
	LD_LIBRARY_PATH=MultiNest/lib/:$(LD_LIBRARY_PATH) python multinest_marginals.py chains/3-
	LD_LIBRARY_PATH=MultiNest/lib/:$(LD_LIBRARY_PATH) python multinest_marginals_corner.py chains/3-

.PHONY: all help install test
