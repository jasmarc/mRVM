CC = g++
CFLAGS = -g -Wall -L/opt/local/lib -I/opt/local/include -I${GTEST_DIR}/include -I${GTEST_DIR}
GSLFLAGS = -lgsl -lgslcblas -lm -I/opt/local/include/gsl
EXEC = mRVM
GTEST_DIR = ./gtest-1.5.0

alltest: cleantest test lint reseed runtest

test:
	$(CC) $(CFLAGS) -c ${GTEST_DIR}/src/gtest-all.cc
	ar -rv libgtest.a gtest-all.o
	$(CC) $(CFLAGS) $(GSLFLAGS) ${GTEST_DIR}/src/gtest_main.cc test.cc libgtest.a -o test

gsl_poly_orth:
	$(CC) -L/opt/local/lib -I/opt/local/include -I${GTEST_DIR}/include -I${GTEST_DIR} $(GSLFLAGS) gsl_poly_orth.c -o gsl_poly_orth.o

lint:
	python cpplint.py test.cc

reseed:
	GSL_RNG_SEED=`date +%s`
	echo $(GSL_RNG_SEED)

cleantest:
	-rm test
	-rm *.exe

runtest:
	./test

$(EXEC): 
	$(CC) $(CFLAGS) -o $(EXEC) main.cc

clean:
	-rm $(EXEC)
	-rm *.exe