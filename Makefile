CC = g++
CFLAGS = -g -Wall -L/usr/local/lib -I/usr/local/include -I$(SRC_DIR) -I${GTEST_DIR}/include -I${GTEST_DIR}
GSLFLAGS = -L/usr/local/Cellar/gsl/1.14/lib -lgsl -lgslcblas -lm -I/usr/local/Cellar/gsl/1.14/include/gsl
EXEC = mRVM
GTEST_DIR = ./lib/gtest-1.5.0
SRC_DIR = ./src
TEST_DIR = ./src/test
OUTPUT_DIR = ./bin
TOOLS_DIR = ./tools

all: clean mRVM

alltest: clean test runtest

test:
	$(CC) $(CFLAGS) -c ${GTEST_DIR}/src/gtest-all.cc -o $(OUTPUT_DIR)/gtest-all.o
	ar -rv $(OUTPUT_DIR)/libgtest.a $(OUTPUT_DIR)/gtest-all.o
	$(CC) $(CFLAGS) $(GSLFLAGS) ${GTEST_DIR}/src/gtest_main.cc $(TEST_DIR)/test.cc $(OUTPUT_DIR)/libgtest.a -o $(OUTPUT_DIR)/test

lint:
	python $(TOOLS_DIR)/cpplint.py $(TEST_DIR)/test.cc

reseed:
	GSL_RNG_SEED=`date +%s`
	echo $(GSL_RNG_SEED)

clean:
	-rm $(OUTPUT_DIR)/test
	-rm $(OUTPUT_DIR)/*.exe

runtest:
	./$(OUTPUT_DIR)/test --gtest_filter=*update_y*

$(EXEC): 
	$(CC) $(CFLAGS) $(GSLFLAGS) -o $(OUTPUT_DIR)/$(EXEC) \
		$(SRC_DIR)/lib/Vector.cc \
		$(SRC_DIR)/lib/Matrix.cc \
		$(SRC_DIR)/lib/Trainer.cc \
		$(SRC_DIR)/lib/Predictor.cc \
		$(SRC_DIR)/lib/Kernel.cc \
		$(SRC_DIR)/lib/LinearKernel.cc \
		$(SRC_DIR)/lib/PolynomialKernel.cc \
		$(SRC_DIR)/lib/GaussianKernel.cc \
		$(SRC_DIR)/lib/RandomNumberGenerator.cc \
		$(SRC_DIR)/lib/GaussHermiteQuadrature.cc \
		$(SRC_DIR)/main.cc