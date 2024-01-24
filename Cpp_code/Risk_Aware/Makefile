CC = g++
CXXFLAGS = -std=c++17 -O3 -fopenmp
BOOSTPATH = ../boost_1_79_0
INCLUDE = -I $(BOOSTPATH) -I ./include/
LIB = -L $(BOOSTPATH)
SRC = ./src/
TEST =

CPPFILES = $(SRC)main.cpp \
					$(SRC)CGrid_2D.cpp \
					$(SRC)CENO.cpp \
					$(SRC)CMarch_Switch.cpp \
					$(SRC)CSpecial_Treatment.cpp \
					$(SRC)Sailing_Threshold_Solver.cpp \

main: $(CPPFILES)
	$(CC) $(CXXFLAGS) -o $@ $^ $(INCLUDE) $(LIB)

run: main
		./main $(TEST)

debug: $(CPPFILES)
	g++ -o1 -g -std=c++17 $(INCLUDE) $^ $(LIB)
	gdb -ex=r --args 1 $(TEST)

clean:
	rm *.dat main
