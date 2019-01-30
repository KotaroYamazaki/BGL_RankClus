#program name and objects file name
PROGRAM = rankclus.out
OBJS = main.o ranking.o graph.o clustering.o
OPTIMIZE = -O2
#define macro
CXX = clang++
CXXFLAG = $(OPTIMIZE) -std=c++11

all:real
real: $(OBJS)
	$(CXX) $(CXXFLAG) $(OBJS) -o $(PROGRAM)
.cpp.o:
	$(CXX) $(CXXFLAG) -c $< -Wall
*.o: graph.hpp
.PHONY: clean
clean:
	$(RM) *.o
