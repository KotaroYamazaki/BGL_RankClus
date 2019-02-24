#program name and objects file name
PROGRAM = rankclus.out
OBJS = main.o graph.o clustering.o
OPTIMIZE = -O2
#define macro
CXX = clang++
CXXFLAG = $(OPTIMIZE) -std=c++11

all:authority
authority: $(OBJS) ranking.o
	$(CXX) $(CXXFLAG) $(OBJS) ranking.o -o $(PROGRAM)
pagerank: $(OBJS) single_pagerank.o
	$(CXX) $(CXXFLAG) $(OBJS) single_pagerank.o -o $(PROGRAM)
.cpp.o:
	$(CXX) $(CXXFLAG) -c $< -Wall
*.o: graph.hpp
.PHONY: clean
clean:
	$(RM) *.o *test*
