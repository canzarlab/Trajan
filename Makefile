CXX = g++
CFLAGS = -O2 -std=c++17 -pthread -march=native
INCL = -I.
VPATH = src
BINARIES = trajan bgen convert fpt
GENO_OBJS = main.o geno/augmentedLagrangian.o geno/lbfgsb.o geno/lineSearch.o

build: $(BINARIES) 

clean:
	rm -f $(BINARIES) $(TEST_BINARIES)
	rm -f *.o
	rm -f geno/*.o
	rm -f *~
	rm -f core

geno/%.o: geno/%.cpp
	$(CXX) -c $< -o $@ $(CFLAGS) $(INCL)

%.o: %.cpp $(HEADERS)
	$(CXX) -c $< $(CFLAGS) $(INCL)

trajan: $(GENO_OBJS) Graph.o Greedy.o Solver.o LP.o AntichainConstraint.o Constraint.o IndependentSetConstraint.o CrossingConstraint.o newick.o BnB.o BnG.o Similarity.o LPInt.o LPCP.o LPFInt.o Parallel.o EditDist.o read_csv.o
	$(CXX) -o $@ $^ $(CFLAGS)

fpt: match.o read_csv.o
	$(CXX) -o $@ $^ $(CFLAGS)

bgen: generator.o newick.o
	$(CXX) -o $@ $^ $(CFLAGS)

convert: convert.o
	$(CXX) -o $@ $^ $(CFLAGS) -lboost_graph
