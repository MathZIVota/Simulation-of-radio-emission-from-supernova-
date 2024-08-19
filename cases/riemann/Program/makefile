.PHONY: program clean

program: main.o AnalyticSolver.o loadTofile.o Tube.o
	g++ main.o AnalyticSolver.o loadTofile.o Tube.o -o program 

Tube.o: Tube.cpp Tube.h
	g++ -c Tube.cpp

loadTofile.o: loadTofile.cpp Tube.h
	g++ -c loadTofile.cpp

AnalyticSolver.o: AnalyticSolver.cpp AnalyticSolver.h
	g++ -c AnalyticSolver.cpp

main.o: main.cpp 
	g++ -c main.cpp

clean:
	rm *.o program