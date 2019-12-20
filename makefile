#SYSTEM = ParabolaLine
#SYSTEM = TwoEllipses
#SYSTEM = WormLoop
#SYSTEM = Trimer
#SYSTEM = Plane
#SYSTEM = ParabolaLine0
#SYSTEM = Trimer0
#SYSTEM = Polymer
#SYSTEM = PolymerWall
SYSTEM = EllipsoidSurface
PROGRAM = test$(SYSTEM)
EQUATIONFILE = Equations$(SYSTEM)
PROPOSALFILE = Proposals


CC=g++ -std=gnu++11 
OPTFLAGS = -O3 -I /home/holmes/Programs/
EXECUTABLE=./$(PROGRAM) 

all: $(EXECUTABLE) run 


$(PROGRAM): $(PROGRAM).o Stratification.o Point.o SampleStrat.o Move.o Proposals.o \
		Equations.o
	$(CC) -o $(PROGRAM) $^

$(PROGRAM).o: $(PROGRAM).cpp Stratification.hpp Point.hpp SampleStrat.hpp Move.hpp  \
        $(EQUATIONFILE).cpp Proposals.hpp \
		makefile
	$(CC) -c $(OPTFLAGS) $(PROGRAM).cpp 

Stratification.o: Stratification.cpp Stratification.hpp Point.hpp 
	$(CC) -c $(OPTFLAGS) Stratification.cpp

Equations.o: $(EQUATIONFILE).cpp Stratification.hpp SampleStrat.hpp Move.hpp Point.hpp \
             makefile
	$(CC) -c $(OPTFLAGS) $(EQUATIONFILE).cpp -o Equations.o

Point.o: Point.cpp Point.hpp Stratification.hpp  
	$(CC) -c $(OPTFLAGS) Point.cpp

SampleStrat.o: SampleStrat.cpp Point.hpp Move.hpp Stratification.hpp Proposals.hpp \
		makefile
	$(CC) -c $(OPTFLAGS) SampleStrat.cpp 

Move.o: Move.cpp Point.hpp Stratification.hpp
	$(CC) -c $(OPTFLAGS) Move.cpp 

Proposals.o: $(PROPOSALFILE).cpp Proposals.hpp Point.hpp Move.hpp Stratification.hpp  makefile
	$(CC) -c $(OPTFLAGS) $(PROPOSALFILE).cpp -o Proposals.o


run: $(EXECUTABLE)
	$(EXECUTABLE) $(ARGS) 

clean:
	rm -f *.o 


