CXXFLAGS= -std=gnu++0x -Wall -pedantic 
LIBSFLAG = -02
LIBS=-lgmpxx -lgmp -lmpfr

all: factor

factor: main.cpp qs.o
	$(CXX) -c $(CXXFLAGS) -o main.o main.cpp
	$(CXX) $(LIBSFLAGS) -o factor main.o qs.o $(LIBS)

qs.o: qs.cpp qs.h mathHelper.h 
	$(CXX) -c $(CXXFLAGS) -o qs.o qs.cpp
