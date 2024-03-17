# Makefile
CXX = g++
CXXFLAGS = -O2

objects = main.o Energy.o Energyquick.o pressure.o unewquick.o uquick.o uvnew.o uvupwind.o vnewquick.o vquick.o vupwind.o vvnew.o

diffHeatedCavity: $(objects)
	$(CXX) -o diffHeatedCavity $(objects)

main.o: main.cpp uvupwind.h uquick.h uvnew.h vquick.h vupwind.h vvnew.h pressure.h unewquick.h vnewquick.h Energy.h Energyquick.h
	$(CXX) $(CXXFLAGS) -c main.cpp

uvupwind.o: uvupwind.cpp 
	$(CXX) $(CXXFLAGS) -c uvupwind.cpp

vupwind.o: vupwind.cpp
	$(CXX) $(CXXFLAGS) -c vupwind.cpp

uquick.o: uquick.cpp 
	$(CXX) $(CXXFLAGS) -c uquick.cpp

vquick.o: vquick.cpp
	$(CXX) $(CXXFLAGS) -c vquick.cpp

pressure.o: pressure.cpp 
	$(CXX) $(CXXFLAGS) -c pressure.cpp

uvnew.o: uvnew.cpp 
	$(CXX) $(CXXFLAGS) -c uvnew.cpp

vvnew.o: vvnew.cpp
	$(CXX) $(CXXFLAGS) -c vvnew.cpp

unewquick.o: unewquick.cpp
	$(CXX) $(CXXFLAGS) -c unewquick.cpp

vnewquick.o: vnewquick.cpp
	$(CXX) $(CXXFLAGS) -c vnewquick.cpp

Energy.o: Energy.cpp
	$(CXX) $(CXXFLAGS) -c Energy.cpp

Energyquick.o: Energyquick.cpp
	$(CXX) $(CXXFLAGS) -c Energyquick.cpp

clean:
	rm diffHeatedCavity $(objects)




