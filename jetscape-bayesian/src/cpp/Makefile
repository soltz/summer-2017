#
#       Makefile for summer-2017/jetscape-bayesian/src/cpp programs
#
CXX=g++
CXXFLAGS=-std=c++11 -c `fastjet-config --cxxflags` 
LDFLAGS=--std=c++11 -lHEpMC -lpythia8 `fastjet-config --libs` 
obj_files=.obj/Utilities.o .obj/Event.o .obj/Background.o

all: GenerateDataOriginal GenerateDataWithBackground GenerateDataWithoutBackground

GenerateDataOriginal: .obj/GenerateDataOriginal.o $(obj_files) 
	$(CXX) $^ -o $@ $(LDFLAGS) 

GenerateDataWithBackground: .obj/GenerateDataWithBackground.o $(obj_files)
	$(CXX) $^ -o $@ $(LDFLAGS)

GenerateDataWithoutBackground: .obj/GenerateDataWithoutBackground.o $(obj_files)
	$(CXX) $^ -o $@ $(LDFLAGS)

# Add includes
.obj/Utilities.o : Utilities.cpp Utilities.h

.obj/Event.o : Event.cpp Event.h

.obj/Background.o : Background.cpp Background.h

# Use pattern rule to generate .obj/object files
.obj/%.o : %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@
