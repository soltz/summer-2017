cxxflags=`fastjet-config --cxxflags`
# compile
obj_files=`echo `
g++ -std=c++11 -c -o .obj/Utilities.o Utilities.cpp $cxxflags
g++ -std=c++11 -c -o .obj/Event.o Event.cpp $cxxflags
g++ -std=c++11 -c -o .obj/Background.o Background.cpp $cxxflags
g++ -std=c++11 -c -o .obj/GenerateDataOriginal.o GenerateDataOriginal.cpp $cxxflags -I.
g++ -std=c++11 -c -o .obj/GenerateDataWithBackground.o GenerateDataWithBackground.cpp $cxxflags -I.
g++ -std=c++11 -c -o .obj/GenerateDataWithoutBackground.o GenerateDataWithoutBackground.cpp $cxxflags -I.
#g++ -std=c++11 -c -o .obj/GenerateDataWithBackgroundWithBackgroundSubtraction.o GenerateDataWithBackgroundWithBackgroundSubtraction.cpp $cxxflags -I.
# link
g++ -std=c++11 -o GenerateDataOriginal .obj/GenerateDataOriginal.o .obj/Utilities.o .obj/Event.o .obj/Background.o $cxxflags `fastjet-config --libs` -lHepMC -lpythia8
g++ -std=c++11 -o GenerateDataWithBackground .obj/GenerateDataWithBackground.o .obj/Utilities.o .obj/Event.o .obj/Background.o $obj_files $cxxflags `fastjet-config --libs` -lHepMC -lpythia8
g++ -std=c++11 -o GenerateDataWithoutBackground .obj/GenerateDataWithoutBackground.o .obj/Utilities.o .obj/Event.o .obj/Background.o $obj_files $cxxflags `fastjet-config --libs` -lHepMC -lpythia8
#g++ -std=c++11 -o GenerateDataWithBackgroundWithBackgroundSubtraction .obj/GenerateDataWithBackgroundWithBackgroundSubtraction.o .obj/Utilities.o .obj/Event.o .obj/Background.o $obj_files $cxxflags `fastjet-config --libs` -lHepMC -lpythia8
