INCLUDES := atomsParam.h  sphere_lebedev_rule.h rhoMod.h sharedkld.h
SRCLIBS :=  atomsParam.cpp  main.cpp rhoMod.cpp sharedkld.cpp sphere_lebedev_rule.cpp 
OBJECTS  := $(patsubst %.cpp, %.o, $(SRCLIBS))

#FOR LINUX
LIBS := -I/opt/intel/mkl/include/intel64 -I/opt/intel/mkl/lib/intel64 
LIBS_FRAMEWORK := /opt/intel/mkl/lib/intel64
CXX := icpc
CXXFLAGS := -std=c++14 -fopenmp -lmkl_rt -lpthread 
# END LINUX

#FOR MAC
#LIBS := -framework Accelerate -I/opt/intel/mkl/include -I/opt/intel/mkl/lib -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -I/Library/Developer/CommandLineTools/usr/include/c++/v1
#LIBS_FRAMEWORK := /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework
#CXX := icpc
#CXXFLAGS := -std=c++14 -fopenmp -lLAPACK
# END MAC

CFLAGS := -g
APPName := kld1d

all: clean main

main: atomsParam.o main.o rhoMod.o sharedkld.o sphere_lebedev_rule.o
	$(CXX) -o $(APPName) $(OBJECTS) -L$(LIBS_FRAMEWORK)  $(CXXFLAGS)
	$(chmod +x $(APPName))
atomsParam.o:
	$(CXX) -c $(CFLAGS) $(LIBS) $(CXXFLAGS) -MMD -MP -MF "atomsParam.o.d" -o atomsParam.o atomsParam.cpp
main.o:
	$(CXX) -c $(CFLAGS) $(LIBS) $(CXXFLAGS) -MMD -MP -MF "main.o.d" -o main.o main.cpp
rhoMod.o:
	$(CXX) -c $(CFLAGS) $(LIBS) $(CXXFLAGS) -MMD -MP -MF "rhoMod.o.d" -o rhoMod.o rhoMod.cpp
sharedkld.o:
	$(CXX) -c $(CFLAGS) $(LIBS) $(CXXFLAGS) -MMD -MP -MF "sharedkld.o.d" -o sharedkld.o sharedkld.cpp
sphere_lebedev_rule.o:
	$(CXX) -c $(CFLAGS) $(LIBS) $(CXXFLAGS) -MMD -MP -MF "sphere_lebedev_rule.o.d" -o sphere_lebedev_rule.o sphere_lebedev_rule.cpp


clean:
	rm -rf *.o
	rm -rf *.o.d
	rm -rf $(APPName)
