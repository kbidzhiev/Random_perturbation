#################################################################
#################################################################
#################################################################


include ./this_dir.mk
include ./options.mk


#Define Flags ----------

TENSOR_HEADERS=$(PREFIX)/itensor/core.h
CCFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(CPPFLAGS) $(OPTIMIZATIONS)
CCGFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS)
LIBFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBFLAGS)
LIBGFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBGFLAGS)



#Rules ------------------

%.o: %.cpp $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cpp $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

CPP_FILES = Random_perturbation.cpp tech_param.cpp time_evolution.cpp Ham_Folded_XXZ.cpp XXZ.cpp observables.cpp quantum_measurement.cpp

Random_perturbation.exe: ${CPP_FILES} $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) ${CPP_FILES} -o $@ $(LIBFLAGS)


all: Random_perturbation.exe
clean:
	rm *.exe *.o

