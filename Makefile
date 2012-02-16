HDF5LIB = /usr/lib/libhdf5.a -lz -lpthread
CXXFLAGS = -Wall -O3 -I/home/bruno/Documents/vendors/silo-4.8-bsd/include
LIBFLAGS = -L/home/bruno/Documents/vendors/silo-4.8-bsd/lib -lsiloh5 $(HDF5LIB) -lm

apollo : apollo.o POST_PROCESSING.o
	$(CXX) $(CXXFLAGS) -o apollo apollo.o POST_PROCESSING.o $(LIBFLAGS) 

apollo.o : apollo.cc
	$(CXX) $(CXXFLAGS) -c apollo.cc

POST_PROCESSING.o : POST_PROCESSING.cc POST_PROCESSING.hh
	$(CXX) $(CXXFLAGS) -c POST_PROCESSING.cc

.PHONY : clean
clean :
	-rm *.o
	-rm apollo
