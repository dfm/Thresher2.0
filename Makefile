CPP = g++
CPPFLAGS = -I/usr/local/include/eigen3 -I/Users/dfm/src/ceres-solver/include
CPPLIB = -L/Users/dfm/src/ceres-solver/ceres-bin/internal/ceres -lceres -lglog -lcxsparse -lcholmod -lamd -lblas -llapack -lcolamd -lprotobuf

SRC = src

TARGETS = $(SRC)/test.o

.cc.o:
	$(CPP) $(CPPFLAGS) -c $< -o $@

thresher: $(TARGETS)
	$(CPP) $(CPPLIB) $(CPPFLAGS) $(TARGETS) -o thresher

clean:
	rm -rf $(TARGETS) thresher
