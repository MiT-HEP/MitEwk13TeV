GCC=g++
CXXFLAGS=`root-config --libs --cflags` -O2 -fPIC  -I./ -std=c++11
## to use RooUnfold

lxplus=$(findstring lxplus, $(shell hostname -f) )
mit=$(findstring mit.edu, $(shell hostname -f) )

###################### determine if you are on lxplus or mit or else
ifneq ($(strip $(lxplus)),)
$(info You are on lxplus)
#ROOUNFOLD=$(HOME)/work/RooUnfold/
ROOUNFOLD=/afs/cern.ch/user/a/amarini/public/RooUnfold/
else ifneq ($(strip $(mit)),)
		## ON MIT
$(info You are on a mit server)
ROOUNFOLD=/home/amarini/RooUnfold/
else
$(info You are Not on lxplus nor mit)
ROOUNFOLD=/home/kfiekas/CMSSW_7_4_14/src/RooUnfold
CXXFLAGS += -I /home/kfiekas/CMSSW_7_4_14/src/BootStrap -L /home/kfiekas/CMSSW_7_4_14/src/BootStrap/bin -l BootStrap
CXXFLAGES += -Wl,-rpath=/home/kfiekas/CMSSW_7_4_14/src/BootStrap/bin
endif

#CXXFLAGS += -L$(PWD)/../NeroProducer/Core/bin -lBare -Wl,-rpath=$(PWD)/../NeroProducer/Core/bin -ggdb
CXXFLAGS += -L$(ROOUNFOLD)/ -lRooUnfold -ggdb
CXXFLAGS += -I$(ROOUNFOLD)/src/
DICTFLAGS+= -I$(ROOUNFOLD)/src/
SOFLAGS=-shared

SRCDIR=src
BINDIR=bin
HPPDIR=interface

SRC=$(wildcard $(SRCDIR)/*.cxx)
BIN=$(patsubst $(SRCDIR)/%.cxx, $(BINDIR)/% , $(SRC)  )

$(info COMPILING $(SRC))

.PHONY: all
all: $(BIN)


# check if CMSSW is defined
ifndef CMSSW_BASE
$(info No CMSSSW !!!!)
$(info I ll sleep 3s to let you acknowledge it.)
$(info To avoid sleeping 'touch NO_CMSSW')
$(shell [ -f NO_CMSSW ] || sleep 3s)
else
$(info CMSSW found: $(CMSSW_BASE) )
endif

## if on mac add the -p to the  DICTFLAGS
UNAME=$(shell uname)
ifeq ($(UNAME),Darwin)
$(info You are compiling on mac)
DICTFLAGS += -p
else 
$(info Your are on a linux machine)
CXXFLAGS +=  -Wl,-rpath=$(ROOUNFOLD) 
endif

# check if Combine is present and compiled 

#$(OBJ) : $(BINDIR)/%.o : $(SRCDIR)/%.cpp interface/%.hpp | $(BINDIR)
#	$(GCC) $(CXXFLAGS) -c -o $(BINDIR)/$*.o $<

$(BIN) : $(BINDIR)/% : $(SRCDIR)/%.cxx | $(BINDIR)
	$(GCC) $(CXXFLAGS) -o $(BINDIR)/$* $<

$(BINDIR):
	mkdir -p $(BINDIR)

.PHONY: clean
clean:
	-rm $(BIN)
	-rmdir $(BINDIR)
