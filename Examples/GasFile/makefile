OBJDIR = $(GARFIELD_HOME)/Object
SRCDIR = $(GARFIELD_HOME)/Source
INCDIR = $(GARFIELD_HOME)/Include
HEEDDIR = $(GARFIELD_HOME)/Heed
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
CXX = `root-config --cxx`
CFLAGS = `root-config --cflags` \
	-O3 -W -Wall -Wextra -Wno-long-long \
	-fno-common \
	-I$(INCDIR) -I$(HEEDDIR)

# Debug flags
# CFLAGS += -g

LDFLAGS = -L$(LIBDIR) -lGarfield
LDFLAGS += `root-config --glibs` -lGeom -lgfortran -lm

generate: generate.C
	$(CXX) $(CFLAGS) -c generate.C
	$(CXX) $(CFLAGS) -o generate generate.o $(LDFLAGS)
	rm generate.o

read: read.C
	$(CXX) $(CFLAGS) -c read.C
	$(CXX) $(CFLAGS) -o read read.o $(LDFLAGS)
	rm read.o

