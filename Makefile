UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),x86_64)
SYSTEM     = x86-64_linux
#------------------------------------------------------------
# Cplex Directorys
#------------------------------------------------------------
CPLEX         = /opt/ibm/ILOG/CPLEX_Studio1261
# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------
CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD
CCLNFLAGS = -lconcert -lilocplex -lcplex -m64 -lm -lpthread -lgmp
else
SYSTEM     = x86_linux
#------------------------------------------------------------
# Cplex Directorys
#------------------------------------------------------------
CPLEX         = /opt/ibm/ILOG/CPLEX_Studio126
# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------
CCOPT = -m32 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD
CCLNFLAGS = -lconcert -lilocplex -lcplex -m32 -lm -lpthread -lgmp
endif
LIBFORMAT  = static_pic

#------------------------------------------------------------
# Cplex & Concert Directorys
#------------------------------------------------------------

CPLEXDIR      = $(CPLEX)/cplex
CONCERTDIR    = $(CPLEX)/concert

# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++ -O0

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)

CONCERTINCDIR	:= $(CONCERTDIR)/include
CPLEXINCDIR	:= $(CPLEXDIR)/include
INCLUDE		:= include

CCFLAGS	:= $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I$(INCLUDE)
#CCFLAGS:= $(CCOPT) -I$(INCLUDE)

# ---------------------------------------------------------------------
#----------------------------------------------------------------------

SRCDIR	:= src
BUILDDIR:= build
SRCEXT	:= cpp
SOURCES	:= $(wildcard $(SRCDIR)/*.$(SRCEXT))
OBJS	:= $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
TARGET	:= bin/SQM

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

$(TARGET): test/main.cpp $(OBJS)
	@echo	" Linking..."
	@echo	"$(CCC) $(CCFLAGS) $(CCLNDIRS) $^ -o $@ $(CCLNFLAGS)"; $(CCC) $(CCFLAGS) $(CCLNDIRS) $^ -o $@ $(CCLNFLAGS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp $(INCLUDE)/%.h
	@mkdir	-p $(BUILDDIR)
	@echo	"$(CCC) $(CCFLAGS) -c -o $@ $(CCLNFLAGS) $<";	$(CCC) $(CCFLAGS) -c $< -o $@ $(CCLNFLAGS)

# Tests
tester: bin/tester

bin/tester: test/tester.cpp $(OBJS)
	$(CCC) $(CCFLAGS) $^ -o $@ $(CCLNDIRS) $(CCLNFLAGS)

clean:
	@echo	" Cleaning..."
	@echo	" rm -r $(BUILDDIR) $(TARGET)"; rm -r $(BUILDDIR) $(TARGET)
# ------------------------------------------------------------
.PHONY: clean


