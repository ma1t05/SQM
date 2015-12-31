SYSTEM     = x86_linux
LIBFORMAT  = static_pic

#------------------------------------------------------------
# Cplex & Concert Directorys
#------------------------------------------------------------

CPLEX         = /opt/ibm/ILOG/CPLEX_Studio126
CPLEXDIR      = $(CPLEX)/cplex
CONCERTDIR    = $(CPLEX)/concert

# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++ -O0

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

CCOPT = -m32 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CCLNFLAGS = -lconcert -lilocplex -lcplex -m32 -lm -lpthread -lgmp
#CCLNFLAGS = -m32 -lm -lpthread -lgmp

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
INCLUDE     = include

EXDIR         = $(CPLEXDIR)/examples
EXINC         = $(EXDIR)/include

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I$(INCLUDE)
#CCFLAGS = $(CCOPT) -I$(INCLUDE)

#------------------------------------------------------------

SRCDIR = src
SRCEXT = cpp
SOURCES = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
BUILDDIR= build
OBJS:=$(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
TARGET = bin/SQM

# ------------------------------------------------------------

$(TARGET): $(OBJS)
	@echo	" Linking..."
	@echo	"$(CCC) $(CCLNDIRS) $^ -o $@ $(CCLNFLAGS)"; $(CCC) $(CCLNDIRS) $^ -o $@ $(CCLNFLAGS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp $(INCLUDE)/%.h
	@mkdir	-p $(BUILDDIR)
	@echo	"$(CCC) $(CCFLAGS) -c -o $@ $(CCLNFLAGS) $<";	$(CCC) $(CCFLAGS) -c $< -o $@ $(CCLNFLAGS)

clean:
	@echo	" Cleaning..."
	@echo	" rm -r $(BUILDDIR) $(TARGET)"; rm -r $(BUILDDIR) $(TARGET)
# ------------------------------------------------------------
.PHONY: clean
