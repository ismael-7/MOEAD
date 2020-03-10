# path and directories
SRCDIR		= ./
OBJDIR		= tmp

CXX		= g++ 
LD		= g++
MKDIR	= mkdir -p
RM		= -rm -f

# Source
SOURCEFILES	= SNP.cpp

CPPFILES	= $(SOURCEFILES:%.cpp= $(SRCDIR)/%.cpp)  

# objects needed to generate target
OBJFILES	= $(SOURCEFILES:%.cpp= $(OBJDIR)/%.o)

# options for preprocessor (include directories)
CPPFLAGS	:=-I./


# options for compiler
CXXFLAGS	:= -fopenmp

# options for linker
LDFLAGS	:= -fopenmp
LDLIBS  :=

# make targets
all : dir $(OBJFILES) moea

dir :
	@$(MKDIR) $(OBJDIR)

moea : moea.cpp
	${CXX} ${CXXFLAGS} ${CPPFLAGS} -c moea.cpp -o $(OBJDIR)/moea.o
	${LD} $(LDFLAGS)  -static -o moea $(OBJDIR)/moea.o $(OBJFILES) $(LDLIBS)

# compile .cpp files
$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	${CXX} ${CXXFLAGS} ${CPPFLAGS} -c $< -o $@

.PHONY:all dir clean

# clean up
clean:
		$(RM) $(OBJDIR)/* moea
