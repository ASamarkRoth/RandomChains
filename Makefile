CC=g++ -std=gnu++11
CFLAGS=-g -c -Wall
INCLUDES=`root-config --cflags`
LIBDIRS= `root-config --libs --glibs`
ROOTLIBS='-lRooFit -lHtml -lMinuit -lRooFitCore -lRooStats -lHistFactory'
LIBRARY= -L ${ROOTSYS}/lib 
LDFLAGS=
SOURCES=RandomChains.cc 
DEPS=RandomChains.h Style.code
OBJECTS=$(SOURCES:.cc=.o)
EXECUTABLE=run_file


#"executes" dependencies $(SOURCES) and target $(EXECUTABLE)
all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) $(DEPS) 
	$(CC) -o $@ $(OBJECTS) $(LIBDIRS)

.cc.o:
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@
	#$(CC) $(CFLAGS) $< -o $@

#Tells make not to confuse possible clean and help files with the targets with the same names
.PHONY: clean help

help:
	@ echo "Makefile to use with ROOT routines to compile"

clean:
	rm -f $(EXECUTABLE) $(OBJECTS)


#Target which allows you to print variables as "make print-VARIABLE"
print-%  : ; @echo $* = $($*)
