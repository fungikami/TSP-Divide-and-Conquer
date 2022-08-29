# Makefile Proyecto 1

KC=	kotlinc
KFLAGS=	-include-runtime
PROG= TestDACTSP.jar
SRC= DACTSP.kt

all:	
	$(KC) $(SRC) $(KFLAGS) -d $(PROG)

.PHONY : clean

clean :
	rm -rf $(PROG) 
