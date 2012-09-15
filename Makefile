include Makefile.inc
EXEC=genvar
$(EXEC): main.cc *h
	$(CC) -o $(EXEC)  main.cc $(LDFLAGS) $(DEBUGFLAGS)

run:    $(EXEC) 
	 ./$(EXEC) $(SEED) $(R0) > out

submit: $(EXEC) logscale
	 bash submit.sh

logscale: logscale.cc
	 $(CC) -o logscale  logscale.cc $(LDFLAGS) $(DEBUGFLAGS)

bak:

clean:
	rm -f $(EXEC)

	
