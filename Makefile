include Makefile.inc
EXEC=genvar
$(EXEC): main.cc *h
	$(CC) -o $(EXEC)  main.cc $(LDFLAGS) $(DEBUGFLAGS)


run:	$(EXEC) 
	./$(EXEC) $(SEED) > out


bak:

clean:
	rm -f $(EXEC)

	
