CC=g++
CPPFLAGS= -O2 -lm -Wall
OBJS=lab06.o mgmres.o
OUT=lab06 

$(OUT): $(OBJS)
	$(CC)  -o $@  $^ $(CPPFLAGS)
		
.c.o:
	$(CC) -c -o $@ $< $(CPPFLAGS) 
	
.PHONY: clean	
clean:
	rm -f *o
	rm -f *~
	rm -f $(OUT)
	rm -f out.dat
run:
	./lab06
