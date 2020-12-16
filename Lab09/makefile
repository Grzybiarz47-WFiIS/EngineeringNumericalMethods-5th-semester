CC=g++
CPPFLAGS= -O2 -lm -Wall -lgsl -lgslcblas
OBJS=lab09.o
OUT=lab09

$(OUT): $(OBJS)
	$(CC)  -o $@  $^ $(CPPFLAGS)
		
.c.o:
	$(CC) -c -o $@ $< $(CPPFLAGS) 
	
.PHONY: clean	
clean:
	rm -f *o
	rm -f *~
	rm -f $(OUT)
	rm -f *dat
run:
	./lab09
