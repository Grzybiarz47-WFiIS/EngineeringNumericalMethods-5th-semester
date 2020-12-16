CC=g++
CPPFLAGS= -O3 -lm -Wall 
OBJS=lab04.o
OUT=lab04

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
	./lab04
