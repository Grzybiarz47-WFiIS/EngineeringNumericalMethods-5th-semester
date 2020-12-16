CC=g++
CPPFLAGS= -O2 -lm -Wall
OBJS=lab08.o
OUT=lab08

$(OUT): $(OBJS)
	$(CC)  -o $@  $^ $(CPPFLAGS)
		
.c.o:
	$(CC) -c -o $@ $< $(CPPFLAGS) 
	
.PHONY: clean	
clean:
	rm -f *o
	rm -f *~
	rm -f $(OUT)
	rm -f out*.dat
run:
	./lab08
