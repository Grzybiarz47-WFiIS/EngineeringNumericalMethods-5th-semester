CC=g++
CPPFLAGS= -O2 -lm -Wall
OBJS=lab01.o
OUT=lab01

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
	./lab01
