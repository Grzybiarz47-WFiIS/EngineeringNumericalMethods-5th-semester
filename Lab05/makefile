CC=g++
CPPFLAGS= -O2 -lm -Wall
OBJS=lab05.o
OUT=lab05

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
	./lab05
