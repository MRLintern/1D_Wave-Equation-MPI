CC=mpicc
CCFLAGS=-O3 -Wall
LDFLAGS=
LIBS=-lm

EXE=main
OBJS=main.o

all: $(EXE)

main.o: main.c

$(OBJS): C_COMPILER := $(CC)

$(EXE): $(OBJS)
	$(CC) $(CCFLAGS) $(OBJS) -o $@ $(LDFLAGS) $(LIBS)
	
%.o: %.c
	$(C_COMPILER) $(CCFLAGS) -c $< -o $@
	
.PHONY: clean
clean:
	-/bin/rm -f $(EXE) a.out *.o *
