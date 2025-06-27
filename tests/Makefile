CC=mpicc
CFLAGS=-Wall -O2
TARGET=test
OBJS=test.o wave.o

.PHONY: all clean results

all: results $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) -lm

test.o: test.c wave.h
	$(CC) $(CFLAGS) -c test.c

wave.o: wave.c wave.h
	$(CC) $(CFLAGS) -c wave.c

clean:
	rm -f *.o $(TARGET)
	