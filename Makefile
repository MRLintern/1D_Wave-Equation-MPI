CC = mpicc
CFLAGS = -std=c11 -Wall -Wextra -O2
SRC = main.c
TARGET = bin/main
RESULTS_DIR = results
TEST_SCRIPT = tests/test.sh

.PHONY: all test clean run

all: test $(TARGET)

$(TARGET): $(SRC)
	mkdir -p bin
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET)

test:
	@echo "Running Unit Tests..."
	chmod +x $(TEST_SCRIPT)
	./$(TEST_SCRIPT)

run: $(TARGET)
	mkdir -p $(RESULTS_DIR)
	mpirun -np 2 $(TARGET) > $(RESULTS_DIR)/results.txt

clean:
	rm -rf bin $(RESULTS_DIR)/results.txt wave_test_output.txt
