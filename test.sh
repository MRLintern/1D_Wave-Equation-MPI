#!/bin/bash

# unit testing script
# using it:
# $ chmod +x test.sh
# $ ./test.sh

echo "[INFO] Compiling main.c with mpicc..."
mpicc -o main main.c -lm || { echo "[FAIL] Compilation failed."; exit 1; }

# test parameters
NUM_PROCS=2 # number of processes
EXECUTABLE="./main"
OUTPUT_FILE="wave_test_output.txt"
EXPECTED_STRINGS=("Normal end of execution" "Elapsed wallclock time" "Using 2 processes.")

# clean previous output
rm -f $OUTPUT_FILE

# run the program using MPI with 2 processes
echo "[INFO] Running the Solver with $NUM_PROCS MPI processes..."
mpirun -np $NUM_PROCS $EXECUTABLE > $OUTPUT_FILE
EXIT_CODE=$?

# check exit status
if [$EXIT_CODE -ne 0]; then
	echo "[FAIL] Program exited with code $EXIT_CODE"
	cat $OUTPUT_FILE
	exit 1
fi

# check expected output
echo "[INFO] Validating output..."
for expected in "${EXPECTED_STRINGS[@]}"; do
	if ! grep -q "$expected" "$OUTPUT_FILE"; then
		echo "[FAIL] Expected string not found: '$expected'"
		cat $OUTPUT_FILE
		exit 1
	fi
done

# check for CFL stability
if grep -q "Computation will not be Stable" "$OUTPUT_FILE"; then
	echo "[WARN CFL Condition violated (as expected in some tests)."
fi

echo "[PASS] All tests passed successfully!"
exit 0



