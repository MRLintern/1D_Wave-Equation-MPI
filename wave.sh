#!/bash/bin

#opt. flag to be changed here
mpicc -c -O2 -Wall main.c

if [$? -ne 0]; then
	echo "Compile Error"
	exit
fi

mpicc main.o -lm

if [$? -ne 0]; then
	echo "Load Error"
	exit
fi

rm main.o

mv a.out main_w

#change number of processors here; currently set to 4.
#run the application. The results are printed to results.txt
mpirun -np 4 ./main_w > results.txt

if [$? -ne 0]; then
	echo "Run-Time Error"
	exit
fi

rm main_w

echo "Execution Success!"

	
