a.out: main.o HISTORY.o SETRHO.o FIELD.o ACCEL.o MOVE.o
	g++ -g main.o HISTORY.o SETRHO.o FIELD.o ACCEL.o MOVE.o -o a.out -llua -ldl
main.o: main.cpp HISTORY.h SETRHO.h FIELD.h ACCEL.h MOVE.h
	g++ -g -c main.cpp
HISTORY.o: HISTORY.cpp HISTORY.h
	g++ -g -c HISTORY.cpp
SETRHO.o: SETRHO.cpp SETRHO.h
	g++ -g -c SETRHO.cpp
FIELD.o: FIELD.cpp FIELD.h FFT.h
	g++ -g -c FIELD.cpp
ACCEL.o: ACCEL.cpp ACCEL.h
	g++ -g -c ACCEL.cpp
MOVE.o: MOVE.cpp MOVE.h
	g++ -g -c MOVE.cpp
clean:
	rm *.o a.out
