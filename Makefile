a.out: main.o HISTORY.o SETRHO.o FIELD.o ACCEL.o MOVE.o INIT.o INPUT.o
	mpic++ -g main.o HISTORY.o SETRHO.o FIELD.o ACCEL.o MOVE.o INIT.o INPUT.o -o a.out ../lib/liblua.a -llua -ldl
main.o: main.cpp HISTORY.h SETRHO.h FIELD.h ACCEL.h MOVE.h INIT.h INPUT.h
	mpic++ -g -c main.cpp
HISTORY.o: HISTORY.cpp HISTORY.h
	mpic++ -g -c HISTORY.cpp
SETRHO.o: SETRHO.cpp SETRHO.h
	mpic++ -g -c SETRHO.cpp
FIELD.o: FIELD.cpp FIELD.h FFT.h
	mpic++ -g -c FIELD.cpp
ACCEL.o: ACCEL.cpp ACCEL.h
	mpic++ -g -c ACCEL.cpp
MOVE.o: MOVE.cpp MOVE.h
	mpic++ -g -c MOVE.cpp
INIT.o: INIT.cpp INIT.h
	mpic++ -g -c INIT.cpp
INPUT.o: INPUT.cpp INPUT.h
	mpic++ -g -c INPUT.cpp
clean:
	rm *.o a.out
