CC = mpic++#.mpich2 #gcc
CFLAGS = -Wall -Werror #-pedantic
.cpp.o:  ; $(CC) -cpp $(CFLAGS) $<

OBJ = 	Monte_Carlo.o\
		main.o

all:  $(OBJ)
	$(CC) $(CFLAGS) -o sim_MC $(OBJ)  -lm

%.o : %.cpp
	$(CC) -c $(CFLAGS) $*.cpp -o $*.o

clean:
	rm $(OBJ) *.mc

Monte_Carlo.o : Monte_Carlo.hpp
 
main.o        : Monte_Carlo.hpp
