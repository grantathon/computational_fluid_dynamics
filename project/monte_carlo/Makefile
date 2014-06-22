CC = mpic++#.mpich2 #gcc
CFLAGS = -Wall -Werror #-pedantic
.cpp.o:  ; $(CC) -cpp $(CFLAGS) $<

OBJ = 	Monte_Carlo.o\
		main.o

all:  $(OBJ)
	$(CC) $(CFLAGS) -o sim $(OBJ)  -lm

%.o : %.cpp
	$(CC) -c $(CFLAGS) $*.cpp -o $*.o

clean:
	rm $(OBJ)

Monte_Carlo.o : Monte_Carlo.hpp
 
main.o        : Monte_Carlo.hpp
