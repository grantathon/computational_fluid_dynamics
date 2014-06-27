CC = mpic++#.mpich2 #gcc
CFLAGS = -Wall -Werror
.cpp.o:  ; $(CC) -cpp $(CFLAGS) $<

OBJ = 	main.o\
		Monte_Carlo.o\
		PBM_File.o

all:  $(OBJ)
	$(CC) $(CFLAGS) -o sim_MC $(OBJ)  -lm

%.o : %.cpp
	$(CC) -c $(CFLAGS) $*.cpp -o $*.o

clean:
	rm $(OBJ) *.mc

Monte_Carlo.o : Monte_Carlo.hpp
PBM_File.o : 	PBM_File.hpp
 
main.o        : Monte_Carlo.hpp PBM_File.hpp
