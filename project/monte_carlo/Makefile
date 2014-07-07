CC = mpicxx
CFLAGS = -Wall -Werror
.cpp.o:  ; $(CC) -cpp $(CFLAGS) $<

INCFLAGS=../

OBJ = 	Monte_Carlo.o\
		Stochastic_Collocations.o\
		helper.o\
		PBM_File.o\
		main.o

all:  $(OBJ)
	$(CC) $(CFLAGS) -I$(INCFLAGS) -o sim_UQ $(OBJ)  -lm

%.o : %.cpp
	$(CC) -c $(CFLAGS) -I$(INCFLAGS) $*.cpp -o $*.o

clean:
	rm $(OBJ) *.mc

Monte_Carlo.o : Monte_Carlo.hpp
Stochastic_Collocations.o : Stochastic_Collocations.hpp
helper.o : helper.hpp
PBM_File.o : PBM_File.hpp
 
main.o        : Monte_Carlo.hpp Stochastic_Collocations.hpp helper.hpp PBM_File.hpp
