CC = gfortran
CFLAGS = -o pdm -O3 -fbackslash
OBJ = stellingwerf.f08 pdm.f08
MOD_PATH = MODULES/

all: main

main: $(OBJ)
	$(CC) -J$(MOD_PATH) $(OBJ) $(CFLAGS)
