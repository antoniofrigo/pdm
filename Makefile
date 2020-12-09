CC = gfortran
CFLAGS = -o pdm -O3 -fbackslash -fcheck=bounds
OBJ = stellingwerf.f08 pdm.f08
MOD_PATH = MODULES/

all: main

main: $(OBJ)
	@mkdir -p $(MOD_PATH)
	$(CC) -J$(MOD_PATH) $(OBJ) $(CFLAGS) 
