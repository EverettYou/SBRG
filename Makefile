EXE = fortran_ext
OBJ = ./build
FC = gfortran
TOUCH = touch
RM = rm
LIB = /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/usr/lib/
CFLAGS = -fdefault-real-8 -fimplicit-none -ffree-line-length-none
LFLAGS = -framework Accelerate 

all: $(EXE)
$(EXE): $(OBJ)/$(EXE).o
	@$(FC) -o $(EXE) $(OBJ)/$(EXE).o -L $(LIB) $(LFLAGS) -I $(OBJ)
$(OBJ)/$(EXE).o: $(EXE).f90
	@$(FC) -c $(EXE).f90 -o $(OBJ)/$(EXE).o -L $(LIB) $(CFLAGS) -J $(OBJ)
clean: 
	$(RM) build/*
	$(RM) $(EXE)
