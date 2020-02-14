CC = gcc
CFLAGS = -Wall -g -O3 -ffast-math -fomit-frame-pointer
LDLIBS = -lm
SDIR = ../../src
IDIR = ../../src
ODIR = ../../src
LDFLAGS = -I $(IDIR)
SRC = $(filter-out $(SDIR)/ljmd.c, $(wildcard $(SDIR)/*.c))
SOBJ = $(SRC:.c=.o)

FIL = $(wildcard *.c)
EXE = $(FIL:.c=.x)
OBJ = $(FIL:.c=.o)

all: $(EXE)

#$(SOBJ): $(SRC)
#	$(CC) $(CFLAGS) -c $^

#$(ODIR)/%.o: $(SDIR)/%.c
#	$(CC) $(CFLAGS) -c $<

#$(OBJ): $(FIL)
#	$(CC) $(CFLAGS) $(LDFLAGS) $(LDLIBS) -c $^

#$(EXE): $(OBJ)
#	$(CC) $(CFLAGS) $(SOBJ) $(OBJ) -o $@ $(LDLIBS) $(LDFLAGS)

$(EXE):
	$(CC) $(CFLAGS) $(FIL) $(SRC) -o $(EXE) $(LDLIBS) $(LDFLAGS)
