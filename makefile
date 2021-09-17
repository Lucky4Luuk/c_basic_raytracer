CC=gcc
CFLAGS=-std=c99

CFILES = $(wildcard *.c)
OBJ = $(CFILES:.c=.o)
ODIR = obj

.PHONY: clean

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $(ODIR)/$@ $(CFLAGS)

all: $(OBJ)
	@$(CC) $(CFLAGS) -o program $(patsubst %,$(ODIR)/%, $(OBJ))

clean:
	rm -f $(ODIR)/$(OBJ) program
