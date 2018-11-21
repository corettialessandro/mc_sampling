CC	= gcc

CFLAGS	= -Wall -O4
LIBS	= -lm

SRC	= mc_sampling.c
EXEC	= mc_sampling

$(EXEC)	: $(SRC)
	$(CC) $(CFLAGS) $(LIBS) $(SRC) -o $(EXEC)

.PHONY	: clean

clean	:
	/bin/rm -f *.out *.gp *~ $(EXEC)
