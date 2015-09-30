CFLAGS = -Wall -pedantic -g -O2 -Wno-long-long
INCLUDES = -I/usr/include/gdal
LIBS = -lgdal -lm

SRCS = main.c initialize.c demloader.c

OBJ = $(SRCS:.c=.o)

MAIN = rdrshifter

.PHONY: depend clean install uninstall

all:	$(MAIN)
		@echo "*** $(MAIN) has been compiled. ***"

$(MAIN): $(OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $+ $(LIBS)

%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

$(OBJ): structs.h

clean:
	$(RM) *.o *~ $(MAIN)

#depend: $(SRCS)
#	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
