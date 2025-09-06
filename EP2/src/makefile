CFLAGS=-I. -Wall -Wextra
LIBS=-lm

ODIR=obj
_OBJ=main.o matrix.o
OBJ=$(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	@mkdir -p $(@D)
	gcc -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
