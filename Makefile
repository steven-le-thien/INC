
CC := /usr/local/bin/gcc-8 -Wall -MP -MD -fopenmp
SRC := src/c

SOURCES := $(wildcard $(SRC)/*.c)
OBJECTS := $(patsubst $(SRC)/%.c,  $(SRC)/%.o, $(SOURCES))

all: $(OBJECTS)
	$(CC) $^ -lm -o ml

$(OBJ)/%.o: $(OBJECTS)
	$(CC) -I$(SRC) -c $< -o $@

clean:
	rm $(SRC)/*.o 
	rm $(SRC)/*.d
	rm ml
