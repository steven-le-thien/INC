
CC := gcc -Wall -MP -MD
SRC := src/c

SOURCES := $(wildcard $(SRC)/*.c)
OBJECTS := $(patsubst $(SRC)/%.c,  $(SRC)/%.o, $(SOURCES))

all: $(OBJECTS)
	$(CC) $^ -lm -o ml

$(OBJ)/%.o: $(OBJECTS)
	$(CC) -I$(SRC) -c $< -o $@

clean:
	rm $(SRC)/*.o 
	rm ml
