
CC := /usr/local/bin/gcc-8 -Wall -MP -MD -fopenmp
INCMLFLG := -D INC_ML_CMPL
TESTFLG := -D TEST
DIR := src/c

SPECCMPL := 

SOURCES := $(wildcard $(DIR)/*.c)
OBJECTS := $(patsubst $(DIR)/%.c,  $(DIR)/%.o, $(SOURCES))

.PHONY: all clean
all: SPECCMPL = -D INC_ML_CMPL

all: $(OBJECTS)
	$(CC) $(SPECCMPL) $^ -lm -o ml

$(DIR)/%.o: $(DIR)/%.c
	$(CC) $(SPECCMPL) -I$(DIR) -c $< -o $@

clean:
	rm $(DIR)/*.o 
	rm $(DIR)/*.d
	rm ml
