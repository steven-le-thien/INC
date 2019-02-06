
CC := /usr/local/bin/gcc-8 -Wall -MP -MD -fopenmp
INCMLFLG := -D INC_ML_CMPL
TESTFLG := -D TEST
DIR := src/c

SPECCMPL := 

SOURCES := $(wildcard $(DIR)/*.c)
OBJECTS := $(patsubst $(DIR)/%.c,  $(DIR)/%.o, $(SOURCES))

.PHONY: ml inc clean
ml: SPECCMPL = -D INC_ML_CMPL
inc: SPECCMPL = -D INC_CMPL 

ml inc: $(OBJECTS)
	$(CC) $(SPECCMPL) $^ -lm -o $@

$(DIR)/%.o: $(DIR)/%.c
	$(CC) $(SPECCMPL) -I$(DIR) -c $< -o $@

clean:
	rm -f $(DIR)/*.o 
	rm -f $(DIR)/*.d
	rm -f ml
	rm -f inc
