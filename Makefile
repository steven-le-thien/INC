
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
	rm $(SRC)/*.gch
	# rm constraint_inc
	rm ml

# ml: inc_ml.o c_inc.o dist.o prim.o options.o tools.o tree.o traversal.o msa.o quartet.o
# 	gcc -Wall inc_ml.o c_inc.o dist.o prim.o options.o tools.o tree.o traversal.o msa.o quartet.o -lm -o ml 
# #	gcc -Wall c_inc.o dist.o prim.o options.o tools.o tree.o traversal.o msa.o quartet.o -o constraint_inc

# c_inc.o: c_inc.c tree.h utilities.h options.h tools.h dist.h prim.h c_inc.h
# 	gcc -Wall -c c_inc.c tree.h options.h tools.h dist.h prim.h traversal.h c_inc.h utilities.h 

# inc_ml.o: inc_ml.c options.h tools.h c_inc.h utilities.h
# 	gcc -Wall -c inc_ml.c options.h tools.h c_inc.h utilities.h 

# dist.o: dist.c dist.h utilities.h options.h tools.h c_inc.h 
# 	gcc -Wall -c dist.c dist.h options.h tools.h c_inc.h utilities.h

# prim.o:  prim.c prim.h options.h utilities.h c_inc.h
# 	gcc -Wall -c prim.c prim.h options.h utilities.h c_inc.h

# options.o: options.c options.h c_inc.h
# 	gcc -Wall -c options.c options.h utilities.h c_inc.h

# tools.o: tools.c tools.h utilities.h options.h c_inc.h
# 	gcc -Wall -c tools.c tools.h options.h utilities.h c_inc.h

# tree.o: tree.c tree.h utilities.h options.h tools.h c_inc.h
# 	gcc -Wall -c tree.c utilities.h options.h tree.h tools.h c_inc.h

# traversal.o: traversal.c traversal.h quartet.h utilities.h options.h c_inc.h
# 	gcc -Wall -c traversal.c traversal.h quartet.h utilities.h options.h c_inc.h

# msa.o: msa.c msa.h utilities.h
# 	gcc -Wall -c msa.c msa.h utilities.h

# quartet.o: quartet.c quartet.h utilities.h
# 	gcc -Wall -c quartet.c quartet.h utilities.h

# clean:
# 	rm *.o 
# 	rm *.gch
# 	# rm constraint_inc
# 	rm ml
