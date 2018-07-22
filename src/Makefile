
decide: main.o msa.o tree.o options.o tools.o stat.o
	gcc -Wall main.o msa.o tree.o options.o tools.o stat.o -o decide

main.o: main.c msa.h tree.h options.h tools.h utilities.h
	gcc -Wall -c main.c msa.h tree.h options.h tools.h utilities.h

msa.o:  msa.c msa.h
	gcc -Wall -c msa.c msa.h tree.h utilities.h

options.o: options.c options.h
	gcc -Wall -c options.c options.h utilities.h

tools.o: tools.c tools.h
	gcc -Wall -c tools.c tools.h

tree.o: tree.c tree.h
	gcc -Wall -c tree.c tree.h

stat.o: stat.c stat.h msa.h utilities.h
	gcc -Wall -c stat.c stat.h msa.h utilities.h

clean:
	rm *.o 
	rm *.gch
	rm decide


new:
	rm double_*
	rm single_*
	rm fasttree.out
