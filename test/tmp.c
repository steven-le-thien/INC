#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int main(){
	char command[100000];
	float ret;
	float tmp;
	FILE * f;
	ret = 0.0;
	for(int i = 0; i < 10; i++){
		sprintf(command, "FastTree -nt -gtr  ~/INDELibleV1.03/src/outputname_%d.fas > ~/ft_out%d", i + 1, i + 1);
		system(command);

		sprintf(command, "python ~/phylogeny_utilities/compare_trees2.py ~/ft_out%d ~/Downloads/trees_100_leaves.txt > ~/ft_compare%d",  i + 1, i + 1);
		system(command);
		
		sprintf(command, "/Users/lethien_96/ft_compare%d", i + 1);

		printf("%s\n", command);
		f = fopen(command, "r");
		fscanf(f, "%f", &tmp);
		ret += tmp;
		fclose(f);
	}

	printf("final result is %f\n", ret / 10);
	return 0;
}
