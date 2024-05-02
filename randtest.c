#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char *argv[])
{	
	int tests = 1 << 10;
	int count[100] = {0};
	
	for(int i=0; i<tests; i++) {
		int summs = 20;
		float value = (float) rand() / (RAND_MAX - 1);
		
		//value /= summs;
		
		count[(int) (value * 100)]++;
	}
	
	for(int y=0; y<20; y++) {
		for(int x=0; x<100; x++) {
			int vl = roundf(20.0f * count[x] * 50 / tests);
			
			printf(vl >= (20 - y) ? "Û" : " ");
		}
		printf("\n");
	}
	
	printf("Tests: %d\n", tests);
	
	return 0;
}