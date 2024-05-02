#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <x86intrin.h>

int main() {
	
	#define ARR_SIZE (65536*256*10)
	
	//short* arr = malloc(ARR_SIZE * sizeof(short));
	short* arr = (short*) _aligned_malloc(ARR_SIZE * sizeof(short), 16);
	printf("%ld\n", arr);
	for(size_t i = 0; i < ARR_SIZE; ++i) arr[i] = rand() & 65535;
	printf("filled\n");
	
	size_t cnt = 0;
	short v = 31561;
	
	printf("start\n");
	clock_t start = clock();
	
	__m128i sseVal = _mm_set1_epi16(v);
	for(size_t i = 0; i < ARR_SIZE; i+=64) {
		//if(arr[i] == v) ++cnt;
		/*__m128i sseArr = *(__m128i*) &arr[i];
		cnt += _popcnt32(_mm_movemask_epi8(_mm_cmpeq_epi16(sseVal, sseArr)));*/
		
		__m128i vals = _mm_set_epi32(
			_popcnt32(_mm_movemask_epi8(_mm_cmpeq_epi16(sseVal, *(__m128i*) &arr[i]))),
			_popcnt32(_mm_movemask_epi8(_mm_cmpeq_epi16(sseVal, *(__m128i*) &arr[i + 8]))),
			_popcnt32(_mm_movemask_epi8(_mm_cmpeq_epi16(sseVal, *(__m128i*) &arr[i + 16]))),
			_popcnt32(_mm_movemask_epi8(_mm_cmpeq_epi16(sseVal, *(__m128i*) &arr[i + 24])))
		);
		
		__m128i vals2 = _mm_set_epi32(
			_popcnt32(_mm_movemask_epi8(_mm_cmpeq_epi16(sseVal, *(__m128i*) &arr[i + 32]))),
			_popcnt32(_mm_movemask_epi8(_mm_cmpeq_epi16(sseVal, *(__m128i*) &arr[i + 40]))),
			_popcnt32(_mm_movemask_epi8(_mm_cmpeq_epi16(sseVal, *(__m128i*) &arr[i + 48]))),
			_popcnt32(_mm_movemask_epi8(_mm_cmpeq_epi16(sseVal, *(__m128i*) &arr[i + 56])))
		);
		
		vals = _mm_hadd_epi32(vals, vals2); //8 -> 4
		vals = _mm_hadd_epi32(vals, vals); //4 -> 2
		vals = _mm_hadd_epi32(vals, vals); //2 -> 1
		cnt += _mm_cvtsi128_si32(vals);
	}
	
	cnt /= 2;
	
	long tSpent = (clock() - start);
	printf("measured\n");
	printf("%ld %ld\n", cnt, tSpent);
	
	return 0;
}