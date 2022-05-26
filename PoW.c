/* Copyright 2016-2018 The Ulord Core Foundation */

#include "PoW.h"

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
// #include <omp.h>

#include "my_time.h"
#include "common.h"
#include "my_rand48_r.h"
#include "oneWayFunction.h"

// #define SSE_VERSION

void initSpeedyWorkMemory(uint8_t *input, uint32_t inputLen, uint8_t *Maddr, const uint32_t K) {
	uint32_t i, j;
	// uint8_t a[OUTPUT_LEN], b[OUTPUT_LEN];
	// funcInfor[0].func(input, inputLen, a);
	// printf("len is %d , input is : ",inputLen);
	// for(int tt=0;tt<T;tt++)printf("%d,",input[tt*inputLen]);
	// printf("\n");
	uint8_t a[T*OUTPUT_LEN], b[T*OUTPUT_LEN];
	
	for(int tt=0;tt<T;tt++)funcInfor[0].func(&input[tt*inputLen], inputLen, &a[tt*OUTPUT_LEN]);
	// printf("a is : ");
	// for(int tt=0;tt<T;tt++)printf("%d,",input[tt*inputLen]);
	// printf("\n");
// 	uint64_t randSeed[4] = {0, 0, 0, 0};
// #ifndef SSE_VERSION
// 	struct my_rand48_data randBuffer[4];
// #else
// 	struct vrand48_data randBuffer[2];
// #endif

	
	uint64_t randSeed[T*4] = {0};
	struct my_rand48_data randBuffer[T*4];


	const uint32_t iterNum = WORK_MEMORY_SIZE >> 5;
	for (i = 0; i < iterNum; ++i) {
		if (i % K) {
// #ifndef SSE_VERSION
// 			uint64_t num = 0;
// 			for (j = 0; j < 4; ++j) {
// 				my_rand64_r(&randBuffer[j], &num);
// 				memcpy(b + (j << 3), (uint8_t *)&num, 8*sizeof(uint8_t));
// 			}
// #elsev
// 			vrand64(b, randBuffer);
// #endif
			uint64_t num[T] = {0};
			for (j = 0; j < 4; ++j) {
				for(int tt=0;tt<T;tt++)my_rand64_r(&randBuffer[tt*4+j], &num[tt]);
				// printf("i%Kfor_num:");
				// for(int tt=0;tt<T;tt++)printf("%d,",num[tt]);
				// printf("\n");
				for(int tt=0;tt<T;tt++)memcpy(&b[tt*OUTPUT_LEN] + (j << 3), (uint8_t *)&num[tt], 8*sizeof(uint8_t));
				//printf("i%Kfor:");
				//for(int tt=0;tt<T;tt++)printf("%d,",b[tt*OUTPUT_LEN]);
				//printf("\n");
			}

			//memcpy(&b[0] + (0 << 3), (uint8_t *)&num[0], T*4*8*sizeof(uint8_t));
			
			uint8_t shift_num;
			uint8_t result[T*OUTPUT_LEN];
			reduce_bit((uint8_t *)&i, 4, (uint8_t *)&shift_num, 8);
			for(int tt=0;tt<T;tt++)rrs(&b[tt*OUTPUT_LEN], OUTPUT_LEN, &result[tt*OUTPUT_LEN], shift_num);

			for(int tt=0;tt<T;tt++)memcpy(Maddr + tt*WORK_MEMORY_SIZE + (i << 5), &result[tt*OUTPUT_LEN], OUTPUT_LEN*sizeof(uint8_t));
			for (j = 0; j < T*32; ++j) {
				a[j] ^= result[j];
			}
		} else {
			uint8_t t[T] = {0}, shift_num = 0;
			for(int tt=0;tt<T;tt++)reduce_bit(&a[tt*OUTPUT_LEN], 32, (uint8_t *)&t[tt], 8);
			for(int tt=0;tt<T;tt++)t[tt] = (t[tt] & 0x0f) ^ (t[tt] >> 4);
			reduce_bit((uint8_t *)&i, 4, (uint8_t *)&shift_num, 8);
			
			uint8_t a_rrs[T*INPUT_LEN];
			for(int tt=0;tt<T;tt++)rrs(&a[tt*OUTPUT_LEN], OUTPUT_LEN, &a_rrs[tt*INPUT_LEN], shift_num);
			for(int tt=0;tt<T;tt++)funcInfor[t[tt]].func(&a_rrs[tt*INPUT_LEN], 32, &a[tt*OUTPUT_LEN]);
			
			// printf("IK0,t is : ");
			// for(int tt=0;tt<T;tt++)printf("%d,",t[tt]);
			// printf("\n");
			
			for(int tt=0;tt<T;tt++)reduce_bit(&a[tt*OUTPUT_LEN],      8, (uint8_t *)&randSeed[tt*4+0], 48);
			for(int tt=0;tt<T;tt++)reduce_bit(&a[tt*OUTPUT_LEN] +  8, 8, (uint8_t *)&randSeed[tt*4+1], 48);
			for(int tt=0;tt<T;tt++)reduce_bit(&a[tt*OUTPUT_LEN] + 16, 8, (uint8_t *)&randSeed[tt*4+2], 48);
			for(int tt=0;tt<T;tt++)reduce_bit(&a[tt*OUTPUT_LEN] + 24, 8, (uint8_t *)&randSeed[tt*4+3], 48);
#ifndef SSE_VERSION
			for(int tt=0;tt<T;tt++)my_seed48_r(randSeed[tt*4+0], &randBuffer[tt*4+0]);
			for(int tt=0;tt<T;tt++)my_seed48_r(randSeed[tt*4+1], &randBuffer[tt*4+1]);
			for(int tt=0;tt<T;tt++)my_seed48_r(randSeed[tt*4+2], &randBuffer[tt*4+2]);
			for(int tt=0;tt<T;tt++)my_seed48_r(randSeed[tt*4+3], &randBuffer[tt*4+3]);
#else
			vseed48(randSeed    , &randBuffer[0]);
			vseed48(randSeed + 2, &randBuffer[1]);
#endif
			for(int tt=0;tt<T;tt++)memcpy(Maddr + tt*WORK_MEMORY_SIZE + (i << 5), &a[tt*OUTPUT_LEN], 32*sizeof(uint8_t));
		}
	}
}
/* 
 * Step 1: Initialize working memory.
*/
void initWorkMemory(uint8_t *input, uint32_t inputLen, uint8_t *Maddr, const uint32_t K) {
	uint32_t i, j;
	uint8_t a[OUTPUT_LEN], b[OUTPUT_LEN];

	funcInfor[0].func(input, inputLen, a);

	uint64_t randSeed[4] = {0, 0, 0, 0};
#ifndef SSE_VERSION
	struct my_rand48_data randBuffer[4];
#else
	struct vrand48_data randBuffer[2];
#endif

	const uint32_t iterNum = WORK_MEMORY_SIZE >> 5;
	for (i = 0; i < iterNum; ++i) {
		if (i % K) {
#ifndef SSE_VERSION
			uint64_t num = 0;
			for (j = 0; j < 4; ++j) {
				my_rand64_r(&randBuffer[j], &num);
				memcpy(b + (j << 3), (uint8_t *)&num, 8*sizeof(uint8_t));
			}
#else
			vrand64(b, randBuffer);
#endif
			
			uint8_t shift_num;
			uint8_t result[OUTPUT_LEN];
			reduce_bit((uint8_t *)&i, 4, (uint8_t *)&shift_num, 8);
			rrs(b, OUTPUT_LEN, result, shift_num);

			memcpy(Maddr + (i << 5), result, OUTPUT_LEN*sizeof(uint8_t));
			for (j = 0; j < 32; ++j) {
				a[j] ^= result[j];
			}
		} else {
			uint8_t t = 0, shift_num = 0;
			reduce_bit(a, 32, (uint8_t *)&t, 8);
			t = (t & 0x0f) ^ (t >> 4);
			reduce_bit((uint8_t *)&i, 4, (uint8_t *)&shift_num, 8);
			
			uint8_t a_rrs[INPUT_LEN];
			rrs(a, OUTPUT_LEN, a_rrs, shift_num);
			funcInfor[t].func(a_rrs, 32, a);
			
			reduce_bit(a,      8, (uint8_t *)&randSeed[0], 48);
			reduce_bit(a +  8, 8, (uint8_t *)&randSeed[1], 48);
			reduce_bit(a + 16, 8, (uint8_t *)&randSeed[2], 48);
			reduce_bit(a + 24, 8, (uint8_t *)&randSeed[3], 48);
#ifndef SSE_VERSION
			my_seed48_r(randSeed[0], &randBuffer[0]);
			my_seed48_r(randSeed[1], &randBuffer[1]);
			my_seed48_r(randSeed[2], &randBuffer[2]);
			my_seed48_r(randSeed[3], &randBuffer[3]);
#else
			vseed48(randSeed    , &randBuffer[0]);
			vseed48(randSeed + 2, &randBuffer[1]);
#endif
			memcpy(Maddr + (i << 5), a, 32*sizeof(uint8_t));
		}
	}
}

/* 
 * Step 2: Modify the working memory contents.
*/
void modifySpeedyWorkMemory(uint8_t *Maddr, const uint32_t L, const uint32_t C,
		uint8_t *result) {
	uint32_t i, j;
	uint8_t a[T*OUTPUT_LEN], b[T*64];
	
	for(int tt=0;tt<T;tt++)funcInfor[0].func(Maddr + tt * WORK_MEMORY_SIZE + WORK_MEMORY_SIZE - 32, 32, &a[tt*OUTPUT_LEN]);
	for(int tt=0;tt<T;tt++)memcpy(&result[tt*OUTPUT_LEN], &a[tt*OUTPUT_LEN], OUTPUT_LEN*sizeof(uint8_t));
	
	uint64_t r[T] = {0};
	for(int tt=0;tt<T;tt++)reduce_bit(&a[tt*OUTPUT_LEN], 32, (uint8_t *)&r[tt], 64);
	
	const uint32_t iterNum = L << 6;
	for (i = 0; i < C; ++i) {
		// uint64_t randSeed = 0;
		// reduce_bit(a, 32, (uint8_t *)&randSeed, 48);
		
		// struct my_rand48_data randBuffer;
		// my_seed48_r(randSeed, &randBuffer);
		
		// uint8_t t1, t2, s;
		// uint64_t randNum = 0, base = 0;
		uint64_t randSeed[T] = {0};
		for(int tt=0;tt<T;tt++)reduce_bit(&a[tt*OUTPUT_LEN], 32, (uint8_t *)&randSeed[tt], 48);
		
		struct my_rand48_data randBuffer[T];
		for(int tt=0;tt<T;tt++)my_seed48_r(randSeed[tt], &randBuffer[tt]);
		
		uint8_t t1[T], t2[T], s[T];
		uint64_t randNum[T] = {0}, base[T] = {0};
		for (j = 0; j < iterNum; ++j) {
			// my_rand48_r(&randBuffer, &randNum);
			// base = randNum + r;
			for(int tt=0;tt<T;tt++)my_rand48_r(&randBuffer[tt], &randNum[tt]);
			for(int tt=0;tt<T;tt++)base[tt] = randNum[tt] + r[tt];
			
			// uint64_t offset = 0;
			// reduce_bit((uint8_t *)&r, 8, (uint8_t *)&offset, 8);
			// offset = (offset << 8) + 1;
			uint64_t offset[T] = {0};
			for(int tt=0;tt<T;tt++)reduce_bit((uint8_t *)&r[tt], 8, (uint8_t *)&offset[tt], 8);
			for(int tt=0;tt<T;tt++)offset[tt] = (offset[tt] << 8) + 1;
			
			// uint64_t addr1 = (base + WORK_MEMORY_SIZE - offset) % WORK_MEMORY_SIZE;
			// uint64_t addr2 = (base + offset) % WORK_MEMORY_SIZE;
			uint64_t addr1[T],addr2[T];
			for(int tt=0;tt<T;tt++)addr1[tt] = (base[tt] + WORK_MEMORY_SIZE - offset[tt]) % WORK_MEMORY_SIZE + tt*WORK_MEMORY_SIZE;
			for(int tt=0;tt<T;tt++)addr2[tt] = (base[tt] + offset[tt]) % WORK_MEMORY_SIZE + tt*WORK_MEMORY_SIZE;
			
			// t1 = Maddr[addr1];
			// t2 = Maddr[addr2]; 
			// s = a[j & 0x1f];
			for(int tt=0;tt<T;tt++)t1[tt] = Maddr[addr1[tt]];
			for(int tt=0;tt<T;tt++)t2[tt] = Maddr[addr2[tt]]; 
			for(int tt=0;tt<T;tt++)s[tt] = a[tt*OUTPUT_LEN + (j & 0x1f)];
			
			// Maddr[addr1] = t2 ^ s;
			// Maddr[addr2] = t1 ^ s;
			// b[j & 0x3f] = t1 ^ t2;
			for(int tt=0;tt<T;tt++)Maddr[addr1[tt]] = t2[tt] ^ s[tt];
			for(int tt=0;tt<T;tt++)Maddr[addr2[tt]] = t1[tt] ^ s[tt];
			for(int tt=0;tt<T;tt++)b[tt*64 + (j & 0x3f)] = t1[tt] ^ t2[tt];
			
			// r = r + s + t1 + t2;
			for(int tt=0;tt<T;tt++)r[tt] = r[tt] + s[tt] + t1[tt] + t2[tt];
		}

		// uint8_t t = 0;
		// reduce_bit((uint8_t *)&r, 8, (uint8_t *)&t, 8);
		// t = (t & 0x0f) ^ (t >> 4);
		uint8_t t[T] = {0};
		for(int tt=0;tt<T;tt++)reduce_bit((uint8_t *)&r[tt], 8, (uint8_t *)&t[tt], 8);
		for(int tt=0;tt<T;tt++)t[tt] = (t[tt] & 0x0f) ^ (t[tt] >> 4);
		
		// reduce_bit(b, 64, a, 256);
		for(int tt=0;tt<T;tt++)reduce_bit(&b[tt*64], 64, &a[tt*OUTPUT_LEN], 256);
		
		// uint8_t shift_num = 0;
		// uint64_t ir = r + i;
		// reduce_bit((uint8_t *)&ir, 8, (uint8_t *)&shift_num, 8);
		uint8_t shift_num[T] = {0};
		uint64_t ir[T];
		for(int tt=0;tt<T;tt++)ir[tt] = r[tt] + i;
		for(int tt=0;tt<T;tt++)reduce_bit((uint8_t *)&ir[tt], 8, (uint8_t *)&shift_num[tt], 8);

		// uint8_t a_rrs[INPUT_LEN];
		// rrs(a, OUTPUT_LEN, a_rrs, shift_num);
		// funcInfor[t].func(a_rrs, 32, a);
		uint8_t a_rrs[T*INPUT_LEN];
		for(int tt=0;tt<T;tt++)rrs(&a[tt*OUTPUT_LEN], OUTPUT_LEN, &a_rrs[tt*INPUT_LEN], shift_num[tt]);
		for(int tt=0;tt<T;tt++)funcInfor[t[tt]].func(&a_rrs[tt*INPUT_LEN], 32, &a[tt*OUTPUT_LEN]);
		
		// for (j = 0; j < OUTPUT_LEN; ++j) {
		// 	result[j] ^= a[j];
		// }
		for (j = 0; j < T*OUTPUT_LEN; ++j) {
			result[j] ^= a[j];
		}
	}
}
/* 
 * Step 2: Modify the working memory contents.
*/
void modifyWorkMemory(uint8_t *Maddr, const uint32_t L, const uint32_t C,
		uint8_t *result) {
	uint32_t i, j;
	uint8_t a[OUTPUT_LEN], b[64];
	
	funcInfor[0].func(Maddr + WORK_MEMORY_SIZE - 32, 32, a);
	memcpy(result, a, OUTPUT_LEN*sizeof(uint8_t));
	
	uint64_t r = 0;
	reduce_bit(a, 32, (uint8_t *)&r, 64);
	
	const uint32_t iterNum = L << 6;
	for (i = 0; i < C; ++i) {
		uint64_t randSeed = 0;
		reduce_bit(a, 32, (uint8_t *)&randSeed, 48);
		
		struct my_rand48_data randBuffer;
		my_seed48_r(randSeed, &randBuffer);
		
		uint8_t t1, t2, s;
		uint64_t randNum = 0, base = 0;
		for (j = 0; j < iterNum; ++j) {
			my_rand48_r(&randBuffer, &randNum);
			base = randNum + r;
			
			uint64_t offset = 0;
			reduce_bit((uint8_t *)&r, 8, (uint8_t *)&offset, 8);
			offset = (offset << 8) + 1;
			
			uint64_t addr1 = (base + WORK_MEMORY_SIZE - offset) % WORK_MEMORY_SIZE;
			uint64_t addr2 = (base + offset) % WORK_MEMORY_SIZE;
			
			t1 = Maddr[addr1];
			t2 = Maddr[addr2]; 
			s = a[j & 0x1f];
			
			Maddr[addr1] = t2 ^ s;
			Maddr[addr2] = t1 ^ s;
			b[j & 0x3f] = t1 ^ t2;
			
			r = r + s + t1 + t2;
		}
		
		uint8_t t = 0;
		reduce_bit((uint8_t *)&r, 8, (uint8_t *)&t, 8);
		t = (t & 0x0f) ^ (t >> 4);
		
		reduce_bit(b, 64, a, 256);
		
		uint8_t shift_num = 0;
		uint64_t ir = r + i;
		reduce_bit((uint8_t *)&ir, 8, (uint8_t *)&shift_num, 8);

		uint8_t a_rrs[INPUT_LEN];
		rrs(a, OUTPUT_LEN, a_rrs, shift_num);
		funcInfor[t].func(a_rrs, 32, a);
		
		for (j = 0; j < OUTPUT_LEN; ++j) {
			result[j] ^= a[j];
		}
	}
}

/* 
 * Step 3: Calculate the final result.
*/
void calculateFinalResult(uint8_t *Maddr, uint8_t *c, const uint32_t D, uint8_t *result) {
	uint32_t i = 0, j = 0, k = 0;
	memcpy(result, c, OUTPUT_LEN*sizeof(uint8_t));
	
	const uint32_t num = (WORK_MEMORY_SIZE >> 5) - 1;
	
	uint32_t it = 0;
	uint8_t result_rrs[OUTPUT_LEN];
	while(1) {
		uint8_t t = 0, shift_num = 0;
		uint32_t d = 0;
		reduce_bit(result, 32, (uint8_t *)&t, 8);
		t = (t & 0x0f) ^ (t >> 4);
		
		reduce_bit(result, 32, (uint8_t *)&d, D);
		++d;
		
		for (j = 0; j < d; ++j) {
			uint32_t index = i << 5;
			for (k = 0; k < 32; ++k) {
				result[k] ^= Maddr[index + k];
			}
			++i;

			if (i == num) {
				it = i + t;
				reduce_bit((uint8_t *)&it, 4, (uint8_t *)&shift_num, 8);

				rrs(result, OUTPUT_LEN, result_rrs, shift_num);
				funcInfor[0].func(result_rrs, 32, result);
				
				return;
			}
		}
		it = t + i;
		reduce_bit((uint8_t *)&it, 4, (uint8_t *)&shift_num, 8);

		rrs(result, OUTPUT_LEN, result_rrs, shift_num);
		funcInfor[t].func(result_rrs, 32, result);
	}
}
                                                                                                                                                                                                                                                                                                       
/* 
 * Correctness & Performance test for Proof of work
*/
/*
void testPowFunction(uint8_t *mess, uint32_t messLen, const int64_t iterNum) {
	int64_t j;

	uint32_t inputLen = messLen;
	uint8_t input[INPUT_LEN], output[OUTPUT_LEN];
	memset(input, 0, INPUT_LEN*sizeof(uint8_t));
	memcpy(input, mess, messLen*sizeof(char));

	// Init all one-way function
	initOneWayFunction();
	
	uint8_t *Maddr = (uint8_t *)malloc(64 * WORK_MEMORY_SIZE*sizeof(uint8_t));
	assert(NULL != Maddr);
	memset(Maddr, 0, 64 * WORK_MEMORY_SIZE*sizeof(uint8_t));
	
	
	printf("****************************** Correctness test (PoW function) ******************************\n");
	printf("Test message: %s\n", mess);
	powFunction(input, inputLen, Maddr, output);
	view_data_u8("PoW", output, OUTPUT_LEN);
	printf("*********************************************************************************************\n");
	
	printf("*************************************************** Performance test (PoW function) ***************************************************\n");
	uint8_t *result = (uint8_t *)malloc(iterNum * OUTPUT_LEN * sizeof(uint8_t));
	assert(NULL != result);
	memset(result, 0, iterNum * OUTPUT_LEN * sizeof(uint8_t));
	
	uint32_t threadNumArr[] = {1, 4, 8, 12, 16, 20, 24, 32, 48, 64};
	uint32_t threadNumTypes = sizeof(threadNumArr) / sizeof(uint32_t);
	printf("   %-18s", "Algorithm");
	for (uint32_t ix = 0; ix < threadNumTypes; ++ix)
		printf("%12d", threadNumArr[ix]);
	printf("\n");
	
	printf("00 %-18s\t", "PoW");
	for (uint32_t ix = 0; ix < threadNumTypes; ++ix) {
		omp_set_num_threads(threadNumArr[ix]);
		double startTime = get_wall_time();
		if (threadNumArr[ix] == 1) {
			for (j = 0; j < iterNum; ++j) {
				powFunction(input, inputLen, Maddr, result + j * OUTPUT_LEN);
			}
		} else {
			#pragma omp parallel for firstprivate(input), private(j) shared(result)
			for (j = 0; j < iterNum; ++j) {
				powFunction(input, inputLen, Maddr + omp_get_thread_num() * WORK_MEMORY_SIZE, result + j * OUTPUT_LEN);
			}
		}
		double endTime = get_wall_time();
		double costTime = endTime - startTime;
		printf("%5.0f bps   ", iterNum / costTime); fflush(stdout);
		
		// Check result
		for (j = 0; j < iterNum; j += 1) {
			if (memcmp(output, result + j * OUTPUT_LEN, OUTPUT_LEN)) {
				printf("Thread num: %d, j: %ld\n", threadNumArr[ix], j);
				view_data_u8("output", output, OUTPUT_LEN);
				view_data_u8("result", result + j * OUTPUT_LEN, OUTPUT_LEN);
				abort();
			}
		}
	}
	printf("\n");
	printf("***************************************************************************************************************************************\n");
	
	if (NULL != result) {
		free(result);
		result = NULL;
	}
	
	if (NULL != Maddr) {
		free(Maddr);
		Maddr = NULL;
	}
}
*/

#define OUTPUT_BUFFER_SIZE	(32 * 1024UL * 1024UL)
#define MAX_TEST_INPUT_LEN		140
#define MAX_OUT_FILE_NAME_LEN	25

const char testInputCase[][MAX_TEST_INPUT_LEN] = {
	"",
	"HelloWorld",
	"0123456789"
};

void powNistTest(const char *outFileName) {
	const uint64_t iterNum = 1024UL * 1024UL;
	// const uint64_t iterNum = 1024UL;
	
	uint8_t *outputBuffer = (uint8_t *)malloc(OUTPUT_BUFFER_SIZE * sizeof(uint8_t));
	assert(NULL != outputBuffer);
	memset(outputBuffer, 0, OUTPUT_BUFFER_SIZE * sizeof(uint8_t));
	
	uint8_t *Maddr = (uint8_t *)malloc(WORK_MEMORY_SIZE*sizeof(uint8_t));
	assert(NULL != Maddr);
	memset(Maddr, 0, WORK_MEMORY_SIZE*sizeof(uint8_t));
	
	initOneWayFunction();
	
	uint32_t testInputCaseNum = sizeof(testInputCase) / sizeof(const char [MAX_TEST_INPUT_LEN]);
	for (uint32_t testCaseIx = 0; testCaseIx < testInputCaseNum; ++testCaseIx) {
		char curOutFileName[MAX_OUT_FILE_NAME_LEN] = "";
		sprintf(curOutFileName, "%s-%u.txt", outFileName, testCaseIx);
		
		FILE *fp = NULL;
		if (NULL != (fp = fopen(curOutFileName, "wb"))) {
			const uint32_t testInputCaseLen = strlen((char *)testInputCase[testCaseIx]);
			
			uint8_t input[MAX_TEST_INPUT_LEN];
			memset(input, 0, MAX_TEST_INPUT_LEN*sizeof(uint8_t));
			memcpy(input, testInputCase[testCaseIx], testInputCaseLen*sizeof(uint8_t));

			double startTime = get_wall_time();
			powFunction(input, testInputCaseLen, Maddr, outputBuffer);
			for (uint64_t i = 1, j = 0; i < iterNum; ++i) {
				memcpy(input, outputBuffer +  j, OUTPUT_LEN * sizeof(uint32_t));
				j += OUTPUT_LEN;
				
				powFunction(input, OUTPUT_LEN, Maddr, outputBuffer + j);
				
				/* if (j == OUTPUT_BUFFER_SIZE) {
					fwrite(outputBuffer, sizeof(uint8_t), OUTPUT_BUFFER_SIZE / sizeof(uint8_t), fp);
					j = 0;
				} */
			}
			double endTime = get_wall_time();
			double costTime = endTime - startTime;
			fprintf(stdout, "TestCaseIx: %d, Input: %s, IterNum: %lu, Time: %4.2f, Performance: %5.2f bps\n", testCaseIx, \
				testInputCase[testCaseIx], iterNum, costTime, ((double)(iterNum * OUTPUT_LEN)) / costTime); fflush(stdout);

			fwrite(outputBuffer, sizeof(uint8_t), OUTPUT_BUFFER_SIZE / sizeof(uint8_t), fp);

			fclose(fp);
		} else {
			fprintf(stderr, "Error: Open %s failed!\n", curOutFileName);
			abort();
		}
	}
	
	if (NULL != outputBuffer) {
		free(outputBuffer);
		outputBuffer = NULL;
	}

	if (NULL != Maddr) {
		free(Maddr);
		Maddr = NULL;
	}
}

