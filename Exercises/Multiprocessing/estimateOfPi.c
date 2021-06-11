#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>

struct params{
	int numberOfPoints;
	int numberInsideCircle;
	unsigned int seed;
};

void* pointThrow(void* args){
	struct params* p = (struct params*) args;
	p->numberInsideCircle = 0;
	double x, y;
	for (int i = 0; i < p->numberOfPoints; i++) {
		x = (double) rand_r(&(p->seed))/RAND_MAX;
		y = (double) rand_r(&(p->seed))/RAND_MAX;
		if (x*x + y*y <= 1) {
			p->numberInsideCircle++;
		}
	}
	return NULL;
}

int main(void){
	printf("=============== Exercise A ===============\n");
	
	int numberOfPoints = 1e8;
	
	struct params paramsThread1 = {.numberOfPoints = numberOfPoints/3.0, .numberInsideCircle = 0, .seed = 3};
	struct params paramsThread2 = {.numberOfPoints = numberOfPoints/3.0, .numberInsideCircle = 0, .seed = 7};
	struct params paramsThread3 = {.numberOfPoints = numberOfPoints/3.0, .numberInsideCircle = 0, .seed = 13};

	pthread_t thread1, thread2, thread3;
	pthread_create(&thread1, NULL, pointThrow, (void*) &paramsThread1);
	pthread_create(&thread2, NULL, pointThrow, (void*) &paramsThread2);
	pthread_create(&thread3, NULL, pointThrow, (void*) &paramsThread3);

	pthread_join(thread1, NULL);
	pthread_join(thread2, NULL);
	pthread_join(thread3, NULL);

	double numberInsideCircle = paramsThread1.numberInsideCircle + paramsThread2.numberInsideCircle + paramsThread3.numberInsideCircle;
	double pi = 4*(double)numberInsideCircle/numberOfPoints;

	printf("Estimate of pi = %g.\n", pi);
	
	return 0;
}
