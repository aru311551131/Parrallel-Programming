#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>

//#include "123/WFLCG_c.h"

long thread_count;
long long n;
double sum;
pthread_mutex_t lock;
double max = RAND_MAX / 2;

void* Thread_sum(void* thread) {
  int thr = (int) thread;
  long long pass = 0;
  long long cnt = n / thread_count;
  if (thr == 0)
    cnt += n % thread_count;
  
  /*
  WFLCG_c rng;
  WFLCG_c_init_1_seed(&rng, thr);
  */
  for(int i = 0; i < cnt; i++)
  {
    /*
    float f1 = 2*(WFLCG_c_get_float(&rng) - 1.5f);
    float f2 = 2*(WFLCG_c_get_float(&rng) - 1.5f);
    */
    
    float f1 = rand_r(&thr) / max - 1.f;
    float f2 = rand_r(&thr) / max - 1.f;
    
    if (f1 * f1 + f2 * f2 <= 1)
      pass++;
  }

  pthread_mutex_lock(&lock);
  sum += pass;
  pthread_mutex_unlock(&lock);

}

int main(int argc, char* argv[]) {


  long thread; /* Use long in case of a 64-bit system */
  pthread_t* thread_handles;
  thread_count = strtol(argv[1], NULL, 10);
  n = strtoll(argv[2], NULL, 10);
  thread_handles = (pthread_t*) malloc (thread_count*sizeof(pthread_t)); 
  sum = 0.0;
  pthread_mutex_init(&lock, NULL);

  for (int thread = 0; thread < thread_count; thread++) 
    pthread_create(&thread_handles[thread], NULL, Thread_sum, (void*)thread); 
  
  for (thread = 0; thread < thread_count; thread++) 
    pthread_join(thread_handles[thread], NULL); 
  
  free(thread_handles);
  pthread_mutex_destroy(&lock);

  float pi = 4 * sum / ((float) n);
    printf("%lf\n", pi);

  return 0;
} /* main */
