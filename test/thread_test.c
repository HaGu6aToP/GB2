#include "headers.h"
#include <pthread.h>
#include <stdlib.h>

#define N 8
int threads_count = 0;
int executed = 0;

struct Data{
    int id;
    int executed;
    int number;
    pthread_mutex_t* mutex;
    pthread_cond_t* cond;
};

typedef struct Data Data;

void* thread_func(void* params){
    printf("thread_func start \n");

    int k = 0;
    int* p = (int*)params;
    (*p)++;
    for (int i = 0; i < 100000; i++){
        k++;
    }
    (*p)--;
    

    printf("thread_func end, threads_count=%d\n", *p);
}

void* thread_func2(void* params){
    
    Data* data = (Data*)params;
    printf("thread_func2 thread-%d start \n", data->id);

    int wait = 200000 + rand()%200000;
    
    data->number = wait;

    for (int i = 0; i < wait; i++){
        

    }

    // pthread_mutex_lock(data->mutex);
        
    //     if (executed >= 2){
    //         pthread_cond_wait(data->cond, data->mutex);
    //     }

    // pthread_mutex_unlock(data->mutex);

    if (executed >= 2)
        pthread_exit(NULL);

    data->executed = 1;
    executed++;
    printf("thread_func2 thread-%d end, executed=%d\n", data->id, executed);
}



void thread_test(){
    pthread_t* threads = calloc(N, sizeof(pthread_t));
    pthread_attr_t* attr = calloc(N, sizeof(pthread_attr_t));
    Data* data = calloc(N, sizeof(Data));
    pthread_mutex_t lock;
    pthread_cond_t cond;

    pthread_mutex_init(&lock, NULL);
    pthread_cond_init(&cond, NULL);

    printf("thread create\n");
    for (int i = 0; i < N; i++){
        pthread_attr_init(&attr[i]);
        pthread_create(&threads[i], &attr[i], thread_func, &threads_count);
    }

    for(int i = 0; i < N; i++)
        pthread_join(threads[i], NULL);


    printf("------------------------------------------------------\n");

    for (int i = 0; i < N; i++){
        data[i].executed = 0;
        data[i].id = i;
        data[i].number = 0;
        data[i].mutex = &lock;
        data[i].cond = &cond;
        pthread_create(&threads[i], &attr[i], thread_func2, &data[i]);
    }

    
    while(executed < 2){}

    
    printf("end\n");

    for(int i = 0; i < N; i++){
        if (data[i].executed){
            printf("thread-%d is complete. Number is %d\n", data[i].id, data[i].number);
        }
    }

    // return;
    // while(1){}

    free(data);
    free(attr);
    free(threads);
}