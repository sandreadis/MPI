        /*KNN-Algorithm*/
/***Find the k nearest neighbors***/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>

typedef struct{
    double* distances;
    int* ids;
}knnresult;

void swap(double *a, double *b);
void swap_int(int *a, int *b);
int partition(double *array, int* index_array, int low, int high);
void quickSort(double *array, int* index_array, int low, int high);
double** find_distances(double** array_1, int size_1, double** array_2, int size_2, int dimensions);
void create_X(double** X, int size, int dimensions);
void routine(double **query, int query_size, double **corpus, int corpus_id,  int corpus_size, double** di, int* index_array, knnresult *kNN, int k, int size, int right_size, int processes);
void write_to_file(knnresult *knn_result, int size, int k);
void read_file(double** X, int size, int dimensions);

int main(int argc, char** argv){

    struct timeval start;	/* starting time */
	struct timeval end;	/* ending time */
	unsigned long e_usec;	/* elapsed microseconds */
    gettimeofday(&start, 0);

    double **X;
    double **di;
    int* index_array;
    double* distances_copy;
    int* ids_copy;
    knnresult *knn_result;
   

    int size, processes, dimensions, right_size, false_size, k;
    //test
    processes = 4;
    if(atoi(argv[1]) == 2){
        size = 10000; 
        dimensions = 2; 
        k = 9;
    }
    else if (atoi(argv[1]) == 3)
    {
        size = 1000; 
        dimensions = 3; 
        k = 27;
    }
    else if (atoi(argv[1]) == 4)
    {
        size = 100000000; 
        dimensions = 4; 
        k = 81;
    }
    else{
        printf("Didn't give a valid dimension\n");
        exit(1);
    }
    
    knn_result = malloc(size * sizeof(knnresult));
    X = malloc(size * sizeof(double*));
    for(int i = 0; i < size; i++){
        X[i] = malloc(dimensions * sizeof(double));
        knn_result[i].distances = malloc(k * sizeof(double));
        knn_result[i].ids = malloc(k * sizeof(int));
    }
    //create_X(X, size, dimensions);
    read_file(X, size, dimensions);

    right_size = size/processes;
    false_size = size/processes + size%processes;

    int repetition_x = 0;
    int repetition_y = 0;
    int start_x;
    int end_x;
    int start_y;
    int end_y;
    int size_x;
    int size_y;

    int counter;
    int id_1, id_2;

    repetition_x = 0;
    while(repetition_x < processes){
        if(repetition_x == 0){
            start_x = 0;
            end_x = false_size;
        }
        else{
            start_x = false_size + (repetition_x - 1) * right_size;
            end_x = start_x + right_size;
        }
        repetition_y = 0;
        size_x = end_x - start_x;
        while(repetition_y < processes){
            if(repetition_y == 0){
                start_y = 0;
                end_y = false_size;
            }
            else{
                start_y = false_size + (repetition_y - 1) * right_size;
                end_y = start_y + right_size;
            }
            size_y = end_y - start_y;
            di = malloc(size_x * sizeof(double*));
            for(int i = 0; i < size_x; i++){
                di[i] = calloc(size_y, sizeof(double));
            }
            index_array = malloc(size_y * sizeof(int)); 
            
            for(int i = 0; i < size_x; i++){
                for(int j = 0; j < size_y; j++){
                    for(int l = 0; l < dimensions; l++){
                        di[i][j] += X[i + start_x][l] * X[i + start_x][l] + X[j + start_y][l] * X[j + start_y][l] - 2 * (X[i + start_x][l]*X[j + start_y][l]);
                    }
                    di[i][j] = sqrt(di[i][j]);
                }
            }
            if(repetition_y == 0){   
                for(int i = 0; i < size_x; i++){
                    for(int j = 0; j < size_y; j++){
                        index_array[j] = start_y + j;
                    }
                    quickSort(di[i], index_array, 0, size_y - 1);
                    for(int l = 0; l < k; l++){
                        knn_result[i + start_x].distances[l] = di[i][l];
                        knn_result[i + start_x].ids[l] = index_array[l];
                    }
                }
            }
            else{
                for(int i = 0; i < size_x; i++){
                    for(int j = 0; j < size_y; j++){
                        index_array[j] = start_y + j;
                    }
                    quickSort(di[i], index_array, 0, size_y - 1);
                    distances_copy = knn_result[i + start_x].distances;
                    ids_copy = knn_result[i + start_x].ids;
                    knn_result[i + start_x].ids = NULL;
                    knn_result[i + start_x].distances = NULL;
                    knn_result[i + start_x].ids = malloc(k * sizeof(int));
                    knn_result[i + start_x].distances = malloc(k * sizeof(double));
                    counter = 0;
                    id_1 = 0;
                    id_2 = 0;
                    while(counter < k){
                        if(distances_copy[id_2] > di[i][id_1]){
                            knn_result[i + start_x].distances[counter] = di[i][id_1];
                            knn_result[i + start_x].ids[counter] = index_array[id_1];
                            id_1++;
                            counter++;
                        }
                        else{
                            knn_result[i + start_x].distances[counter] = distances_copy[id_2];
                            knn_result[i + start_x].ids[counter] = ids_copy[id_2];
                            id_2++;
                            counter++;
                        }
                    }
                    free(distances_copy);
                    free(ids_copy);
                }
            }
            for(int i = 0; i < size_x; i++) free(di[i]);
            free(di);
            free(index_array);
            repetition_y++;
        }
        repetition_x++;
    }

    write_to_file(knn_result, size, k);
    gettimeofday(&end, 0);
    e_usec = ((end.tv_sec * 1000000) + end.tv_usec) - ((start.tv_sec * 1000000) + start.tv_usec);
    printf("time: %lf\n", e_usec / 1000000.0);

    return 0; 

}


void read_file(double** X, int size, int dimensions){
    
    FILE *fptr;
    switch (dimensions){
        case 2:
            fptr = fopen("data/regular_grid_2d_10000.txt", "r");
            break;
        case 3:
            fptr = fopen("data/regular_grid_3d.txt", "r");
            break;
        case 4:
            fptr = fopen("data/regular_grid_4d.txt", "r");
            break;
    }

    if(fptr == NULL){
        printf("Error! opening file\n");

        // Program exits if the file pointer returns NULL.
        exit(1);
    }
    for(int i = 0; i < size; i++){
        for(int d = 0; d < dimensions; d++){
            fscanf(fptr,"%lf ", &X[i][d]);
        }
        fscanf(fptr, "\n");
    }
}
void write_to_file(knnresult *knn_result, int size, int k){

    FILE *fptr;
    fptr = fopen("database/results.txt", "w");

    if(fptr == NULL){
        printf("Error! opening file\n");

        // Program exits if the file pointer returns NULL.
        exit(1);
    }
    for(int i = 0; i < size; i++){
        fprintf(fptr, "point%d neighbors:\n", i);
        for(int d = 0; d < k; d++){
            fprintf(fptr,"%d ", knn_result[i].ids[d]);
        }
        fprintf(fptr, "\n");
    }

    fclose(fptr);

}
void routine(double **query, int query_size, double **corpus, int corpus_id,  int corpus_size, double** di, int* index_array, knnresult *kNN, int k, int size, int right_size, int processes){

    index_array = malloc(corpus_size * sizeof(int));
    int counter;
    int id_1, id_2, id_total;
    double* distances_copy;
    int* ids_copy;

    for(int i = 0; i < query_size; i++){
        counter = 0;
        id_1 = 0;
        id_2 = 0;       
        id_total = 0;
        for(int j = 0; j < corpus_size; j++){
            if(corpus_id == 0){
                index_array[j] = j;
            }
            else{
                index_array[j] = j + size/processes + size%processes + (corpus_id - 1) * right_size;
            }
        }
        quickSort(di[i], index_array, 0, corpus_size - 1);
        distances_copy = kNN[i].distances;
        ids_copy = kNN[i].ids;
        kNN[i].distances = NULL;
        kNN[i].ids = NULL;
        kNN[i].ids = malloc(k * sizeof(int));
        kNN[i].distances = malloc(k * sizeof(double));
        while(counter < k){
            if(distances_copy[id_2] > di[i][id_1]){
                kNN[i].distances[id_total] = di[i][id_1];
                kNN[i].ids[id_total] = index_array[id_1];
                id_1++;
                id_total++;
            }
            else{
                kNN[i].distances[id_total] = distances_copy[id_2];
                kNN[i].ids[id_total] = ids_copy[id_2];
                id_2++;
                id_total++;
            }
            counter++;
        }
        free(distances_copy);
        free(ids_copy);
    }
    for(int i = 0; i < query_size; i++) free(di[i]);
    free(di);
    free(index_array);

}
void create_X(double** X, int size, int dimensions){

    // X = malloc(size * sizeof(double*));
    // for(int i = 0; i < size; i++){
    //     X[i] = malloc(dimensions * sizeof(double));
    // }


    int range = cbrt(size);
    int step = 1;
    int id = 0;


    for(int i = 0; i < range; i++){
        for(int j = 0; j < range; j++){
            for(int l = 0; l < range; l++){
                X[id][0] = i;
                X[id][1] = j; 
                X[id][2] = l;
                id++;
            }
        }
    }

    // for(int i = 0; i < size; i++){
    //     for(int j = 0; j < dimensions; j++)
    //         printf("%lf ", X[i][j]);
    //     printf("\n");
    // }


}
double** find_distances(double** array_1, int size_1, double** array_2, int size_2, int dimensions){

    double** d = malloc(size_1 * sizeof(double*)); 
    for(int i = 0; i < size_1; i++){
        d[i] = calloc(size_2, sizeof(double));
    }

    for(int i = 0; i < size_1; i++){
        for(int j = 0; j < size_2; j++){
            for(int l = 0; l < dimensions; l++){
                d[i][j] += array_1[i][l] * array_1[i][l] + array_2[j][l] * array_2[j][l] - 2 * (array_1[i][l]*array_2[j][l]);
            }
            d[i][j] = sqrt(d[i][j]);
        }
    }
    return d;
}
void swap(double *a, double *b){
    double t = *a;
    *a = *b;
    *b = t;
}
void swap_int(int *a, int *b){
    int t = *a;
    *a = *b;
    *b = t;
}
int partition(double *array, int* index_array, int low, int high){

    double pivot = array[high];
    //printf("pivot: %lf\n", pivot);
    int i = low - 1;

    for(int j = low; j < high; j++) {
        if (array[j] <= pivot) {
            i++;
            //printf("i, j: %d %d\n", i, j);
            swap(&array[i], &array[j]);
            swap_int(&index_array[i], &index_array[j]);
        }
    }
    //print_array(array, index_array, high + 1);
    //printf("i+1, high: %d %d\n", i+1, high);
    swap(&array[i + 1], &array[high]);
    swap_int(&index_array[i + 1], &index_array[high]);


    return (i + 1); 
}
void quickSort(double *array, int* index_array, int low, int high) {
    if (low < high) {
        int pi = partition(array, index_array, low, high);
        //print_array(array, index_array, high);
        //printf("pi: %d\n", pi);
        quickSort(array, index_array, low, pi - 1);
        quickSort(array, index_array, pi + 1, high);
    }
}

