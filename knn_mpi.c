#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
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

    //mpi-variables
    int rank, processes, next, previous, corpus_id;
    //X kNN
    int size, dimensions, k;
    //sizes
    int query_size, corpus_size, new_corpus_size, right_size, false_size;
    
    //matrices
    double **X, **xi, **yi, **zi, **di;
    int* index_array;   
    knnresult *kNN, *kNN_result;
    
    //test
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
        size = 4096; 
        dimensions = 4; 
        k = 81;
    }
    else{
        printf("Didn't give a valid dimension\n");
        exit(1);
    }
    

    //mpi-start
    MPI_Init(&argc, &argv);

    //init variables
    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    next = (rank + 1) % processes;
    previous = (rank + processes - 1) % processes;
    right_size = size/processes;
    if(rank == 0){
        query_size = size/processes + size%processes;
    }
    else{
        query_size = size/processes;
    }
    //printf("%d of %d\nsize:%d\ndimensions:%d\nk:%d\nprevious:%d\nnext:%d\nquery_size:%d\n\n", rank, processes, size, dimensions, k,previous, next, query_size);
    
    //initialize kNN and xi
    kNN = malloc(query_size * sizeof(knnresult));
    xi = malloc(query_size * sizeof(double*));
    //error checking
    if(!xi || !kNN){
        printf("Error in memory allocation.\n");
        return(1);
    }
    for(int i = 0; i < query_size; i++){
        kNN[i].distances = malloc(k * sizeof(double));
        kNN[i].ids = malloc(k * sizeof(int));
        xi[i] = malloc(dimensions * sizeof(double));
        if(!xi[i] || !kNN[i].distances || !kNN[i].ids){
            printf("Error in memory allocation.\n");
            return(1);
        }
    }

    //initialize xi - split data
    if(rank == 0){
        
        // create X
        X = malloc(size * sizeof(double*));
        if(!X){
            printf("Error in memory allocation.\n");
            return(1);
        }
        for(int i = 0; i < size; i++){
            X[i] = malloc(dimensions * sizeof(double));
            if(!X[i]){
                printf("Error in memory allocation.\n");
                return(1);
            }
        }

        //test
        //create_X(X, size, dimensions);
        read_file(X, size, dimensions);
 

        //MPI Send
        int id;
        for(int j = 1; j < processes; j++){
            for(int i = 0; i < right_size; i++){
                id = query_size + (j-1) * right_size;
                MPI_Send(X[id + i], dimensions, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
            }
        }
        //init xi0
        for(int i = 0; i < query_size; i++){
            for(int j = 0; j < dimensions; j++)    
                xi[i][j] = X[i][j];
        }
        //find k with my database for process 0
        di = find_distances(xi, query_size, xi, query_size, dimensions);
        corpus_size = query_size;
        index_array = malloc(corpus_size * sizeof(int)); //keep indices
        for(int i = 0; i < query_size; i++){
            for(int j = 0; j < corpus_size; j++){
                index_array[j] = j;
            }
            quickSort(di[i], index_array, 0, corpus_size - 1);   
            for(int l = 0; l < k; l++){
                kNN[i].distances[l] = di[i][l];
                kNN[i].ids[l] = index_array[l];
            }
        }
        //prepare for ring
        zi = xi;    //send corpus
        corpus_id = rank;   //know your id without sending matrix

        //free
        for(int i = 0; i < query_size; i++) free(di[i]);
        free(di);
        free(index_array);
        for(int i = 0; i < size; i++) free(X[i]);
        free(X);
    }
    else{

        //MPI Receive
        for(int j = 0; j < query_size; j++){
            MPI_Recv(xi[j], dimensions, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        
        //find k with my database for processes except 0
        di = find_distances(xi, query_size, xi, query_size, dimensions);
        corpus_size = query_size;
        index_array = malloc(corpus_size * sizeof(int)); //keep indices
        for(int i = 0; i < query_size; i++){
            for(int j = 0; j < corpus_size; j++){
                index_array[j] = (rank-1) * right_size + 
                (right_size + size%processes) + j;  //global index
            }
            quickSort(di[i], index_array, 0, corpus_size - 1);
            for(int l = 0; l < k; l++){
                kNN[i].distances[l] = di[i][l];
                kNN[i].ids[l] = index_array[l];
            }
        }
        //prepare for ring
        zi = xi; //send corpus
        corpus_id = rank;

        //free
        for(int i = 0; i < query_size; i++) free(di[i]);
        free(di);
        free(index_array);
    }


    //ring-start
    int repetition = 0;
    while(repetition < processes - 1){  //we already did the first
        
        //Send info for reception
        int temp_corpus_id = corpus_id;
        MPI_Send(&temp_corpus_id, 1, MPI_INT, next, 0, MPI_COMM_WORLD);
        MPI_Recv(&corpus_id, 1, MPI_INT, previous, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&corpus_size, 1, MPI_INT, next, 0, MPI_COMM_WORLD);    
        MPI_Recv(&new_corpus_size, 1, MPI_INT, previous, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        //Send corpus to next
        for(int i = 0; i < corpus_size; i++){
            MPI_Send(zi[i], dimensions, MPI_DOUBLE, next, 0, MPI_COMM_WORLD);
        
        }
        //Receive new corpus
        corpus_size = new_corpus_size;
        yi = malloc(corpus_size * sizeof(double*));
        if(!yi){
            printf("Error in memory allocation.\n");
            return(1);
        }
        for(int i = 0; i < corpus_size; i++){
            yi[i] = malloc(dimensions * sizeof(double));
            if(!yi[i]){
                printf("Error in memory allocation.\n");
                return(1);
            }   
        }
        for(int i = 0; i < corpus_size; i++){
            MPI_Recv(yi[i], dimensions, MPI_DOUBLE, previous, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        di = find_distances(xi, query_size, yi, corpus_size, dimensions);
        //routine
        routine(xi, query_size, yi, corpus_id, corpus_size, di, index_array, kNN, k, size, right_size, processes);

        //change pointers for next routine
        zi = yi;
        yi = NULL;

        repetition++;
    }
    for(int i = 0; i < corpus_size; i++){
        free(zi[i]);
    }
    free(zi);
    for(int i = 0; i < query_size; i++){
        free(xi[i]);
    }
    free(xi);


    if(rank != 0){
        //Send results to process 0
        for(int i = 0; i < query_size; i++){
            MPI_Send(kNN[i].distances, k, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(kNN[i].ids, k, MPI_INT, 0, 1, MPI_COMM_WORLD);
        }
    }
    else{
        //Receive result and create the results
        kNN_result = malloc(size * sizeof(knnresult));
        if(!kNN_result){
            printf("Error in memory allocation.\n");
            return(1);           
        }
        for(int i = 0; i < size; i++){
            for(int j = 0; j < size; j++){
                kNN_result[i].distances = malloc(k * sizeof(double));
                kNN_result[i].ids = malloc(k * sizeof(int));
                if(!kNN_result[i].distances || !kNN_result[i].ids){
                    printf("Error in memory allocation.\n");
                    return(1);
                }
            }
        }
        for(int i = 0; i < query_size; i++){
            kNN_result[i].distances = kNN[i].distances;
            kNN_result[i].ids = kNN[i].ids;
        }
        for(int i = 1; i < processes; i++){
            for(int j = 0; j < right_size; j++){
                MPI_Recv(kNN_result[j + (i - 1) * right_size + query_size].distances, k, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(kNN_result[j + (i - 1) * right_size + query_size].ids, k, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

    
    //print results
    if(rank == 0){ 
        write_to_file(kNN_result, size, k);
        gettimeofday(&end, 0);
        e_usec = ((end.tv_sec * 1000000) + end.tv_usec) - ((start.tv_sec * 1000000) + start.tv_usec);
        printf("time: %lf\n", e_usec / 1000000.0);

    }

    MPI_Finalize();


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
    int i = low - 1;

    for(int j = low; j < high; j++) {
        if (array[j] <= pivot) {
            i++;
            swap(&array[i], &array[j]);
            swap_int(&index_array[i], &index_array[j]);
        }
    }
    swap(&array[i + 1], &array[high]);
    swap_int(&index_array[i + 1], &index_array[high]);


    return (i + 1); 
}
void quickSort(double *array, int* index_array, int low, int high) {
    if (low < high) {
        int pi = partition(array, index_array, low, high);
        quickSort(array, index_array, low, pi - 1);
        quickSort(array, index_array, pi + 1, high);
    }
}