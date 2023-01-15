#include <stdio.h>
#include <stdlib.h>


void create_X(int** X, int size, int dimensions);

int main(int argc, char** argv){

    
    int size = 100000000;
    int dimensions = 4;
    int** X = malloc(size * sizeof(int*));
    for(int i = 0; i < size; i++){
        X[i] = malloc(dimensions * sizeof(int));
    }
    create_X(X, size, dimensions);
    FILE *fptr;
    fptr = fopen("data/regular_grid_4d_100M.txt", "w");
    if(fptr == NULL)
    {
      printf("Error!\n");   
      exit(1);             
    }

    for(int i = 0; i < size; i++){
        for(int d = 0; d < dimensions; d++){
            fprintf(fptr,"%d ", X[i][d]);
        }
        fprintf(fptr, "\n");
    }

    fclose(fptr);
    return 0;
}


void create_X(int** X, int size, int dimensions){

    // X = malloc(size * sizeof(double*));
    // for(int i = 0; i < size; i++){
    //     X[i] = malloc(dimensions * sizeof(double));
    // }


    int range = 100;
    int step = 1;
    int id = 0;


    for(int i = 0; i < range; i++){
        for(int j = 0; j < range; j++){
            for(int l = 0; l < range; l++){
                for(int d = 0; d < range; d++){
                    
                    X[id][0] = i;
                    X[id][1] = j; 
                    X[id][2] = l;
                    X[id][3] = d;
                    id++;
                }
            }
        }
    }

    // for(int i = 0; i < size; i++){
    //     for(int j = 0; j < dimensions; j++)
    //         printf("%lf ", X[i][j]);
    //     printf("\n");
    // }


}