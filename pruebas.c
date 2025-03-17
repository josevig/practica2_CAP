#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

/* Función para inicializar una matriz con valores aleatorios en [0,1) */
void init_matrix(double *M, int rows, int cols) {
    if (rows <= 0 || cols <= 0) {
        printf("Las dimensiones de la matriz deben ser un valor superior a 0\n");
        exit(1);
    }
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            M[i * cols + j] = (double)rand() / (double)RAND_MAX;
        }
    }
}

/* Función para poner en cero una matriz */
void zero_matrix(double *M, int rows, int cols) {
    if (rows <= 0 || cols <= 0) {
        printf("Las dimensiones de la matriz deben ser un valor superior a 0\n");
        exit(1);
    }
    memset(M, 0, rows * cols * sizeof(double));
}

/* Función para imprimir una matriz (solo si es pequeña) */
void print_matrix(double *M, int rows, int cols) {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            printf("%6.3f ", M[i * cols + j]);
        }
        printf("\n");
    }
    printf("\n");
}

/* Producto de matrices usando la ordenación ijk: C = A * B */
void mult(double *A, double *B, double *C, int m, int k, int n) {
    int i, j, l;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            for (l = 0; l < k; l++) {
                C[i * n + j] += A[i * k + l] * B[l * n + j];
            }
        }
    }
}

/* Suma de matrices: C = A + B */
void sum_matrix(double *A, double *B, double *C, int rows, int cols) {
    int i, j;
    for(i = 0; i < rows; i++){
       for(j = 0; j < cols; j++){
           C[i * cols + j] = A[i * cols + j] + B[i * cols + j];
       }
    }
}

int main(void) {
    int dim;
    printf("Ingrese la dimensión de las matrices cuadradas: ");
    scanf("%d", &dim);

    /* Reservar memoria para las matrices A, B, C, D y E */
    double *A = (double *)malloc(dim * dim * sizeof(double));
    double *B = (double *)malloc(dim * dim * sizeof(double));
    double *C = (double *)malloc(dim * dim * sizeof(double));
    double *D = (double *)malloc(dim * dim * sizeof(double));
    double *E = (double *)malloc(dim * dim * sizeof(double));
    
    if (!A || !B || !C || !D || !E) {
        printf("Error al reservar memoria.\n");
        exit(1);
    }

    /* Inicializar la semilla aleatoria */
    srand((unsigned int)time(NULL));

    /* Inicializar las matrices A y B con valores aleatorios */
    init_matrix(A, dim, dim);
    init_matrix(B, dim, dim);

    /* Inicializar las matrices C, D y E en cero */
    zero_matrix(C, dim, dim);
    zero_matrix(D, dim, dim);
    zero_matrix(E, dim, dim);

    /* Realizar las operaciones:
       C = A * B
       D = C * B
       E = D + C */
    /* Variables para medir el tiempo */
    clock_t start, end;
    double t_mult_AB, t_mult_CB, t_sum;

    /* Medir tiempo para C = A * B */
    start = clock();
    mult(A, B, C, dim, dim, dim);
    end = clock();
    t_mult_AB = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion de C = A * B: %f segundos\n", t_mult_AB);

    /* Medir tiempo para D = C * B */
    start = clock();
    mult(C, B, D, dim, dim, dim);
    end = clock();
    t_mult_CB = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion de D = C * B: %f segundos\n", t_mult_CB);

    /* Medir tiempo para E = D + C */
    start = clock();
    sum_matrix(D, C, E, dim, dim);
    end = clock();
    t_sum = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecucion de E = D + C: %f segundos\n", t_sum);

    /* Imprimir resultados 
    printf("\nMatrix A:\n");
    print_matrix(A, dim, dim);

    printf("Matrix B:\n");
    print_matrix(B, dim, dim);

    printf("Matrix C = A * B:\n");
    print_matrix(C, dim, dim);

    printf("Matrix D = C * B:\n");
    print_matrix(D, dim, dim);

    printf("Matrix E = D + C:\n");
    print_matrix(E, dim, dim); */

    /* Liberar la memoria reservada */
    free(A);
    free(B);
    free(C);
    free(D);
    free(E);

    return 0;
}
