#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>

#pragma warning(disable : 4996)
#pragma warning(disable : 4700)
using namespace std;
#define NUM_DIMS 1
#define SIZE 5

int recv(int source) {
    int value;
    MPI_Status status;
    MPI_Recv(&value, 1, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    return value;
}

void send(int dest, int val) {
    MPI_Send(&val, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
}

void sendTo0Process(const int sendData, int *recvbuf, int rank) {
    if (rank == 0) {
        recvbuf[0] = sendData;
        for (int i = 1; i < SIZE; i++) {
            recvbuf[i] = recv(i);
        }

    } else {
        send(0, sendData);
    }

}

void getFrom0Process(const int *sendData, int count, int *recvbuf, int rank) {
    if (rank == 0) {
        for (int i = 0; i < SIZE; i++) {
            for (int j = i * count; j < (i + 1) * count; j++) {
                send(i, sendData[j]);
            }
        }

    }
    for (int i = 0; i < count; i++) {
        recvbuf[i] = recv(0);
    }

}


void ring(int *matrix, int *vector, int *result, int rank, MPI_Comm comm) {
    int source, dest;
    int matrix_str[SIZE];
    int vector_elem;
    char elements;
    getFrom0Process(matrix, SIZE, matrix_str, rank);
    getFrom0Process(vector, 1, &vector_elem, rank);
    MPI_Cart_shift(comm, 0, 1, &source, &dest);
    int element = 0;
    for (int i = 0; i < SIZE; i++) {
        element += matrix_str[(rank + SIZE - i) % SIZE] * vector_elem;
        if (rank == 0) {
            send(dest, vector_elem);
            vector_elem = recv(SIZE - 1);
        } else {
            int temp = vector_elem;
            vector_elem = recv(source);
            send(dest, temp);
        }
    }
    sendTo0Process(element, result, rank);
    if (rank == 0) {
        printf("Результирующий вектор при использовании 'кольцевой' топологии:\n");
        for (int i = 0; i < SIZE; ++i) {
            printf("%d\t", result[i]);
        }
        printf("\n");
    }
}

void line(int *matrix, int *vector, int *result, int rank, MPI_Comm comm) {
    int source, dest;
    int matrix_str[SIZE];

    getFrom0Process(matrix, SIZE, matrix_str, rank);

    MPI_Cart_shift(comm, 0, 1, &source, &dest);

    int element = 0;
    for (int i = 0; i < SIZE; ++i) {
        int vector_elem = rank == 0 ? vector[i] : recv(source);
        element += matrix_str[i] * vector_elem;

        if (dest > 0) {
            send(dest, vector_elem);
        }
    }

    sendTo0Process(element, result, rank);

    if (rank == 0) {
        printf("Результирующий вектор при использовании 'линейной' топологии:\n");
        for (int i = 0; i < SIZE; ++i) {
            printf("%d\t", result[i]);
        }
        printf("\n");
    }
}


int main(int argc, char **argv) {
    int rank, size, dims[NUM_DIMS];
    int periods[NUM_DIMS];
    int reorder = 0;
    int vector[SIZE] = {1, 2, 3, 4, 5};
    int result[SIZE];
    int matrix[
            SIZE * SIZE] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};

    MPI_Comm comm;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {

        printf("Исходная матрица:\n");
        for (int i = 0; i < SIZE; i++) {
            for (int j = 0; j < SIZE; j++) {
                printf("%d\t", matrix[SIZE * i + j]);
            }
            printf("\n");
        }

        printf("Исходный вектор:\n");
        for (int i = 0; i < SIZE; ++i) {
            printf("%d\t", vector[i]);
        }
        printf("\n");
    }
    bool isRingMode = atoi(argv[1]);

    if (isRingMode) {
        for (int i = 0; i < NUM_DIMS; i++) {
            dims[i] = SIZE;
            periods[i] = 1;
        }
    } else {
        for (int i = 0; i < NUM_DIMS; i++) {
            dims[i] = SIZE;
            periods[i] = 0;
        }
    }

    MPI_Dims_create(size, NUM_DIMS, dims);

    MPI_Cart_create(MPI_COMM_WORLD, NUM_DIMS, dims, periods, reorder,&comm);
    if (isRingMode) {
        ring(matrix, vector, result, rank, comm);
    } else {
        line(matrix, vector, result, rank, comm);
    }

    MPI_Comm_free(&comm);
    MPI_Finalize();
    return 0;
}