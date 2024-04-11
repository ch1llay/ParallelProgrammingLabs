#include <stdio.h>
#include <mpi/mpi.h>
#include <limits.h>
#include <cstdlib>
#include <iostream>
#include <csignal>

#pragma warning(disable : 4996)

using namespace std;

#define SIZE 5
#define DIMS_DG 2
#define DIMS_SFG 1

int receive(int source);

void send(int destination, int array);

void sendTo0Process(const int sendData, int *recvbuf, int rank) {
    if (rank == 0) {
        recvbuf[0] = sendData;
        for (int i = 1; i < SIZE; i++) {
            recvbuf[i] = receive(i);
        }

    } else {
        send(0, sendData);
    }

}


int *sortByDependencyGraph(int inputArray[SIZE], int rank, MPI_Comm comm) {
    int min, max;
    int coords[2];

    MPI_Cart_coords(comm, rank, DIMS_DG, coords);

    if (coords[0] <= coords[1]) {
        int ySource, yDest, xSource, xDest;

        MPI_Cart_shift(comm, 0, 1, &ySource, &yDest);
        MPI_Cart_shift(comm, 1, 1, &xSource, &xDest);

        min = coords[0] == coords[1] ? INT_MAX : receive(xSource);

        if (coords[0] == 0) {
            if (coords[1] == 0) {
                max = inputArray[0];
                for (int i = 1; i < SIZE; ++i)
                    send(i, inputArray[i]);
            } else {
                max = receive(0);
            }
        } else {
            max = receive(ySource);
        }

        if (min > max) {
            int temp = min;
            min = max;
            max = temp;
        }

        if (xDest < 0) {
            send(0, min);
        } else {
            send(xDest, min);
        }

        if (coords[0] != coords[1])
            send(yDest, max);

        if (rank == 0) {
            int *result = (int *) malloc(sizeof(int) * SIZE);
            int target = 0;

            for (int i = 0; i < SIZE; ++i) {
                if (i == 0)
                    target += SIZE - 1;
                else
                    target += SIZE;

                int currentItem = receive(target);
                result[i] = currentItem;
            }
            return result;
        }
    }

    return NULL;
}

int *sortBySignalsGraph(int inputArray[SIZE], int rank, MPI_Comm comm) {
    int min = INT_MAX;
    int source, dest;

    MPI_Cart_shift(comm, 0, 1, &source, &dest);


    for (int i = 0; i < SIZE; ++i) {
        int max = rank == 0 ? inputArray[i] : receive(source);

        if (min > max) {
            int temp = min;
            min = max;
            max = temp;
        }

        if (dest > 0) {
            send(dest, max);
        }
    }

    int *result = rank == 0 ? (int *) malloc(sizeof(int) * SIZE) : NULL;
    sendTo0Process(min, result, rank);
    return result;
}

int main(int argc, char **argv) {

    int rank;
    int size;
    MPI_Comm comm;

    MPI_Init(&argc, &argv);


//    int i = 0;
//    while (!i)
//        sleep(5);


    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    bool isGraphDependency = true;
    if (size <= SIZE) {
        isGraphDependency = false;
    }


    int inputArray[SIZE] = {2, 5, 3, 1, 4};
    int *result;
    if (rank == 0) {
        cout << "Исходный массив:" << endl;;
        for (int i = 0; i < SIZE; ++i) {
            cout << inputArray[i] << "\t";
        }
        cout << endl;
    }

    if (isGraphDependency) {
        int dims[DIMS_DG];
        int periods[DIMS_DG];

        for (int i = 0; i < DIMS_DG; i++) {
            dims[i] = 0;
            periods[i] = 0;
        }
        MPI_Dims_create(size, DIMS_DG, dims);

        MPI_Cart_create(MPI_COMM_WORLD, DIMS_DG, dims, periods, 1, &comm);

        if (rank == 0) {
            cout << "Граф зависимостей" << endl;
        }

        result = sortByDependencyGraph(inputArray, rank, comm);

        if (rank == 0) {
            cout << "Результат:" << endl;
            for (int i = 0; i < SIZE; ++i) {
                cout << result[i] << "\t";
            }
            cout << endl;;
        }
    } else {
        int dims[DIMS_SFG];
        int periods[DIMS_SFG];


        for (int i = 0; i < DIMS_SFG; i++) {
            dims[i] = 0;
            periods[i] = 0;
        }
        MPI_Dims_create(size, DIMS_SFG, dims);
        MPI_Cart_create(MPI_COMM_WORLD, DIMS_SFG, dims, periods, 1, &comm);

        if (rank == 0) {
            cout << "Граф потока сигналов" << endl;
        }

        result = sortBySignalsGraph(inputArray, rank, comm);

        if (rank == 0) {
            cout << "Результат:" << endl;;
            for (int i = 0; i < SIZE; ++i) {
                cout << result[i] << "\t";
            }
            cout << endl;
        }
    }
    free(result);

    MPI_Finalize();
    return 0;
}

int receive(int source) {
    int value;
    MPI_Status status;
    MPI_Recv(&value, 1, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    return value;
}

void send(int destination, int X) {
    MPI_Send(&X, 1, MPI_INT, destination, 0, MPI_COMM_WORLD);
}