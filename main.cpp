#include <iostream>
#include "mpi.h"


void task1(){
    int rank, size;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("Hello World from %d procces out of %d\n", rank, size);

    MPI_Finalize();

}

void task2(){
    const int n = 20;
    int a[n], max, maxp, rank, size, total_max;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (!rank)
    {
        for (int i = 0; i < n; i++) {
            a[i] = rand() % 50;
            if (i == 0) maxp = a[i];
            else if (maxp < a[i]) 
                maxp = a[i];
        }
    }

    MPI_Bcast(a, n, MPI_INT, 0, MPI_COMM_WORLD);
    max = a[0];

    int k = n / size;
    int i1 = k * (rank);
    int i2 = k * (rank + 1);
    if (rank == size - 1) {
        i2 = n;
    }

    for (int i = i1; i < i2; i++) {
        if (a[i] > max) max = a[i];
    }
    MPI_Reduce(&max, &total_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    if (!rank) {
        printf("\nTotal_max = %d, maxp = %d\n", total_max, maxp);
    }
    MPI_Finalize();

}

bool IsInCircle(double x, double y, double R) {
    return ((x * x + y * y) < R * R);
}

void task3(){
    const int NPoints = 10000;
    int NPointsInCircle = 0, correctPoints = 0;
    int rank, size;

    double PointX, PointY, pi;


    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int k = NPoints / size;
    int i1 = k * (rank);
    int i2 = k * (rank + 1);
    if (rank == size - 1)
        i2 = NPoints;
    

    for (int i = i1; i < i2; i++) {
        PointX = (rand() % 200 - 100) / 100.;
        PointY = (rand() % 200 - 100) / 100.;    
        
        if (IsInCircle(PointX, PointY, 1.)) {
            correctPoints++;
        }
    }
    MPI_Reduce(&correctPoints, &NPointsInCircle, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (!rank) {
        pi = NPointsInCircle * 4.0 / NPoints;
        printf("\nPi = %f\n", pi);
    }

    MPI_Finalize();
}

void task4(){
    const int n = 34;
    int a[n], count = 0, sum = 0, rest = n, sum_total = 0, count_total = 0, checks = 0,checkc = 0;
    int rank, size;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int *len = new int [size];
    int *ind = new int [size];
    int k = rest / size;
    len[0] = k;
    ind[0] = 0;

    for (int i = 1; i < size; i++) {
        rest -= k;
        k = rest / (size - i);
        len[i] = k;
        ind[i] = ind[i - 1] + len[i - 1];
    }

    if (rank == 0) {
        for (int i = 0; i < n; i++) {
            a[i] = rand() % 100 - 50;
            if (a[i] >= 0) { checkc++; checks += a[i]; }
        }
    }

    int sizeProcA = len[rank];
    int* procA = new int[sizeProcA];

    MPI_Scatterv(a, len, ind, MPI_INT, procA, sizeProcA, MPI_INT, 0, MPI_COMM_WORLD);

    for (int i = 0; i < sizeProcA; i++) {
        if (procA[i] >= 0) {
            count++;
            sum += procA[i];
        }
    }

    MPI_Reduce(&sum, &sum_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&count, &count_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (!rank){
        printf("\ncount_total = %d, sum_total = %d, avg = %f, avg_check = %f\n", count_total, sum_total, ((double)sum_total/count_total), ((double)checks/checkc));
    }
   
    MPI_Finalize();
}

void task5(){
    const int n = 34;
    int a[n], b[n], sum = 0, rest = n, sum_total = 0;
    int rank, size;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int* len = new int[size];
    int* ind = new int[size];
    int k = rest / size;
    len[0] = k;
    ind[0] = 0;

    for (int i = 1; i < size; i++) {
        rest -= k;
        k = rest / (size - i);
        len[i] = k;
        ind[i] = ind[i - 1] + len[i - 1];
    }

    if (!rank) {
        for (int i = 0; i < n; i++) {
            a[i] = rand() % 20;
            b[i] = rand() % 20;
        }
    }

    int sizeProc = len[rank];
    int* procA = new int[sizeProc];
    int* procB = new int[sizeProc];

    MPI_Scatterv(a, len, ind, MPI_INT, procA, sizeProc, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(b, len, ind, MPI_INT, procB, sizeProc, MPI_INT, 0, MPI_COMM_WORLD);

    for (int i = 0; i < sizeProc; i++) {
        sum += procA[i] * procB[i];
    }

    MPI_Reduce(&sum, &sum_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0)
        printf("Prod = %d\n", sum_total);
    

    MPI_Finalize();
}

void task6(){
   const int n = 10; 
   int *pMatrix = new int[n * n]; 
   int *pProcRows;
   int RowNum;
   int size, rank, tot_maxmin, tot_minmax;

   MPI_Init(NULL, NULL);

   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if (rank == 0) {
       for (int i = 0; i < n*n; ++i) {
           pMatrix[i] = rand() % 100 - 50;
       }
   }

   int* len = new int[size]; 
   int* ind = new int[size]; 
   int RestRows = n;

   RowNum = (n / size);
   len[0] = RowNum * n;
   ind[0] = 0;

   for (int i = 1; i < size; i++) {
       RestRows -= RowNum;
       RowNum = RestRows / (size - i);
       len[i] = RowNum * n;
       ind[i] = ind[i - 1] + len[i - 1];
   }

   pProcRows = new int[len[rank]];

   MPI_Scatterv(pMatrix, len, ind, MPI_INT, pProcRows, len[rank], MPI_INT, 0, MPI_COMM_WORLD);

   int partialMaxmin = -100000;
   int partialMinmax = 100000;

   int s = len[rank] / n;
   int min, max;

   for (int i = 0; i < s; i++) {

       min = 100000;
       max = -100000;

       for (int j = 0; j < n; j++)
       {
           if (pProcRows[i * n + j] < min)
               min = pProcRows[i * n + j];

           if (pProcRows[i * n + j] > max)
               max = pProcRows[i * n + j];
       }
       if (min > partialMaxmin) partialMaxmin = min;
       if (max < partialMinmax) partialMinmax = max;
   }


   MPI_Reduce(&partialMinmax, &tot_minmax, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
   MPI_Reduce(&partialMaxmin, &tot_maxmin, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

   if (rank == 0)
       printf("\nMaxmin = %d, minmax = %d\n", tot_maxmin, tot_minmax);
   
   MPI_Finalize();
}

void task10(){
    int size, rank;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int N = 10000;
    int* mass = new int[N];
    int* masr = new int[N];
    double st1, rt1, st2, rt2, st3, rt3, st4, rt4;
    MPI_Status st;
    if (rank == 0) {

        for (int i = 0; i < N; i++) {
            mass[i] = rand() % 20 - 7;
        }

        st1 = MPI_Wtime();
        MPI_Send(mass, N, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(masr, N, MPI_INT, 1, 0, MPI_COMM_WORLD, &st);
        rt1 = MPI_Wtime();
        printf("Time for Send = %f\n", rt1 - st1);

        st2 = MPI_Wtime();
        MPI_Ssend(mass, N, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(masr, N, MPI_INT, 1, 0, MPI_COMM_WORLD, &st);
        rt2 = MPI_Wtime();
        printf("Time for Ssend = %f\n", rt2 - st2);

        st3 = MPI_Wtime();
        MPI_Rsend(mass, N, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(masr, N, MPI_INT, 1, 0, MPI_COMM_WORLD, &st);
        rt3 = MPI_Wtime();
        printf("Time for Rsend = %f\n", rt3 - st3);

        int message_buffer_size = N * sizeof(int) + MPI_BSEND_OVERHEAD;
        double* message_buffer = (double*)malloc(message_buffer_size);
        st4 = MPI_Wtime();
        MPI_Buffer_attach(message_buffer, message_buffer_size);
        MPI_Bsend(mass, N, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(masr, N, MPI_INT, 1, 0, MPI_COMM_WORLD, &st);
        MPI_Buffer_detach(message_buffer, &message_buffer_size);
        free(message_buffer);
        rt4 = MPI_Wtime();
        printf("Time for Bsend = %f\n", rt4 - st4);
    }
    else {
        MPI_Recv(masr, N, MPI_INT, 0, 0, MPI_COMM_WORLD, &st);
        MPI_Send(mass, N, MPI_INT, 0, 0, MPI_COMM_WORLD);

        MPI_Recv(masr, N, MPI_INT, 0, 0, MPI_COMM_WORLD, &st);
        MPI_Ssend(mass, N, MPI_INT, 0, 0, MPI_COMM_WORLD);

        MPI_Recv(masr, N, MPI_INT, 0, 0, MPI_COMM_WORLD, &st);
        MPI_Rsend(mass, N, MPI_INT, 0, 0, MPI_COMM_WORLD);

        MPI_Recv(masr, N, MPI_INT, 0, 0, MPI_COMM_WORLD, &st);
        int message_buffer_size = N * sizeof(int) + MPI_BSEND_OVERHEAD;
        double* message_buffer = (double*)malloc(message_buffer_size);
        MPI_Buffer_attach(message_buffer, message_buffer_size);
        MPI_Bsend(mass, N, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Buffer_detach(message_buffer, &message_buffer_size);
    }

    MPI_Finalize();
}

void task11(){
    int  size, rank;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int N = 0;
    int sn = N;
    int rn;
    MPI_Status st;

    if (rank == 0) {
        MPI_Rsend(&sn,1,MPI_INT,1,0,MPI_COMM_WORLD);
    }

    for (int i = 0; i < size-1; i++) {
        if (i == rank) {
            MPI_Recv(&rn,1,MPI_INT,i-1,0,MPI_COMM_WORLD, &st);
            sn = rn + 1;
            MPI_Rsend(&sn, 1, MPI_INT, i+1, 0, MPI_COMM_WORLD);
        }        
    }

    if (rank == size - 1) {
        MPI_Recv(&rn, 1, MPI_INT, size-2, 0, MPI_COMM_WORLD, &st);
        sn = rn + 1;
        MPI_Rsend(&sn, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        MPI_Recv(&rn, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD, &st);
        N = rn+1;
        printf("Result: N = %d\n", N);
    }

    MPI_Finalize();
}

int main()
{
    // task1();
    // task2();
    // task3();
    // task4();
    // task5();
    // task6();
    // task10();
    // task11();

    return 0;
}

