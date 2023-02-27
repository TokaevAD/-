//
//  main.cpp
//  parallelki2
//
//  Created by Артем on 21.02.2023.
//

#include <iostream>
#include <omp.h>
#include <cstdlib>

using namespace std;

//int main() {
//    #pragma omp parallel num_threads(8)
//    {
//        int id = omp_get_thread_num();
//        int numThreads = omp_get_num_threads();
//
//        printf("Thread %d of %d Says: Hello World!\n",id, numThreads);
//    }
//
//    return 0;
//}



//int main()
//{
//    const int n = 16000;
//    int a[n], b[n];
//    double start_time, end_time;
//
//    // инициализируем массив a
//    #pragma omp parallel for
//    for (int i = 0; i < n; ++i)
//    {
//        a[i] = i;
//    }
//
//    // вычисление массива b с распределением static
//        start_time = omp_get_wtime();
//        #pragma omp parallel for schedule(static)
//        for (int i = 1; i < n-1; i++) {
//            b[i] = (a[i-1] + a[i] + a[i+1])/3.0;
//        }
//        end_time = omp_get_wtime();
//        printf("Static time: %f\n", -(start_time-end_time));
//
//    // вычисление массива b с распределением dynamic
//        start_time = omp_get_wtime();
//        #pragma omp parallel for schedule(dynamic)
//        for (int i = 1; i < n-1; i++) {
//            b[i] = (a[i-1] + a[i] + a[i+1])/3.0;
//        }
//        end_time = omp_get_wtime();
//        printf("Dynamic time: %f\n", -(start_time-end_time));
//
//    // вычисление массива b с распределением auto
//        start_time = omp_get_wtime();
//        #pragma omp parallel for schedule(auto)
//        for (int i = 1; i < n-1; i++) {
//            b[i] = (a[i-1] + a[i] + a[i+1])/3.0;
//        }
//        end_time = omp_get_wtime();
//        printf("Auto time: %f\n", -(start_time-end_time));
//    // выводим результат
////    for (int i = 0; i < n; ++i)
////    {
////        //std::cout << "b[" << i << "] = " << b[i] << std::endl;
////        printf("b[%d] = %d\n",i,b[i]);
////    }
//    printf("The winner is ...");
//
//    return 0;
//}





//const int N = 1000; // размер матриц
//
//void generateMatrix(int **matrix) {
//    //srand(time(nullptr));
//    for (int i = 0; i < N; i++) {
//        for (int j = 0; j < N; j++) {
//            matrix[i][j] = rand() % 10;
//        }
//    }
//}
//
//void sequential_multiply(double **A, double **B, double **C) {
//    for (int i = 0; i < N; i++) {
//        for (int j = 0; j < N; j++) {
//            C[i][j] = 0;
//            for (int k = 0; k < N; k++) {
//                C[i][j] += A[i][k] * B[k][j];
//            }
//        }
//    }
//}
//
//void parallel_multiply(double **A, double **B, double **C, int num_threads) {
//    omp_set_num_threads(num_threads);
//    #pragma omp parallel for
//    for (int i = 0; i < N; i++) {
//        for (int j = 0; j < N; j++) {
//            C[i][j] = 0;
//            for (int k = 0; k < N; k++) {
//                C[i][j] += A[i][k] * B[k][j];
//            }
//        }
//    }
//}
//
//int main() {
//    double **A, **B, **C;
//    A = new double*[N];
//    B = new double*[N];
//    C = new double*[N];
//
//    for (int i = 0; i < N; i++) {
//        A[i] = new double[N];
//        B[i] = new double[N];
//        C[i] = new double[N];
//        for (int j = 0; j < N; j++) {
//            A[i][j] = rand() % 10;
//            B[i][j] = rand() % 10;
//        }
//    }
//
//    double start_time, end_time;
//
//    // sequential multiply
//    start_time = omp_get_wtime();
//    sequential_multiply(A, B, C);
//    end_time = omp_get_wtime();
//
//    printf("Sequential execution time: %f s\n", end_time - start_time);
//    // parallel multiply with 2 threads
//    start_time = omp_get_wtime();
//    parallel_multiply(A, B, C, 2);
//    end_time = omp_get_wtime();
//
//    printf("Parallel execution time with 2 threads: %f s\n", end_time - start_time);
//
//    // parallel multiply with 4 threads
//    start_time = omp_get_wtime();
//    parallel_multiply(A, B, C, 4);
//    end_time = omp_get_wtime();
//    printf("Parallel execution time with 4 threads: %f s\n", end_time - start_time);
//
//    // parallel multiply with 8 threads
//    start_time = omp_get_wtime();
//    parallel_multiply(A, B, C, 8);
//    end_time = omp_get_wtime();
//    printf("Parallel execution time with 8 threads: %f s\n", end_time - start_time);
//
//    // free memory
//    for (int i = 0; i < N; i++) {
//        delete[] A[i];
//        delete[] B[i];
//        delete[] C[i];
//    }
//    delete[] A;
//    delete[] B;
//    delete[] C;
//
//    return 0;
//}



//int main() {
//    #pragma omp parallel num_threads(8)
//    {
//        int thread_id = omp_get_thread_num();
//        int num_threads = omp_get_num_threads();
//
//        #pragma omp ordered
//        {
//
//            //не работает
//                   printf("Hello wolrd! Thread %d out of %d \n", num_threads - thread_id, num_threads);
//        }
//    }
//
//    return 0;
//}



#include <vector>

int main() {
    int num_threads = 8;
//    std::vector<int> thread_ids(num_threads);
//
//    #pragma omp parallel num_threads(num_threads)
//    {
//        int thread_id = omp_get_thread_num();
//        thread_ids[thread_id] = thread_id;
//    }
//
//    for (int i = num_threads - 1; i >= 0; --i) {
//        std::cout << "Hello world! Thread " << thread_ids[i] << " out of " << num_threads << std::endl;
//    }

    printf("\n");
    #pragma omp parallel
    {
        #pragma omp critical
        {
        int thread_id = omp_get_thread_num();
        int num_threads = omp_get_num_threads();
        printf("Hello World from thread %d of %d\n", thread_id, num_threads);
        }
    }


    
    return 0;
}


