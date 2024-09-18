#include <stdio.h>
#include <math.h>

#define SIZE 3
#define EPSILON 1e-10

// Function to perform the Jacobi rotation
void jacobiRotation(double A[SIZE][SIZE], double V[SIZE][SIZE], int p, int q) {
    double theta, t, c, s;
    double apq = A[p][q];
    double app = A[p][p];
    double aqq = A[q][q];

    if (apq != 0.0) {
        theta = 0.5 * atan2(2 * apq, aqq - app);
        t = tan(theta);
        c = 1.0 / sqrt(1 + t * t);
        s = c * t;
    } else {
        c = 1.0;
        s = 0.0;
    }

    // Update matrix A
    for (int i = 0; i < SIZE; i++) {
        double aip = A[i][p];
        double aiq = A[i][q];
        A[i][p] = A[p][i] = c * aip - s * aiq;
        A[i][q] = A[q][i] = c * aiq + s * aip;
    }
    A[p][p] = app * c * c - 2 * apq * c * s + aqq * s * s;
    A[q][q] = app * s * s + 2 * apq * c * s + aqq * c * c;
    A[p][q] = A[q][p] = 0.0;

    // Update eigenvector matrix V
    for (int i = 0; i < SIZE; i++) {
        double vip = V[i][p];
        double viq = V[i][q];
        V[i][p] = c * vip - s * viq;
        V[i][q] = c * viq + s * vip;
    }
}

// Function to perform the Jacobi method for finding eigenvalues and eigenvectors
void jacobiMethod(double A[SIZE][SIZE], double eigenvalues[SIZE], double eigenvectors[SIZE][SIZE]) {
    // Initialize eigenvectors to the identity matrix
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            eigenvectors[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    double maxOffDiagonal;
    int p, q;

    do {
        maxOffDiagonal = 0.0;
        p = 0;
        q = 1;

        // Find the largest off-diagonal element
        for (int i = 0; i < SIZE; i++) {
            for (int j = i + 1; j < SIZE; j++) {
                if (fabs(A[i][j]) > fabs(maxOffDiagonal)) {
                    maxOffDiagonal = A[i][j];
                    p = i;
                    q = j;
                }
            }
        }

        // If the largest off-diagonal element is larger than EPSILON, perform a Jacobi rotation
        if (fabs(maxOffDiagonal) > EPSILON) {
            jacobiRotation(A, eigenvectors, p, q);
        }

    } while (fabs(maxOffDiagonal) > EPSILON);

    // Extract eigenvalues
    for (int i = 0; i < SIZE; i++) {
        eigenvalues[i] = A[i][i];
    }
}

// Function to print a matrix
void printMatrix(double matrix[SIZE][SIZE]) {
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            printf("%10.4f ", matrix[i][j]);
        }
        printf("\n");
    }
}

// Main function
int main() {
    double A[SIZE][SIZE] = {
        {4, 1, 2},
        {1, 3, 1},
        {2, 1, 3}
    };

    double eigenvalues[SIZE];
    double eigenvectors[SIZE][SIZE];

    jacobiMethod(A, eigenvalues, eigenvectors);

    printf("Eigenvalues:\n");
    for (int i = 0; i < SIZE; i++) {
        printf("%10.4f\n", eigenvalues[i]);
    }

    printf("\nEigenvectors:\n");
    printMatrix(eigenvectors);

    return 0;
}
