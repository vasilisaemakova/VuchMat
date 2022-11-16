#include <iostream>
#include <cmath>

using namespace std;

double fieldStrength(double E0, double radius) {
    return E0 / radius;
    //return E0 * radius;
}

void systemFunction(double *X_data, double *K_data, double *FX_data, double (*fieldStrength)(double, double)) {
    FX_data[0] = X_data[3];
    FX_data[1] = X_data[4];
    FX_data[2] = X_data[5];
    double radius = sqrt(X_data[0] * X_data[0] + X_data[2] * X_data[2]);
    double E = fieldStrength(K_data[3], radius);
    FX_data[3] = K_data[0] * X_data[5] * K_data[1] / K_data[2] + K_data[1] * E * X_data[0] / (K_data[2] * radius);
    FX_data[4] = 0;
    FX_data[5] = -K_data[0] * X_data[3] * K_data[1] / K_data[2] + K_data[1] * E * X_data[2] / (K_data[2] * radius);
}

void Euler(double *X_data, double *K_data, double *FX_data,
           void (*systemFunction)(double *, double *, double *, double (*fieldStrength)(double, double)),
           double (*fieldStrength)(double, double), double hop, int dim) {
    systemFunction(X_data, K_data, FX_data, fieldStrength);
    for (int i = 0; i < dim; i++) {
        X_data[i] += hop * FX_data[i];
    }
}

void RK4(double *X_data, double *K_data, double *FX_data,
         void (*systemFunction)(double *, double *, double *, double (*fieldStrength)(double, double)),
         double (*fieldStrength)(double, double), double hop, int dim) {
    double *tmp_X_data = new double[dim];

    systemFunction(X_data, K_data, FX_data, fieldStrength);

    for (int i = 0; i < dim; i++) {
        tmp_X_data[i] = X_data[i] + 0.5 * hop * FX_data[0 * 6 + i];
    }
    systemFunction(tmp_X_data, K_data, FX_data + 1 * 6, fieldStrength);

    for (int i = 0; i < dim; i++) {
        tmp_X_data[i] = X_data[i] + 0.5 * hop * FX_data[1 * 6 + i];
    }
    systemFunction(tmp_X_data, K_data, FX_data + 2 * 6, fieldStrength);

    for (int i = 0; i < dim; i++) {
        tmp_X_data[i] = X_data[i] + 1 * hop * FX_data[2 * 6 + i];
    }
    systemFunction(tmp_X_data, K_data, FX_data + 3 * 6, fieldStrength);

    for (int i = 0; i < dim; i++) {
        X_data[i] += hop * FX_data[0 * 6 + i] / 6 + hop * FX_data[1 * 6 + i] / 3 +
                     hop * FX_data[2 * 6 + i] / 3 + hop * FX_data[3 * 6 + i] / 6;
    }

}

int main() {
    double *X_data = new double[6];
    X_data[0] = 1; //x
    X_data[1] = 0; //y
    X_data[2] = 1; //z
    X_data[3] = 1; //Vx
    X_data[4] = 0; //Vy
    X_data[5] = 1; //Vz

    double *K_data = new double[4];
    K_data[0] = 0.3; //B
    K_data[1] = 1; //q
    K_data[2] = 1; //m
    K_data[3] = 1; //E0

    double *FX_data = new double[6 * 4];
    double hop = 0.01;
    double observationTime = 200;
    //double currentTime = 0;
    for (double t = 0; t < observationTime; t += hop) {
        RK4(X_data, K_data, FX_data, systemFunction, fieldStrength, hop, 6);
        //Euler(X_data, K_data, FX_data, systemFunction, fieldStrength, hop, 6);

        cout << X_data[0] << ' ' << X_data[1] << ' ' << X_data[2] << endl;
    }

}
