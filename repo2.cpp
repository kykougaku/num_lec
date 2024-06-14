#include <bits/stdc++.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <tuple>

using namespace std;
using namespace Eigen;
const int max_iter = 1000;

tuple<int, VectorXd> jacobi(MatrixXd A, MatrixXd b, double tol=1e-6){
    int n = A.rows();
    VectorXd x_new = VectorXd::Zero(n);
    MatrixXd D = MatrixXd::Zero(n,n);
    MatrixXd N = MatrixXd::Zero(n,n);
    VectorXd x = VectorXd::Zero(n);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i==j){
                D(i,j) = A(i,j);
            }
            else{
                N(i,j) = A(i,j);
            }
        }
    }
    MatrixXd inv_D = D.inverse();
    MatrixXd G = -inv_D*N;
    MatrixXd C = inv_D*b;
    double error = 1<<30;
    int epoch = 0;
    while(error>tol){
        epoch++;
        x_new = G*x + C;
        error = (x_new-x).norm();
        x = x_new;
        if(epoch>max_iter) break;
    }
    return make_tuple(epoch, x);
}

tuple<int, VectorXd> gaussseidel(MatrixXd A, MatrixXd b, double tol=1e-6){
    int n = A.rows();
    VectorXd x_new = VectorXd::Zero(n);
    MatrixXd D = MatrixXd::Zero(n,n);
    MatrixXd L = MatrixXd::Zero(n,n);
    MatrixXd U = MatrixXd::Zero(n,n);
    VectorXd x = VectorXd::Zero(n);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i==j){
                D(i,j) = A(i,j);
            }
            else if(i>j){
                L(i,j) = A(i,j);
            }
            else{
                U(i,j) = A(i,j);
            }
        }
    }
    MatrixXd DL = D+L;
    MatrixXd inv_DL = DL.inverse();
    MatrixXd H = -1*inv_DL*U;
    MatrixXd C = inv_DL*b;
    double error = 1<<30;
    int epoch = 0;
    while(error>tol){
        epoch++;
        x_new = H*x + C;
        error = (x_new-x).norm();
        x = x_new;
        if(epoch>max_iter) break;
    }
    return make_tuple(epoch, x);
}

int main(){
    Matrix4d A2;
    A2 << 7,2,4,4,
        6,8,1,4,
        2,-5,10,4,
        4,4,4,4;
    Matrix3d A1;
    A1 << 7,2,4,
        6,8,1,
        2,-5,10;
    Vector3d b1;
    b1 << 23,25,22;


    VectorXd b2;
    b2 << -78,25,32,10,-37;

    double tol = 1e-6;

    cout << "Jacobi Method in EQ1" << endl;
    cout << "tol: " << tol << endl;
    auto [epoch11, x11] = jacobi(A1,b1,tol);
    cout << "Epoch: " << epoch11 << endl;
    for(int i=0; i<3; i++){
        cout << x11(i) << endl;
    }
    cout << endl;

    cout << "Gauss-Seidel Method in EQ1" << endl;
    cout << "tol: " << tol << endl;
    auto [epoch12, x12] = gaussseidel(A1,b1,tol);
    cout << "Epoch: " << epoch12 << endl;
    for(int i=0; i<3; i++){
        cout << x12(i) << endl;
    }
    cout << endl;

    cout << "Jacobi Method in EQ2" << endl;
    cout << "tol: " << tol << endl;
    auto [epoch21, x21] = jacobi(A2,b2,tol);
    cout << "Epoch: " << epoch21 << endl;
    for(int i=0; i<5; i++){
        cout << x21(i) << endl;
    }
    cout << endl;

    cout << "Gauss-Seidel Method in EQ2" << endl;
    cout << "tol: " << tol << endl;
    auto [epoch22, x22] = gaussseidel(A2,b2,tol);
    cout << "Epoch: " << epoch22 << endl;
    for(int i=0; i<5; i++){
        cout << x22(i) << endl;
    }
    return 0;
}