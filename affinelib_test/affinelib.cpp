//   test file for Library for affine transformation
//   by Shizuo KAJI,     Nov. 2013

//#define  EIGEN_USE_MKL_ALL

#include "affinelib.h"
#include <iostream>
#include <ctime>
#include <vector>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Eigen>
#include <Eigen/Geometry>

using namespace Eigen;
using namespace std;
using namespace AffineLib;

int main(void){
    int count=1000000;
    int num=5;
    std::vector<Matrix4d> Aff(count), logSE(count), SE(count), res4(count), prevLogSE(count);
    std::vector<Matrix3d> GL(count), U(count),R(count),S(count),logS(count),logR(count);
    std::vector<Matrix3d> resS(count), resR(count), res3(count), rLogS(count);
    std::vector<Vector3d> s(count), L(count);
    MatrixXf resXS,resXR;

    // prepare random matrices
    for(int j=0;j<count;j++){
        Matrix3d A=Matrix3d::Random();
        logS[j] = A+A.transpose();
        logR[j] = A-A.transpose();
        S[j] = expSym(logS[j]);
        R[j] = expSO(logR[j]);
        GL[j] = S[j] * R[j];
        L[j] = Vector3d::Random();
        SE[j]= pad(R[j], L[j]);
        Aff[j]=pad(GL[j],L[j]);
    }
    for(int j=0;j<count;j++){
        prevLogSE[j] = Matrix4d::Zero().eval();
    }
    
    PRINT_MAT(GL[0]);
    PRINT_MAT(S[0]);
    PRINT_MAT(R[0]);
    
    double t, min_t=HUGE_VAL, max_t=0.0, sum_t=0.0;
    for(int i=0;i<num;i++){
        Matrix3d result;
        clock_t clock_start=clock();
        for(int j=0;j<count;j++){
//            res4[j] = X[j].log();
//            res4[j]=logSEc(SE[j], prevTheta[j], prevN[j]);
//            res4[j] = SE[j].log();
//            res4[j] = expSE(logSE[j]);
//            res4[j] = logSE[j].exp();
//            res3[j] = expSym(logS[j]);
//            res3[j] = logS[j].exp();
//            resS[j] = S[j].log();
//            resS[j] = logSymD(expSymD(logS[j]));
//            Vector3f lambda;
//            resS[j] = logSym(S[j], lambda);
            
            // parametrisation
//            polar(X[j].block(0,0,3,3),U[j],s[j],resR[j]);
//            res3[j] = logDiag(U[j],s[j]);
           // parametrisation by spectral decomp
//            parametriseGL(X[j].block(0,0,3,3),res3[j],resR[j]);
//            res4[j]=logSEc(affine(resR[j],L[j]), prevTheta[j], prevN[j]);

//            res4[j] = affine(expSym(logS[j]))*expSE(logSE[0]);
            // polar decomposition
            polarByParam(GL[j],res3[j],resR[j]);
//              polarHigham(GL[j],res3[j],resR[j]);
//            polarBySVD(GL[j],U[j],s[j],resR[j]);
//              polarDiag(GL[j],U[j],s[j],resR[j]);
//            res3[j] = U[j] * s[j].asDiagonal() * U[j].transpose();
//            MatrixXd Y=MatrixXf::Random(4,4);
//            PRINT_MAT(Y);
//            polarN(Y,resXS,resXR);
//            PRINT_MAT(resXS*resXR);
        }
        t=(double)(clock()- clock_start)/CLOCKS_PER_SEC;
        min_t=std::min(t, min_t);
        max_t=std::max(t, max_t);
        sum_t += t;
        cout<< "Elapsed time: "<< t <<endl;
    }
    
    //	    PRINT_MAT(affine(logS[0].exp(),Vector3f(0,0,0))*logSE[0].exp());
    PRINT_MAT(res3[0]);
    PRINT_MAT(resR[0]);

    
	cout << "num of loops:" << count << endl;
    cout<< "Min Time in " << num << " trials: "<< min_t << "s" <<endl;
    cout<< "Max Time in " << num << " trials: "<< max_t << "s" <<endl;
    cout<< "Mean Time in " << num << " trials: "<< sum_t/num << "s" <<endl;

    /**
    MatrixXf A = MatrixXf::Random(1000,1000);
    clock_t clock_start=clock();
    MatrixXf B(A.inverse());
    cout<< "Inverse: " << (double)(clock()- clock_start)/CLOCKS_PER_SEC <<"s" <<endl;

    std::vector<float> w(3);
    std::vector<Matrix3f> mm(3);
    w[0] = 0.33;
    w[1] = 0.33;
    w[2] = 0.33;
    mm[0] = R;
    mm[1] = R;
    mm[2] = R;
    PRINT_MAT(frechetSO(mm, w)-R);
**/    

//    cin.get();
    return(0);
}