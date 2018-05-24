#include<iostream>
#include<eigen3/Eigen/Core>


int main()
{
    int i, j;
    Eigen::VectorXd v;
    Eigen::MatrixXd m;

    v.resize(15);
    for(i=0; i<15; i++)
        v(i) = i;

    Eigen::VectorXd a(5);
    for(i=0; i<5; i++)
        a[i] = 0;
    v.segment(3,5) = a;


    std::cout << v << "\n" << std::endl;


    
    return 0;
}