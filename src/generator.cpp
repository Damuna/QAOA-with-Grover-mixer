//
// Created by Sören Wilkening on 07.06.23.
//
#include <iostream>
#include <random>

//#include <bits/stdc++.h>
using namespace std;

double eps;

int main()
{
#define int long long
//
////    cerr << "n=" << endl;
    int n;
//    cin >> n;
////    cerr << "capacity=" << endl;
    long long cap;
//    cin >> cap;
////    cerr << "classes=" << endl; // >= 2
    int classes;
//    cin >> classes;
//    classes--;
////    cerr << "fraction=" << endl; // 0<=frac<=1, in practice frac should be quite a lot smaller than 1
    double frac;
//    cin >> frac;
////    cerr << "eps=" << endl;
    double eps;
//    cin >> eps;
////    cerr << "small=" << endl;
    long long small;
//    cin >> small;
    FILE *file = fopen("generator_input.txt", "r");
    fscanf(file, "%lld", &n);
    fscanf(file, "%lld", &cap);
    fscanf(file, "%lld", &classes);
    classes--;
    fscanf(file, "%lf", &frac);
    fscanf(file, "%lf", &eps);
    fscanf(file, "%lld", &small);
    fclose(file);


    char filename[1023];
    char folder[1023];
    sprintf(filename, "../instances/n_%lld_c_%lld_g_%lld_f_%g_eps_%g_s_%lld/test.in", n, cap, classes, frac, eps, small);
    sprintf(folder, "mkdir -p ../instances/n_%lld_c_%lld_g_%lld_f_%g_eps_%g_s_%lld", n, cap, classes, frac, eps, small);
    system(folder);
    file = fopen(filename, "w");
    std::random_device device;
    std::mt19937 generator(device());
    std::uniform_int_distribution<int> distribution(1,small);
//    cout << n << endl;
    fprintf(file, "%lld\n", n);
    int amountSmall=n*frac;
    int am1=(n-amountSmall)/classes;
    double denominator=2.0;
    int amountCtr=0;
    for(int j=0; j<classes; j++)
    {
        for(int i=0; i<am1; i++)
        {
            int num1=distribution(generator);
            int num2=distribution(generator);
//            cout << amountCtr << " " << (int)((1/denominator+eps)*cap+num1) << " " << (int)((1/denominator+eps)*cap+num2) << endl;
            fprintf(file, "%lld %lld %lld\n", amountCtr, (int)((1/denominator+eps)*cap+num1), (int)((1/denominator+eps)*cap+num2));
            amountCtr++;
        }
        denominator*=2;
    }
    for(int i=amountCtr; i<n; i++)
    {
        int num1=distribution(generator);
        int num2=distribution(generator);
        fprintf(file, "%lld %lld %lld\n", i, num1, num2);
//        cout << i << " " << num1 << " " << num2 << endl;
    }
    fprintf(file, "%lld\n", cap);
//    cout << cap << endl;
    fclose(file);
    return 0;
}