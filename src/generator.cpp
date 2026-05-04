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
    int n;
    long long cap;
    int classes;
    double frac;
    double eps;
    long long small;

    FILE *file = fopen("executables/generator_input.txt", "r");
    fscanf(file, "%lld", &n);
    fscanf(file, "%lld", &cap);
    fscanf(file, "%lld", &classes);
    classes--;
    fscanf(file, "%lf", &frac);
    fscanf(file, "%lf", &eps);
    fscanf(file, "%lld", &small);
    fclose(file);
    
    std::random_device device;
    std::mt19937 generator(device());
    std::uniform_int_distribution<int> distribution(1, small);
    
    
    for (int i = 0; i < 5; ++i) {
        char filename[1023];
        char folder[1023];
#if defined(_WIN32) || defined(_WIN32)
        sprintf(filename, "..\\instances\\n_%lld_c_%lld_g_%lld_f_%g_eps_%g_s_%lld\\test.in", n, cap, classes, frac, eps, small);
        sprintf(folder, "mkdir ..\\instances\\n_%lld_c_%lld_g_%lld_f_%g_eps_%g_s_%lld 2>nul", n, cap, classes, frac, eps, small);
#else
        sprintf(filename,
                "new_instances/n_%lld_c_%lld_g_%lld_f_%g_eps_%g_s_%lld/test_%d.in",
                n,
                cap,
                classes,
                frac,
                eps,
                small,
                i);
        sprintf(folder,
                "mkdir -p new_instances/n_%lld_c_%lld_g_%lld_f_%g_eps_%g_s_%lld",
                n,
                cap,
                classes,
                frac,
                eps,
                small,
                i);
#endif
        system(folder);
        file = fopen(filename, "w");
//    cout << n << endl;
        fprintf(file, "%lld\n", n);
        int amountSmall = n * frac;
        int am1 = (n - amountSmall) / classes;
        double denominator = 2.0;
        int amountCtr = 0;
        for (int j = 0; j < classes; j++) {
            for (int i = 0; i < am1; i++) {
                int num1 = distribution(generator);
                int num2 = distribution(generator);
//            cout << amountCtr << " " << (int)((1/denominator+eps)*cap+num1) << " " << (int)((1/denominator+eps)*cap+num2) << endl;
                fprintf(file,
                        "%lld %lld %lld\n",
                        amountCtr,
                        (int) ((1 / denominator + eps) * cap + num1),
                        (int) ((1 / denominator + eps) * cap + num2));
                amountCtr++;
            }
            denominator *= 2;
        }
        for (int i = amountCtr; i < n; i++) {
            int num1 = distribution(generator);
            int num2 = distribution(generator);
            fprintf(file, "%lld %lld %lld\n", i, num1, num2);
//        cout << i << " " << num1 << " " << num2 << endl;
        }
        fprintf(file, "%lld\n", cap);
//    cout << cap << endl;
        fclose(file);
//        n += 500;
    }
    return 0;
}