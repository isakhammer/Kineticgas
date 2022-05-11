#pragma once
#include "KineticGas.h"
#include "global_params.h"
#include "math.h"
#include "cmath"
#include <cstdio>

class Fac{
    public:
    int val;

    Fac(int v);
    long long eval();
};

class Product{
    public:
    int isize;
    int dsize;
    int ilist[1000];
    double dlist[1000];

    Product(int const& i);
    Product(double const& d);
    Product(Fac const& f);
    Product(Product const& p);

    double eval();

    Product operator*=(const int& rhs);
    Product operator*=(const double& rhs);
    Product operator*(const Product& rhs);
    Product operator*=(const Product& rhs);
};

class Frac{
    public:
    Product numerator;
    Product denominator;

    Frac(Product num, Product den);
    Frac(Product num);
    Frac(Frac const& f);

    double eval();

    Frac operator*(const Product& rhs);
    Frac operator*(const double& rhs);
    Frac operator*(const int& rhs);
    Frac operator*(const Frac& rhs);
    Frac operator/(const Frac& rhs);
    Frac operator/(const Product& rhs);
};

Product operator*(const Fac& f1, const Fac& f2);
Product operator*(const Fac& lhs, const double& rhs);
Product operator*(const double& lhs, Fac& rhs);
Product operator*(const int& lhs, Product& rhs);
Product operator*(const double& lhs, Product& rhs);
Frac operator/(const Product& lhs, const Product& rhs);
Frac operator/(const double& lhs, const Product& rhs);
Frac operator/(const Product& lhs, const double& rhs);

double operator+(Frac& lhs, Frac& rhs);
double operator+=(double& lhs, const Frac& rhs);

Product ipow(int base, int expo);
int factorial_tests();