#include <iostream>
#include <fstream>
#include<cmath>
#include <iomanip>
#include<string>
double pi = 3.14159265358979323;
double L, Tmax, tau, h, eps, a, u0, t0,ct, aprox;
int indU, n, ind,indD,st;
double U0(double x)
{
    if (indU == 0)
    {
        return std::sin(pi * x);
    }
    if (indU == 1)
    {
        return x*(1-x);
    }
    if (indU == 2)
    {
        return (1-x*x)*std::cos(pi*x);
    }
    if (indU == 3)
    {
        return 0.5*(1+x)*(1+x);
    }
    if (indU == 4)
    {
        if (x > -1. / 3.)
        {
            if (x < 1. / 3.)
            {
                return 1;
            }
        }
        return 0;
    }
    if (indU == 5)
    {
        return 0;
    }
    if (indU == 6)
    {
        return 0;
    }
}
double UT(double x)
{
    if (indU == 0)
    {
        return 0;
    }
    if (indU == 1)
    {
        return 0;
    }
    if (indU == 2)
    {
        return 2*x+0.6;
    }
    if (indU == 3)
    {
        return (x+0.5)*std::cos(pi*x);
    }
    if (indU == 4)
    {
        if (x == -1 / 3)
        {
            return 1 / (h*h* h);
        }
        if (x == 1 / 3)
        {
            return 1 / (h*h*h);
        }
        return 0;
    }
    if (indU == 5)
    {
        if (x >= -1. / 2.)
        {
            if (x <= 1. / 2.)
            {
                return 1-2*std::fabs(x);
            }
        }
        return 0;
    }
    if (indU == 6)
    {
        return 0;
    }
}
double U0xx(double x)
{
    if (indU == 0)
    {
        return -pi*pi*std::sin(pi * x);
    }
    if (indU == 1)
    {
        return -2;
    }
    if (indU == 2)
    {
        return -2*std::cos(pi*x)-pi*pi*(1-x*x)*std::cos(pi*x)+4*pi*x*std::sin(pi*x);
    }
    if (indU == 3)
    {
        return 1;
    }
    if (indU == 4)
    {
        
        return 0;
    }
    if (indU == 5)
    {
        return 0;
    }
    if (indU == 6)
    {
        return 0;
    }
}
double AU0xx(double yl, double y, double yn)
{
    return (yn - 2 * y + yl) / (h * h);
}
double UL(double t)
{
    if (indU == 0)
    {
        return 0;
    }
    if (indU == 1)
    {
        return 0;
    }
    if (indU == 2)
    {
        return 1 + 0.4 * t;
    }
    if (indU == 3)
    {
        return 0.5;
    }
    if (indU ==4)
    {
        return 0;
    }
    if (indU == 5)
    {
        return 0;
    }
    if (indU == 6)
    {
        return std::sin(t);
    }
}
double UR(double t)
{
    if (indU == 0)
    {
        return 0;
    }
    if (indU == 1)
    {
        return 0;
    }
    if (indU == 2)
    {
        return 0;
    }
    if (indU == 3)
    {
        return 2-3*t;
    }
    if (indU == 4)
    {
        return 0;
    }
    if (indU == 5)
    {
        return 0;
    }
    if (indU == 6)
    {
        return 0;
    }
}
void MStep(double* tek, double* x, double* res0, double* res1,double* res2, int n)
{
    for (int i = 1; i < n;i++)
    {
        res2[i] = 2 * res1[i] - res0[i] + (a * a * tau * tau / (h * h)) * (res1[i + 1] - 2 * res1[i] + res1[i - 1]);
    }
    res2[0] = UL(tek[2]);
    res2[n] = UR(tek[2]);
}
void test1()
{
    std::ofstream ans;
    ans.open(std::to_string(indD) + "k" + std::to_string(ct) + "test1.txt");
    ans << std::setprecision(15);
    h = 0.05;
    tau = ct*0.01;
    a = 5;
    t0 = 0;
    u0 = 0;
    L = 1;
    Tmax = 10;
    indU = 0;
    n = int(L / h);
    double* x;
    double* res0;
    double* res1;
    double* res2;
    double* tek;
    res0 = new double[n + 1];
    res1 = new double[n + 1];
    res2 = new double[n + 1];
    x = new double[n + 1];
    tek = new double[3];
    tek[0] = t0;
    tek[1] = tek[0] + tau;
    tek[2] = tek[1] + tau;
    for (int i = 0; i < n + 1; i++)
    {
        x[i] = u0 + i * h;
        res0[i] = U0(x[i]);
    }
    if (indD == 0)
    {
        for (int i = 0; i < n + 1; i++)
        {
            res1[i] = res0[i] + tau * UT(x[i]) + (a * a * tau * tau / 2) * U0xx(x[i]);
        }
    }
    res1[0] = UL(tek[1]);
    res1[n]= UR(tek[1]);
    if (indD == 1)
    {
        for (int i = 1; i < n; i++)
        {
            res1[i] = res0[i] + tau * UT(x[i]) + (a * a * tau * tau / 2) * AU0xx(res0[i-1],res0[i],res0[i+1]);
        }
    }
    ans << h << std::endl;
    ans << tek[0] << ' ';
    for (int i = 0; i < n + 1; i++)
    {
        ans << res0[i] << ' ';
    }
    ans << std::endl;
    ans << tek[1] << ' ';
    for (int i = 0; i < n + 1; i++)
    {
        ans << res1[i] << ' ';
    }
    ans << std::endl;
    while (tek[2] < (Tmax+tau))
    {
        MStep(tek, x, res0, res1, res2, n);
        ans << tek[2] << ' ';
        for (int i = 0; i < n + 1; i++)
        {
            ans << res2[i] << ' ';
        }
        ans << std::endl;
        for (int i = 0; i < n + 1; i++)
        {
            res0[i] = res1[i];
            res1[i] = res2[i];
        }
        tek[0] = tek[0] + tau;
        tek[1] = tek[1] + tau;
        tek[2] = tek[2] + tau;
    }
    ans.close();
    delete[] x;
    delete[] tek;
    delete[] res1;
    delete[] res2;
    delete[] res0;
}
void test2()
{
    std::ofstream ans;
    ans.open(std::to_string(indD) + "k" + std::to_string(ct) + "test2.txt");
    ans << std::setprecision(15);
    h = 0.05;
    tau = ct * 0.01;
    a = 5;
    t0 = 0;
    u0 = 0;
    L = 1;
    Tmax = 10;
    indU = 1;
    n = int(L / h);
    double* x;
    double* res0;
    double* res1;
    double* res2;
    double* tek;
    res0 = new double[n + 1];
    res1 = new double[n + 1];
    res2 = new double[n + 1];
    x = new double[n + 1];
    tek = new double[3];
    tek[0] = t0;
    tek[1] = tek[0] + tau;
    tek[2] = tek[1] + tau;
    for (int i = 0; i < n + 1; i++)
    {
        x[i] = u0 + i * h;
        res0[i] = U0(x[i]);
    }
    if (indD == 0)
    {
        for (int i = 0; i < n + 1; i++)
        {
            res1[i] = res0[i] + tau * UT(x[i]) + (a * a * tau * tau / 2) * U0xx(x[i]);
        }
    }
    res1[0] = UL(tek[1]);
    res1[n] = UR(tek[1]);
    if (indD == 1)
    {
        for (int i = 1; i < n; i++)
        {
            res1[i] = res0[i] + tau * UT(x[i]) + (a * a * tau * tau / 2) * AU0xx(res0[i - 1], res0[i], res0[i + 1]);
        }
    }
    ans << h << std::endl;
    ans << tek[0] << ' ';
    for (int i = 0; i < n + 1; i++)
    {
        ans << res0[i] << ' ';
    }
    ans << std::endl;
    ans << tek[1] << ' ';
    for (int i = 0; i < n + 1; i++)
    {
        ans << res1[i] << ' ';
    }
    ans << std::endl;
    while (tek[2] < (Tmax + tau))
    {
        MStep(tek, x, res0, res1, res2, n);
        ans << tek[2] << ' ';
        for (int i = 0; i < n + 1; i++)
        {
            ans << res2[i] << ' ';
        }
        ans << std::endl;
        for (int i = 0; i < n + 1; i++)
        {
            res0[i] = res1[i];
            res1[i] = res2[i];
        }
        tek[0] = tek[0] + tau;
        tek[1] = tek[1] + tau;
        tek[2] = tek[2] + tau;
    }
    ans.close();
    delete[] x;
    delete[] tek;
    delete[] res1;
    delete[] res2;
    delete[] res0;
}
void test16()
{
    std::ofstream ans;
    ans.open(std::to_string(indD) + "k" + std::to_string(ct) + "v16.txt");
    ans << std::setprecision(15);
    h = 0.05;
    tau = ct * 0.01;
    a = 5;
    t0 = 0;
    u0 = 0;
    L = 1;
    Tmax = 10;
    indU = 2;
    n = int(L / h);
    double* x;
    double* res0;
    double* res1;
    double* res2;
    double* tek;
    res0 = new double[n + 1];
    res1 = new double[n + 1];
    res2 = new double[n + 1];
    x = new double[n + 1];
    tek = new double[3];
    tek[0] = t0;
    tek[1] = tek[0] + tau;
    tek[2] = tek[1] + tau;
    for (int i = 0; i < n + 1; i++)
    {
        x[i] = u0 + i * h;
        res0[i] = U0(x[i]);
    }
    if (indD == 0)
    {
        for (int i = 0; i < n + 1; i++)
        {
            res1[i] = res0[i] + tau * UT(x[i]) + (a * a * tau * tau / 2) * U0xx(x[i]);
        }
    }
    res1[0] = UL(tek[1]);
    res1[n] = UR(tek[1]);
    if (indD == 1)
    {
        for (int i = 1; i < n; i++)
        {
            res1[i] = res0[i] + tau * UT(x[i]) + (a * a * tau * tau / 2) * AU0xx(res0[i - 1], res0[i], res0[i + 1]);
        }
    }
    ans << h << std::endl;
    ans << tek[0] << ' ';
    for (int i = 0; i < n + 1; i++)
    {
        ans << res0[i] << ' ';
    }
    ans << std::endl;
    ans << tek[1] << ' ';
    for (int i = 0; i < n + 1; i++)
    {
        ans << res1[i] << ' ';
    }
    ans << std::endl;
    while (tek[2] < (Tmax + tau))
    {
        MStep(tek, x, res0, res1, res2, n);
        ans << tek[2] << ' ';
        for (int i = 0; i < n + 1; i++)
        {
            ans << res2[i] << ' ';
        }
        ans << std::endl;
        for (int i = 0; i < n + 1; i++)
        {
            res0[i] = res1[i];
            res1[i] = res2[i];
        }
        tek[0] = tek[0] + tau;
        tek[1] = tek[1] + tau;
        tek[2] = tek[2] + tau;
    }
    ans.close();
    delete[] x;
    delete[] tek;
    delete[] res1;
    delete[] res2;
    delete[] res0;
}
void test20()
{
    std::ofstream ans;
    ans.open(std::to_string(indD) + "k" + std::to_string(ct) + "v20.txt");
    ans << std::setprecision(15);
    h = 0.05;
    tau = ct * 0.01;
    a = 5;
    t0 = 0;
    u0 = 0;
    L = 1;
    Tmax = 10;
    indU = 3;
    n = int(L / h);
    double* x;
    double* res0;
    double* res1;
    double* res2;
    double* tek;
    res0 = new double[n + 1];
    res1 = new double[n + 1];
    res2 = new double[n + 1];
    x = new double[n + 1];
    tek = new double[3];
    tek[0] = t0;
    tek[1] = tek[0] + tau;
    tek[2] = tek[1] + tau;
    for (int i = 0; i < n + 1; i++)
    {
        x[i] = u0 + i * h;
        res0[i] = U0(x[i]);
    }
    if (indD == 0)
    {
        for (int i = 0; i < n + 1; i++)
        {
            res1[i] = res0[i] + tau * UT(x[i]) + (a * a * tau * tau / 2) * U0xx(x[i]);
        }
    }
    res1[0] = UL(tek[1]);
    res1[n] = UR(tek[1]);
    if (indD == 1)
    {
        for (int i = 1; i < n; i++)
        {
            res1[i] = res0[i] + tau * UT(x[i]) + (a * a * tau * tau / 2) * AU0xx(res0[i - 1], res0[i], res0[i + 1]);
        }
    }
    ans << h << std::endl;
    ans << tek[0] << ' ';
    for (int i = 0; i < n + 1; i++)
    {
        ans << res0[i] << ' ';
    }
    ans << std::endl;
    ans << tek[1] << ' ';
    for (int i = 0; i < n + 1; i++)
    {
        ans << res1[i] << ' ';
    }
    ans << std::endl;
    while (tek[2] < (Tmax + tau))
    {
        MStep(tek, x, res0, res1, res2, n);
        ans << tek[2] << ' ';
        for (int i = 0; i < n + 1; i++)
        {
            ans << res2[i] << ' ';
        }
        ans << std::endl;
        for (int i = 0; i < n + 1; i++)
        {
            res0[i] = res1[i];
            res1[i] = res2[i];
        }
        tek[0] = tek[0] + tau;
        tek[1] = tek[1] + tau;
        tek[2] = tek[2] + tau;
    }
    ans.close();
    delete[] x;
    delete[] tek;
    delete[] res1;
    delete[] res2;
    delete[] res0;
}
void Dtest1()
{
    std::ofstream ans;
    ans.open(std::to_string(indD) + "k" + std::to_string(ct) + "Dtest1.txt");
    ans << std::setprecision(15);
    h = 0.05;
    tau = 0.01;
    a = 5.;
    t0 = 0;
    u0 = -2;
    L = 4;
    Tmax = 10;
    indU = 4;
    n = int(L / h);
    double* x;
    double* res0;
    double* res1;
    double* res2;
    double* tek;
    res0 = new double[n + 1];
    res1 = new double[n + 1];
    res2 = new double[n + 1];
    x = new double[n + 1];
    tek = new double[3];
    tek[0] = t0;
    tek[1] = tek[0] + tau;
    tek[2] = tek[1] + tau;
    for (int i = 0; i < n + 1; i++)
    {
        x[i] = u0 + i * h;
        res0[i] = U0(x[i]);
        //std::cout << x[i] << ' '<< U0(x[i])<<std::endl;
    }
    if (indD == 0)
    {
        for (int i = 0; i < n + 1; i++)
        {
            res1[i] = res0[i] + tau * UT(x[i]) + (a * a * tau * tau / 2) * U0xx(x[i]);
        }
    }
    res1[0] = UL(tek[1]);
    res1[n] = UR(tek[1]);
    if (indD == 1)
    {
        for (int i = 1; i < n; i++)
        {
            res1[i] = res0[i] + tau * UT(x[i]) + (a * a * tau * tau / 2) * AU0xx(res0[i - 1], res0[i], res0[i + 1]);
        }
    }
    
    ans << h << std::endl;
    ans << tek[0] << ' ';
    for (int i = 0; i < n + 1; i++)
    {
        ans << res0[i] << ' ';
    }
    ans << std::endl;
    ans << tek[1] << ' ';
    for (int i = 0; i < n + 1; i++)
    {
        ans << res1[i] << ' ';
    }
    ans << std::endl;
    while (tek[2] < (Tmax + tau))
    {
        MStep(tek, x, res0, res1, res2, n);
        ans << tek[2] << ' ';
        for (int i = 0; i < n + 1; i++)
        {
            ans << res2[i] << ' ';
        }
        ans << std::endl;
        for (int i = 0; i < n + 1; i++)
        {
            res0[i] = res1[i];
            res1[i] = res2[i];
        }
        tek[0] = tek[0] + tau;
        tek[1] = tek[1] + tau;
        tek[2] = tek[2] + tau;
    }
    ans.close();
    delete[] x;
    delete[] tek;
    delete[] res1;
    delete[] res2;
    delete[] res0;
}
void Dtest2()
{
    std::ofstream ans;
    ans.open(std::to_string(indD) + "k" + std::to_string(ct) + "Dtest2.txt");
    ans << std::setprecision(15);
    h = 0.05;
    tau = ct * 0.01;
    a = 5;
    t0 = 0;
    u0 = -1;
    L = 2;
    Tmax = 10;
    indU = 5;
    n = int(L / h);
    double* x;
    double* res0;
    double* res1;
    double* res2;
    double* tek;
    res0 = new double[n + 1];
    res1 = new double[n + 1];
    res2 = new double[n + 1];
    x = new double[n + 1];
    tek = new double[3];
    tek[0] = t0;
    tek[1] = tek[0] + tau;
    tek[2] = tek[1] + tau;
    for (int i = 0; i < n + 1; i++)
    {
        x[i] = u0 + i * h;
        res0[i] = U0(x[i]);
        //std::cout << x[i] << ' '<< U0(x[i])<<std::endl;
    }
    if (indD == 0)
    {
        for (int i = 0; i < n + 1; i++)
        {
            res1[i] = res0[i] + tau * UT(x[i]) + (a * a * tau * tau / 2) * U0xx(x[i]);
        }
    }
    res1[0] = UL(tek[1]);
    res1[n] = UR(tek[1]);
    if (indD == 1)
    {
        for (int i = 1; i < n; i++)
        {
            res1[i] = res0[i] + tau * UT(x[i]) + (a * a * tau * tau / 2) * AU0xx(res0[i - 1], res0[i], res0[i + 1]);
        }
    }
    ans << h << std::endl;
    ans << tek[0] << ' ';
    for (int i = 0; i < n + 1; i++)
    {
        ans << res0[i] << ' ';
    }
    ans << std::endl;
    ans << tek[1] << ' ';
    for (int i = 0; i < n + 1; i++)
    {
        ans << res1[i] << ' ';
    }
    ans << std::endl;
    while (tek[2] < (Tmax + tau))
    {
        MStep(tek, x, res0, res1, res2, n);
        ans << tek[2] << ' ';
        for (int i = 0; i < n + 1; i++)
        {
            ans << res2[i] << ' ';
        }
        ans << std::endl;
        for (int i = 0; i < n + 1; i++)
        {
            res0[i] = res1[i];
            res1[i] = res2[i];
        }
        tek[0] = tek[0] + tau;
        tek[1] = tek[1] + tau;
        tek[2] = tek[2] + tau;
    }
    ans.close();
    delete[] x;
    delete[] tek;
    delete[] res1;
    delete[] res2;
    delete[] res0;
}
void Dtest3()
{
    std::ofstream ans;
    ans.open(std::to_string(indD) + "k" + std::to_string(ct) + "Dtest3.txt");
    ans << std::setprecision(15);
    h = 0.05;
    tau = ct * 0.01;
    a = 5;
    t0 = 0;
    u0 = 0;
    L = 4*pi;
    Tmax = 10;
    indU = 6;
    n = int(L / h);
    double* x;
    double* res0;
    double* res1;
    double* res2;
    double* tek;
    res0 = new double[n + 1];
    res1 = new double[n + 1];
    res2 = new double[n + 1];
    x = new double[n + 1];
    tek = new double[3];
    tek[0] = t0;
    tek[1] = tek[0] + tau;
    tek[2] = tek[1] + tau;
    for (int i = 0; i < n + 1; i++)
    {
        x[i] = u0 + i * h;
        res0[i] = U0(x[i]);
        //std::cout << x[i] << ' '<< U0(x[i])<<std::endl;
    }
    if (indD == 0)
    {
        for (int i = 0; i < n + 1; i++)
        {
            res1[i] = res0[i] + tau * UT(x[i]) + (a * a * tau * tau / 2) * U0xx(x[i]);
        }
    }
    res1[0] = UL(tek[1]);
    res1[n] = UR(tek[1]);
    if (indD == 1)
    {
        for (int i = 1; i < n; i++)
        {
            res1[i] = res0[i] + tau * UT(x[i]) + (a * a * tau * tau / 2) * AU0xx(res0[i - 1], res0[i], res0[i + 1]);
        }
    }
    ans << h << std::endl;
    ans << tek[0] << ' ';
    for (int i = 0; i < n + 1; i++)
    {
        ans << res0[i] << ' ';
    }
    ans << std::endl;
    ans << tek[1] << ' ';
    for (int i = 0; i < n + 1; i++)
    {
        ans << res1[i] << ' ';
    }
    ans << std::endl;
    while (tek[2] < (Tmax + tau))
    {
        MStep(tek, x, res0, res1, res2, n);
        ans << tek[2] << ' ';
        for (int i = 0; i < n + 1; i++)
        {
            ans << res2[i] << ' ';
        }
        ans << std::endl;
        for (int i = 0; i < n + 1; i++)
        {
            res0[i] = res1[i];
            res1[i] = res2[i];
        }
        tek[0] = tek[0] + tau;
        tek[1] = tek[1] + tau;
        tek[2] = tek[2] + tau;
    }
    ans.close();
    delete[] x;
    delete[] tek;
    delete[] res1;
    delete[] res2;
    delete[] res0;
}
void Aproxtest16()
{
    std::ofstream ans;
    ans.open("test163"+std::to_string(st+1)+".txt");
    ans << std::setprecision(15);
    h = 0.2/aprox;
    tau = 0.2/aprox;
    a = 0.5;
    t0 = 0;
    u0 = 0;
    L = 1;
    Tmax = 1;
    indU = 3;
    n = int(L / h);
    double* x;
    double* res0;
    double* res1;
    double* res2;
    double* tek;
    res0 = new double[n + 1];
    res1 = new double[n + 1];
    res2 = new double[n + 1];
    x = new double[n + 1];
    tek = new double[3];
    tek[0] = t0;
    tek[1] = tek[0] + tau;
    tek[2] = tek[1] + tau;
    for (int i = 0; i < n + 1; i++)
    {
        x[i] = u0 + i * h;
        res0[i] = U0(x[i]);
    }
    if (indD == 0)
    {
        for (int i = 0; i < n + 1; i++)
        {
            res1[i] = res0[i] + tau * UT(x[i]) + (a * a * tau * tau / 2) * U0xx(x[i]);
        }
    }
    res1[0] = UL(tek[1]);
    res1[n] = UR(tek[1]);
    if (indD == 1)
    {
        for (int i = 1; i < n; i++)
        {
            res1[i] = res0[i] + tau * UT(x[i]) + (a * a * tau * tau / 2) * AU0xx(res0[i - 1], res0[i], res0[i + 1]);
        }
    }
    ans << h << std::endl;
    ans << tek[0] << ' ';
    for (int i = 0; i < n + 1; i++)
    {
        ans << res0[i] << ' ';
    }
    ans << std::endl;
    ans << tek[1] << ' ';
    for (int i = 0; i < n + 1; i++)
    {
        ans << res1[i] << ' ';
    }
    ans << std::endl;
    while (tek[2] < (Tmax + tau))
    {
        MStep(tek, x, res0, res1, res2, n);
        ans << tek[2] << ' ';
        for (int i = 0; i < n + 1; i++)
        {
            ans << res2[i] << ' ';
        }
        ans << std::endl;
        for (int i = 0; i < n + 1; i++)
        {
            res0[i] = res1[i];
            res1[i] = res2[i];
        }
        tek[0] = tek[0] + tau;
        tek[1] = tek[1] + tau;
        tek[2] = tek[2] + tau;
    }
    ans.close();
    delete[] x;
    delete[] tek;
    delete[] res1;
    delete[] res2;
    delete[] res0;
}
void Aproxtest1()
{
    std::ofstream ans;
    ans.open("test11" + std::to_string(st + 1) + ".txt");
    ans << std::setprecision(15);
    h = 0.2/aprox;
    tau = 0.1/aprox;
    a = 1;
    t0 = 0;
    u0 = 0;
    L = 1;
    Tmax = 3;
    indU = 0;
    n = int(L / h);
    double* x;
    double* res0;
    double* res1;
    double* res2;
    double* tek;
    res0 = new double[n + 1];
    res1 = new double[n + 1];
    res2 = new double[n + 1];
    x = new double[n + 1];
    tek = new double[3];
    tek[0] = t0;
    tek[1] = tek[0] + tau;
    tek[2] = tek[1] + tau;
    for (int i = 0; i < n + 1; i++)
    {
        x[i] = u0 + i * h;
        res0[i] = U0(x[i]);
    }
    if (indD == 0)
    {
        for (int i = 0; i < n + 1; i++)
        {
            res1[i] = res0[i] + tau * UT(x[i]) + (a * a * tau * tau / 2) * U0xx(x[i]);
        }
    }
    res1[0] = UL(tek[1]);
    res1[n] = UR(tek[1]);
    if (indD == 1)
    {
        for (int i = 1; i < n; i++)
        {
            res1[i] = res0[i] + tau * UT(x[i]) + (a * a * tau * tau / 2) * AU0xx(res0[i - 1], res0[i], res0[i + 1]);
        }
    }
    ans << h << std::endl;
    ans << tek[0] << ' ';
    for (int i = 0; i < n + 1; i++)
    {
        ans << res0[i] << ' ';
    }
    ans << std::endl;
    ans << tek[1] << ' ';
    for (int i = 0; i < n + 1; i++)
    {
        ans << res1[i] << ' ';
    }
    ans << std::endl;
    while (tek[2] < (Tmax + tau))
    {
        MStep(tek, x, res0, res1, res2, n);
        ans << tek[2] << ' ';
        for (int i = 0; i < n + 1; i++)
        {
            ans << res2[i] << ' ';
        }
        ans << std::endl;
        for (int i = 0; i < n + 1; i++)
        {
            res0[i] = res1[i];
            res1[i] = res2[i];
        }
        tek[0] = tek[0] + tau;
        tek[1] = tek[1] + tau;
        tek[2] = tek[2] + tau;
    }
    ans.close();
    delete[] x;
    delete[] tek;
    delete[] res1;
    delete[] res2;
    delete[] res0;
}
int main()
{
    indD =1;
    ct = 1;
    aprox = 1;
    st = 0;
    /*for (int i = 1; i < 8; i++)
    {
        Aproxtest1();
        aprox = aprox * 2;
        st = st + 1;
    }*/
    
    
    //test1();
    //test2();
    //test16();
    //test20();
    Dtest1();
   // Dtest2();
   // Dtest3();
    std::cout << "Hello World!\n";
}

