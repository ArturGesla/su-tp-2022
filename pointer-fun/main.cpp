#include <iostream>

void fun(double *a)
{
    std::cout << a[0] << std::endl;
}

int main()
{
    double a[3] = {1, 9, 3};
    double c = 10;
    double *b = &c;
    fun(b);
    // a[3]=1;
    return 0;
}