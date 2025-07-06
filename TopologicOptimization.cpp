// TopologicOptimization.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <fstream>
#include <cmath>
//#include <iostream>
//#include <string>
#include "rayTrace.cpp"
//#include "ray.cpp"

using namespace std;



 const int N = Num;

double sqr(double x) {
    return x * x;
};



double sqrSum(double* A) {
    double r = 0;
    for (int i = 0;i < N;i++) {
        r += sqr(A[i]);
    }
    return r;
}
double f(double* x);
double kroneker(int a, int b) {//определяем символ кронекера
    if (a == b) return 1; else return 0;
}
double gradf(double* x, double h, int k) {//определяем операцию нахождения градиента в конечных разностях
    double df;
    double xp[N];
    double xn[N];
        for (int j = 0;j < N;j++) { xp[j] = x[j] + h * kroneker(k, j); }
        for (int j = 0;j < N;j++) { xn[j] = x[j] - h * kroneker(k, j); }
        df = (f(xp) - f(xn)) / (2 * h);//определяем градиент функции f 
    return df;
}
double* gradFall(double* x, double h, double b) {
    double A[N];
    for (int i = 0;i < N;i++)  A[i] = gradf(x, h, i);
       
    double NORMA = sqrt(sqrSum(A));

    for (int i = 0; i < N;i++) A[i] = A[i] / NORMA;
    
    static double xCurent[N];
    static double x1dCurent[N];
    static double x2dCurent[N];

    for (int i = 0; i < N;i++) xCurent[i] = x[i];
    for (int i = 0;i < N;i++) x1dCurent[i] = x[i] + A[i] * h;
    for (int i = 0;i < N;i++) x2dCurent[i] = x[i] + A[i] * 2 * h;

    double f1 = f(xCurent);
    double f2 = f(x1dCurent);
    double f3 = f(x2dCurent);

    std::cout << std::to_string(sqrt(sqrSum(A))) << std::endl;
 
    vector R;

    R.set(f1, (4 * f2 - f1 * 3 - f3) / (2 * h), (f2 - f1) / (h * h) - (4 * f2 - f1 * 3 - f3) / (2 * h * h));

    while (sqrt(sqrSum(A)) > b) {



        for (int i = 0;i < N;i++) xCurent[i] -= A[i] * h;
       
        for (int i = 0;i < N;i++) A[i] = gradf(xCurent, h, i);
      
        cout << to_string(f(xCurent)) << endl;
    }
    return xCurent;
}
rayTrace::system mirror;

int main()
{
    setlocale(LC_ALL, "rus");
   // ifstream surfaceData("A.txt");
    
    vector a;
    vector b;
    vector r;
    plane P;
    string T[Num];
    //double t[Num];

    cout << "переменные инициализированы" << endl;

    /*for (int i = 0;i < Num;i++) {
        surfaceData >> T[i];
    }
    surfaceData.close();
    for (int i = 0;i < Num;i++) {
        t[i] = stod(T[i]);
    }
    cout << "стрелки прогиба импротированы" << endl;*/
    a.set(1, 0, 0);
    b.set(0, 1, 0);
    r.set(0, 0, 0);
    P.set(a, b, r);

    //axisSpherSurface S;
    ///*S.setGrid(0, 10);
    //S.setZ(50);
    //S.setT(t);*/
    //S.set(-100, 50);
    //cout << "поверхность собрана" << endl;
    const int numberRay = 12;
    double d = 10.f / (numberRay - 1);
    vector e[numberRay][numberRay];
    for (int i = 0;i < numberRay;i++) {
        for (int j = 0;j < numberRay;j++) {
            e[i][j].set(-5.f + i * d, -5.f + j * d, 0);
            // cout << to_string(i * 10 + j) << endl;
        }
    }
    cout << "сетка входного зрачка посчитана" << endl;
    vector D;
    D.set(0, 0, 1);
    emiter E;
    for (int i = 0;i < numberRay;i++) {
        for (int j = 0;j < numberRay;j++) {
            E.E[i][j].set(D, e[i][j], 1);
        }
    }
    vector h;
    h.set(1, 0, 50);
    ray test;
    test.set(D, h, 1);
    /*for (int i = 0;i < numberRay;i++) {
        for (int j = 0;j < numberRay;j++) {
            E.E[i][j].set(D, h, 1);
        }
    }*/

    bool err;
    cout << "эмиттер входного зрачка посчитан" << endl;
    
    cout << "система ининциализирована" << endl;
    mirror.setEntranceEmiter(E);
    cout << "эмиттер входного зрачка загружен в систему" << endl;
    mirror.setSurface(-100, 200, 0);
    cout << "поверхность зеркала загружена в систему" << endl;
    mirror.imagePlane.set(a, b, r);
    cout << "плоскость изображения загружена в систему" << endl;
    mirror.n[0] = 1;
    mirror.n[1] = 2;
    cout << "показатели преломления загружены в систему" << endl;
    mirror.emiterTrace();
    cout << "лучики протрасировались" << endl;
    
   /* double* res = gradFall(mirror.surface1.t, 0.001, 0.001);*/
    
    ofstream surfaceDataOut("D.txt");
    /*for (int i = 0;i < N;i++) {
        surfaceDataOut << to_string(res[i]) << '\t' << to_string(S.R[i]) << endl;
    }*/
    emiter res = mirror.E[2];

    //cout << to_string(mirror.Surface[0].shift.z) << endl;
    
    for (int i = 0;i < numberRay;i++) {
        for (int j = 0;j < numberRay;j++) {
            //cout << to_string(res.E[i][j].R.x) << '\t' << to_string(res.E[i][j].R.y) << '\t' << to_string(res.E[i][j].R.z) << endl;
            //surfaceDataOut << to_string(res.E[i][j].T.x) << '\t' << to_string(res.E[i][j].T.y) << '\t' << to_string(res.E[i][j].T.z) << endl;
            surfaceDataOut << to_string(res.E[i][j].R.x) << '\t' << to_string(res.E[i][j].R.y) << '\t' << to_string(res.E[i][j].R.z) << endl;
        }
    }
    surfaceDataOut.close();
    
    std::cout << "евреи спиздили наши лучики" << std::endl;

   //Res = 0;
   //Res.zero();
   //std::cout << std::to_string(Res.T.Norm()) << std::endl;
   //std::cout << std::to_string(L.T.x) << '\t' << std::to_string(L.T.y) << '\t' << std::to_string(L.T.z) << std::endl;
   //std::cout << std::to_string(Res.T.x) << '\t' << std::to_string(Res.T.y) << '\t' << std::to_string(Res.T.z) << std::endl;
}   
double f(double* x) {
    //mirror.set1Surface(50, x, 10, 0);
    mirror.emiterTrace();
    double res = 0;
    for (int i = 0;i < mirror.E->numberRay;i++) {
        for (int j = 0;j < mirror.E->numberRay;j++) {
            res += abs(mirror.E[2].E[i][j].R.x) + abs(mirror.E[2].E[i][j].R.y);
        }
    }
    return res;
}
// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
