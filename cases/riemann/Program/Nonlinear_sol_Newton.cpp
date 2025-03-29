#include <iostream>
#include <cmath>
#include <iomanip>
#include <tuple> // Для std::tuple и std::make_tuple

using namespace std;

// Объявления функций
double velocity(double x);
double f(double y);
double g(double y);
double velocity_derivative(double x);

// Функция для решения задачи нахождения u и p с использованием метода Ньютона
tuple<double, double> solve_nonlinear_shock() {
    double initial_guess = 0.5; // Начальное приближение
    double tolerance = 1e-4; // Точность
    //int max_iterations = 100; // Максимальное количество итераций

    
    double x0 = initial_guess; // Начальное приближение
    double x1 = x0;
    double x2 = 0;
    double u = 0.0, p = 0.0; // Инициализируем u и p


    while(abs(x2 - x1) > tolerance) {
        if (x2 > 0) x1 = x2;
        x2 = x1 - velocity(x1) / velocity_derivative(x1);
    }
    u = x2;
    p = f(u);

    if (isnan(u) || isnan(p)) {
        // Если метод не сошелся
        cout << "Newton's method did not converge." << endl; 
        u = 0.0;
        p = 0.0;
    } else {
        cout << "u = " << fixed << setprecision(8) << u << endl; // Выводим u
        cout << "p = " << fixed << setprecision(8) << p << endl; // Выводим p
    }
    return make_tuple(u, p); // Возвращаем кортеж (u, p)

    // Если метод не сходится, возвращаем NaN для u и p
    //return make_tuple(NAN, NAN);
}


// Производная функции velocity (вам нужно вычислить это)
double velocity_derivative(double x) {
    double h = 0.00001;
    return (velocity(x + h) - velocity(x)) / h;
}

// Другие функции (f, g, velocity) остаются прежними
double velocity(double x) {
    return x * x - g(x);
}

// Функция f (давление)
double f(double y) {
    double rho_left = 1.0;
    double p_left = 1.0;
    //double rho_right = 0.125;
    //double p_right = 0.1;
    double gamma = 1.4;

    double c_left = sqrt(gamma * p_left / rho_left);
    double tmp = 1.0 - (gamma - 1.0) * y / (2.0 * c_left);

    return p_left * pow(tmp, (2.0 * gamma) / (gamma - 1.0));
}

// Функция g (связанная со скоростью)
double g(double y) {
    //double rho_left = 1.0;
    //double p_left = 1.0;
    double rho_right = 0.125;
    double p_right = 0.1;
    double gamma = 1.4;

    double tmp1 = pow((f(y) - p_right), 2);
    double tmp2 = (gamma - 1.0) * p_right + (gamma + 1.0) * f(y);

    return 2.0 * tmp1 / rho_right / tmp2;
}

/*
int main() {
    
    tuple<double, double> result = solve_nonlinear_shock(); // Вызываем функцию
    double u = get<0>(result); // Получаем u из кортежа
    double p = get<1>(result); // Получаем p из кортежа

    if (isnan(u) || isnan(p)) {
        cout << "Newton's method did not converge." << endl; // Если метод не сошелся
    } else {
        cout << "u = " << fixed << setprecision(8) << u << endl; // Выводим u
        cout << "p = " << fixed << setprecision(8) << p << endl; // Выводим p
    }

    return 0;
}
*/