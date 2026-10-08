#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

const int MIN_ELEM = -20;
const int MAX_ELEM = 20;
const double EPS = 1e-6;

double random_nonzero();
void generate_run_through_matrix(vector<vector<double>>& A, int n);
void generate_jacobi_matrix(vector<vector<double>>& A, int n);
void generate_random_matrix(vector<vector<double>>& A, int n);

bool is_tridiagonal(vector<vector<double>> A, int n);
bool has_diagonal_dominance(vector<vector<double>> A, int n);

void print_matrix(vector<vector<double>> A, int n);

void gaussian_elimination(vector<vector<double>> A, int n);
bool run_through_method(vector<vector<double>> A, int n);
bool jacobi_method(vector<vector<double>> A, int n);

int main()
{
    setlocale(LC_ALL, "RU");

    int n;

    while (true)
    {
        cout << "Введите порядок СЛАУ: ";
        cin >> n;

        if (cin.fail() || n < 1 || n > 100)
        {
            cin.clear();
            cin.ignore(1000, '\n');
            cout << "Ошибка ввода.\n";
        }
        else break;
    }

    vector<vector<double>> A(n, vector<double>(n + 1));

    srand(time(0));

    cout << "\nВыберите способ задания системы:\n";
    cout << "1 — Под метод прогонки\n";
    cout << "2 — Под метод Якоби\n";
    cout << "3 — Произвольная\n";
    cout << "4 — Ввести вручную\n";

    int mode;
    cin >> mode;
    cin.ignore();

    switch (mode)
    {
        case 1:
            generate_run_through_matrix(A, n);
            break;

        case 2:
            generate_jacobi_matrix(A, n);
            break;

        case 3:
            generate_random_matrix(A, n);
            break;
        case 4:
        default:
        {
            for (int i = 0; i < n; i++) 
            { 
                string line; 
                cout << "Введите строку " << i + 1 << " (" << n << " коэффициентов + свободный член): "; 
                getline(cin, line); 
                if (line.erase(0, line.find_first_not_of(' ')).empty()) 
                { 
                    cout << "Строка сгенерирована: "; 
                    for (int j = 0; j <= n; j++) 
                    { 
                        A[i][j] = rand() % (MAX_ELEM - MIN_ELEM + 1) + MIN_ELEM;
                        cout << A[i][j] << " "; 
                    } 
                    cout << endl; 
                } 
                else 
                { 
                    stringstream ss(line); 
                    for (int j = 0; j < n + 1; j++) 
                    { 
                        ss >> A[i][j]; 
                    } 
                } 
            }
        }
    }

    cout << "\nМатрица системы:\n";
    print_matrix(A, n);

    if (run_through_method(A, n))
        return 0;

    if (jacobi_method(A, n))
        return 0;

    gaussian_elimination(A, n);
}

double random_nonzero()
{
    double x;

    do
    {
        x = rand() % (MAX_ELEM - MIN_ELEM + 1) + MIN_ELEM;
    } while (abs(x) < 1e-12);

    return x;
}

void generate_run_through_matrix(vector<vector<double>>& A, int n)
{
    cout << "\nСгенерирована трехдиагональная матрица\n";

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= n; j++)
            A[i][j] = 0;

        double left = 0;
        double right = 0;

        if (i > 0)
        {
            left = random_nonzero();
            A[i][i - 1] = left;
        }

        if (i < n - 1)
        {
            right = random_nonzero();
            A[i][i + 1] = right;
        }

        double extra = abs(random_nonzero());
        A[i][i] = abs(left) + abs(right) + extra;
        A[i][n] = random_nonzero();
    }
}

void generate_jacobi_matrix(vector<vector<double>>& A, int n)
{
    cout << "\nСгенерирована матрица с диагональным преобладанием\n";

    for (int i = 0; i < n; i++)
    {
        double sum = 0;

        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                A[i][j] = random_nonzero();
                sum += abs(A[i][j]);
            }
        }

        double extra = abs(random_nonzero());
        A[i][i] = sum + extra;
        A[i][n] = random_nonzero();
    }
}

void generate_random_matrix(vector<vector<double>>& A, int n)
{
    cout << "\nСгенерирована произвольная матрица\n";

    for (int i = 0; i < n; i++)
        for (int j = 0; j <= n; j++)
            A[i][j] = rand() % (MAX_ELEM - MIN_ELEM + 1) + MIN_ELEM;
}

bool is_tridiagonal(vector<vector<double>> A, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (abs(i - j) > 1 && abs(A[i][j]) > 1e-12)
                return false;

    return true;
}

bool has_diagonal_dominance(vector<vector<double>> A, int n)
{
    for (int i = 0; i < n; i++)
    {
        double sum = 0;

        for (int j = 0; j < n; j++)
            if (i != j)
                sum += abs(A[i][j]);

        if (abs(A[i][i]) <= sum)
            return false;
    }

    return true;
}

void print_matrix(vector<vector<double>> A, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << left << setw(10) << setfill(' ') << A[i][j];
        }
        cout << "|" << A[i][n] << endl;
    }

    cout << endl;
}

void gaussian_elimination(vector<vector<double>> A, int n)
{
    vector<vector<double>> Test(n, vector<double>(n + 1));
    Test = A;
    vector<double> X(n, 0);

    cout << "\nМЕТОД ПОЛНОГО ИСКЛЮЧЕНИЯ ПЕРЕМЕННЫХ:\n";

    // Прямой ход
    cout << "\nПрямой ход.\n";
    for (int i = 0; i < n; i++)
    {
        // Шаг 0
        int p = i;

        for (int k = i + 1; k < n; k++)
        {
            if (abs(A[k][i]) > abs(A[p][i]))
            {
                p = k;
            }
        }

        if (p != i)
        {
            swap(A[i], A[p]);
        }

        cout << "Шаг 0. Выбор ведущего элемента по столбцу:" << endl;
        print_matrix(A, n);

        // Шаг 1
        double d = A[i][i];

        if (abs(d) > 1e-12)
        {
            // Шаг 2
            for (int j = 0; j < n + 1; j++)
            {
                A[i][j] /= d;
            }

            cout << "Шаги 1, 2. Нормирование ведущего элемента:" << endl;
            print_matrix(A, n);
        }

        // Шаг 3
        for (int k = i + 1; k < n; k++)
        {
            double r = A[k][i];

            for (int j = 0; j < n + 1; j++)
            {
                A[k][j] -= A[i][j] * r;
            }
        }

        cout << "Шаг 3. Исключение переменных под главной диагональю:" << endl;
        print_matrix(A, n);

        // Шаг 4
        int rankA = 0;
        int rankAB = 0;

        for (int i = 0; i < n; i++)
        {
            bool zeroLeft = true;

            for (int j = 0; j < n; j++)
            {
                if (abs(A[i][j]) > 1e-12)
                {
                    zeroLeft = false;
                    break;
                }
            }

            if (!zeroLeft) rankA++;

            bool zeroFull = zeroLeft && abs(A[i][n]) <= 1e-12;

            if (!zeroFull) rankAB++;
        }

        // Если есть строки, где 0 равен const, не равной 0
        if (rankA < rankAB)
        {
            cout << "Система несовместна. Решений нет." << endl;
            return;
        }

        // Число уравнений меньше числа неизвестных переменных
        if (rankA < n)
        {
            cout << "Система имеет бесконечно много решений." << endl;
            return;
        }
    }

    // Обратный ход
    cout << "\nОбратный ход.\n";
    for (int i = n - 1; i >= 0; i--)
    {
        for (int k = i - 1; k >= 0; k--)
        {
            double r = A[k][i];

            for (int j = 0; j < n + 1; j++)
            {
                A[k][j] -= A[i][j] * r;
            }
        }

        cout << "Исключение переменных над главной диагональю:" << endl;
        print_matrix(A, n);
    }

    for (int i = 0; i < n; i++)
    {
        X[i] = A[i][n];
    }

    cout << "\nРешение системы:\n";

    for (int i = 0; i < n; i++)
    {
        cout << "x" << i + 1 << " = " << fixed << X[i] << endl;
    }

    cout << "\nПроверка:\n";

    for (int i = 0; i < n; i++)
    {
        double left = 0;
        for (int j = 0; j < n; j++)
        {
            left += Test[i][j] * X[j];
        }
        cout << left << " = " << fixed << Test[i][n] << endl;
    }
}

bool run_through_method(vector<vector<double>> A, int n)
{
    vector<vector<double>> Test = A;

    if (!is_tridiagonal(A, n))
    {
        cout << "\nМатрица не является трехдиагональной. Метод прогонки не применим.\n";
        return false;
    }

    vector<double> a(n), b(n), c(n), d(n);
    vector<double> P(n), Q(n), X(n);

    // Формирование коэффициентов системы a, b, c, d
    for (int i = 0; i < n; i++)
    {
        a[i] = (i > 0) ? A[i][i - 1] : 0;
        b[i] = A[i][i];
        c[i] = (i < n - 1) ? A[i][i + 1] : 0;
        d[i] = A[i][n];
    }

    for (int i = 1; i < n - 1; i++)
    {
        if (abs(a[i]) < 1e-12 || abs(c[i]) < 1e-12)
        {
            cout << "\nМетод прогонки неустойчив.\n";
            return false;
        }
    }

    for (int i = 0; i < n; i++)
    {
        if (abs(b[i]) < abs(a[i]) + abs(c[i]))
        {
            cout << "\nМетод прогонки неустойчив.\n";
            return false;
        }
    }

    cout << "\nМЕТОД ПРОГОНКИ:\n";

    cout << "\nПрямой ход.\n";

    P[0] = -c[0] / b[0];
    Q[0] = d[0] / b[0];

    cout << "P1 = " << P[0] << endl;
    cout << "Q1 = " << Q[0] << endl;

    for (int i = 1; i < n; i++)
    {
        double denom = b[i] + a[i] * P[i - 1];

        if (abs(denom) < 1e-12)
        {
            cout << "Ошибка метода прогонки: деление на ноль.\n";
            return false;
        }

        P[i] = -c[i] / denom;
        Q[i] = (d[i] - a[i] * Q[i - 1]) / denom;

        cout << "P" << i + 1 << " = " << P[i] << endl;
        cout << "Q" << i + 1 << " = " << Q[i] << endl;
    }

    cout << "\nОбратный ход.\n";

    X[n - 1] = Q[n - 1];
    cout << "x" << n << " = " << X[n - 1] << endl;

    for (int i = n - 2; i >= 0; i--)
    {
        X[i] = P[i] * X[i + 1] + Q[i];
        cout << "x" << i + 1 << " = " << X[i] << endl;
    }

    cout << "\nРешение системы:\n";

    for (int i = 0; i < n; i++)
        cout << "x" << i + 1 << " = " << fixed << X[i] << endl;

    cout << "\nПроверка:\n";

    for (int i = 0; i < n; i++)
    {
        double left = 0;

        for (int j = 0; j < n; j++)
            left += Test[i][j] * X[j];

        cout << left << " = " << fixed << Test[i][n] << endl;
    }

    return true;
}

bool jacobi_method(vector<vector<double>> A, int n)
{
    vector<vector<double>> Test = A;

    // Шаг 0. Перестановка строк для устранения нулей на диагонали
    for (int i = 0; i < n; i++)
    {
        if (abs(A[i][i]) < 1e-12)
        {
            bool swapped = false;

            for (int k = i + 1; k < n; k++)
            {
                if (abs(A[k][i]) > 1e-12)
                {
                    swap(A[i], A[k]);
                    swapped = true;
                    break;
                }
            }

            if (!swapped)
            {
                cout << "Невозможно устранить нулевой диагональный элемент.\n";
                return false;
            }
        }
    }

    // Шаг 1. Проверка условия сходимости
    if (!has_diagonal_dominance(A, n))
    {
        cout << "\nУсловие сходимости метода Якоби не выполняется.\n";
        return false;
    }

    // Шаг 2. Формирование матриц a и b
    vector<vector<double>> a(n, vector<double>(n, 0));
    vector<double> b(n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i != j)
                a[i][j] = -A[i][j] / A[i][i];
        }

        b[i] = A[i][n] / A[i][i];
    }

    // Шаг 3. Начальное приближение
    vector<double> X = b;
    vector<double> Xnew(n);

    const int maxIter = 10000;

    cout << "\nМЕТОД ЯКОБИ:\n";

    // Шаг 4. Итерационный процесс
    for (int iter = 1; iter <= maxIter; iter++)
    {
        cout << "\nИтерация " << iter << endl;

        for (int i = 0; i < n; i++)
        {
            double sum = 0;

            for (int j = 0; j < n; j++)
                sum += a[i][j] * X[j];

            Xnew[i] = sum + b[i];

            cout << "x" << i + 1 << " = " << Xnew[i] << endl;
        }

        // Шаг 5. Проверка критерия остановки
        double norm = 0;

        for (int i = 0; i < n; i++)
            norm = max(norm, abs(Xnew[i] - X[i]));

        if (norm < EPS)
        {
            cout << "\nРешение найдено с точностью " << EPS << endl;
            X = Xnew;
            break;
        }

        X = Xnew;
    }

    cout << "\nРешение системы:\n";

    for (int i = 0; i < n; i++)
        cout << "x" << i + 1 << " = " << fixed << X[i] << endl;

    cout << "\nПроверка:\n";

    for (int i = 0; i < n; i++)
    {
        double left = 0;

        for (int j = 0; j < n; j++)
            left += Test[i][j] * X[j];

        cout << left << " = " << fixed << Test[i][n] << endl;
    }

    return true;
}