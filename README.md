# Metoda Elementów Skończonych (MES) - Potencjał Elektromagnetyczny

Projekt zrealizowany w ramach przedmiotu **Równania Różniczkowe i Różnicowe**. Celem było rozwiązanie równania różniczkowego drugiego rzędu opisującego potencjał elektromagnetyczny przy użyciu metody elementów skończonych (MES).

## 1. Problem obliczeniowy
Rozwiązanie dotyczy problemu nr 4.5:
- **Równanie**: $\frac{d^2\phi}{dx^2} = -\frac{\rho}{\epsilon_r}$
- **Warunki brzegowe**:
  - Robin: $\phi'(0) - \phi(0) = 2$
  - Dirichlet: $\phi(4) = 1$
- **Parametry**: $\rho = 2$, $\epsilon_r$ jest funkcją przedziałami stałą.

## 2. Implementacja
Program został napisany w języku **C++** i realizuje następujące kroki:
1. **Sformułowanie wariacyjne**: Sprowadzenie równania do postaci słabej $b(w, v) = l(v)$.
2. **Kwadratury Gaussa-Legendre'a**: Całkowanie numeryczne (2-punktowe) w celu wyznaczenia macierzy sztywności.
3. **Algorytm Thomasa**: Efektywne rozwiązanie układu równań liniowych z macierzą trójpasmową.
4. **Eksport danych**: Wyniki zapisywane są do pliku `.csv` w celu łatwej wizualizacji.

## 3. Wyniki
Poniżej znajduje się wykres potencjału $\Phi(x)$ wygenerowany dla $n = 50$ elementów:

![Wizualizacja wyników](visualization.jpg)

## 4. Jak uruchomić
1. Skompiluj program (wymagany standard C++11 lub nowszy):
   ```bash
   g++ src/main.cpp -o mes_solver
