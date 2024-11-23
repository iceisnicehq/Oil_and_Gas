####################################################################
# Основы нефтегазового дела: переработка и логистика углеводородов |
#                   Домашнее задание №1                            |
#                       Вариант №52                                |
#__________________________________________________________________|
#                                                         Выполнил:|
#                                           Студент группы КС-22-03|
#                                          Сабиров Данияр Серикович|
#                                                         Проверил:|
#                                   преподаватель кафедры газохимии|
#                                     Кондратенко Андрей Дмитриевич|
#__________________________________________________________________|
#                       Москва, 2024                               |
#__________________________________________________________________|

from math import sqrt, exp, sinh, cosh
from consts import R, Lt, \
    Ki, Gi, Qi, Fi, Ei, Si, Wi, \
    Kij, Gij, Eij, Vij, \
    B0i, C0i, D0i, E0i, F0i, G0i, H0i, I0i, J0i, \
    a, b, c, u, k, g, q, f, s, w, zci

# ИСХОДНЫЕ ДАННЫЕ (52)
natural_gas = {
    "methane":     0.887647,
    "ethane":      0.039515,
    "propane":     0.023270,
    "u_butane":    0.004948,
    "h_butane":    0.002650,
    "neo_pentane": 0.0,
    "u_pentane":   0.0,
    "h_pentane":   0.000595,
    "neo_hexane":  0.0,
    "h_hexane":    0.000753,
    "h_heptane":   0.0,
    "nitrogen":    0.032885,
    "c_dioxide":   0.007254,
    "helium":      0.000483,
    "hydgrogen":   0.0
}

pressure    = 20
celsius     = -1

def main(ri: list[float]) -> float:
    # (1) 
    denominator = sum(ri_val / zci_val for ri_val, zci_val in zip(ri, zci))
    mix = [
        (ri_val / zci_val) / denominator
        for ri_val, zci_val in zip(ri, zci)
    ]

    # (2)
    kelvin = celsius + 273.15

    # (3)
    exponent = 5 / 2

    Sum_xi_ki = sum(
        xi * (ki ** exponent)
        for xi, ki in zip(mix, Ki)
    ) ** 2
    DSum_xi_xj = sum(
        mix[i] * mix[j] * (Kij[i][j] ** 5 - 1) * (Ki[i] * Ki[j]) ** exponent
        for i in range(len(mix) - 1)
        for j in range(i + 1, len(mix))
    )
    Kx = (Sum_xi_ki + 2 * DSum_xi_xj) ** (1 / 5)

    # (4)
    reduced_temp = kelvin / Lt

    # (6)
    p0m = 10 ** (-3) * (Kx ** (-3) * R * Lt)
    
    # (5)
    reduced_pres = pressure / p0m

    # (7)-(17)
    G = sum(x * g for x, g in zip(mix, Gi)) + sum(
        mix[i] * mix[j] * (Gij[i][j] - 1) * (Gi[i] + Gi[j]) 
        for i in range(len(mix) - 1)
        for j in range(i + 1, len(mix))
    )
    Q = sum(x * q for x, q in zip(mix, Qi))
    F = sum((x**2) * f for x, f in zip(mix, Fi))
    V = (
        sum(x * (e ** exponent) for x, e in zip(mix, Ei)) ** 2 +
        2 * sum(
            mix[i] * mix[j] * (Vij[i][j] ** 5 - 1) * ((Ei[i] * Ei[j]) ** exponent)
            for i in range(len(mix) - 1)
            for j in range(i + 1, len(mix))
        )
    ) ** (1 / 5)
    Cn = lambda n: (
        ((G + 1 - g[n]) ** g[n]) *
        ((Q ** 2 + 1 - q[n]) ** q[n]) *
        ((F + 1 - f[n]) ** f[n]) *
        (V ** u[n])
    )
    U = [0] * 12 + [Cn(n) for n in range(12, 58)]
    Bn = lambda n: sum(
        mix[i] * mix[j] * Bnij(n, i, j) * Eijn(i, j) ** u[n] *
        ((Ki[i] * Ki[j]) ** (3 / 2))
        for i in range(len(mix))
        for j in range(len(mix))
    )
    Bnij = lambda n, i, j: (
        ((Gijn(i, j) + 1 - g[n]) ** g[n]) *
        ((Qi[i] * Qi[j] + 1 - q[n]) ** q[n]) *
        ((sqrt(Fi[i] * Fi[j]) + 1 - f[n]) ** f[n]) *
        ((Si[i] * Si[j] + 1 - s[n]) ** s[n]) *
        ((Wi[i] * Wi[j] + 1 - w[n]) ** w[n])
    )
    Eijn = lambda i, j: Eij[i][j] * sqrt(Ei[i] * Ei[j])         
    Gijn = lambda i, j: (Gij[i][j] * (Gi[i] + Gi[j])) / 2
    D = (
        [Bn(n) * (Kx ** -3) for n in range(12)] +
        [Bn(n) * (Kx ** -3) - Cn(n) for n in range(12, 18)] +
        [0] * 40
    )

    # (18)
    delta0 = (10 ** 3 * pressure * Kx ** 3) / (R * kelvin) 

    A = [
        # (21) = A0
        lambda density: sum(
            a[n] * density ** b[n] * reduced_temp ** (-u[n]) * 
            (b[n] * D[n] + (b[n] - c[n] * k[n] * density ** k[n]) * U[n] * exp(-c[n] * density ** k[n]))
            for n in range(58)
        ),
        # (22) = A1
        lambda density: sum(
            a[n] * density ** b[n] * reduced_temp ** (-u[n]) * (
                (b[n] + 1) * b[n] * D[n] + 
                ((b[n] - c[n] * k[n] * density ** k[n]) * 
                (b[n] - c[n] * k[n] * density ** k[n] + 1) - 
                c[n] * k[n] ** 2 * density ** k[n]) 
                * U[n] * exp(-c[n] * density ** k[n])
            )
            for n in range(58)
        ),
        # (25) = A2
        lambda density: sum(
            a[n] * density ** b[n] * reduced_temp ** (-u[n]) * (1 - u[n]) * 
            (b[n] * D[n] + (b[n] - c[n] * k[n] * density ** k[n]) * U[n] * exp(-c[n] * density ** k[n]))
            for n in range(58)
        ),
        # (26) = A3
        lambda density: sum(
            a[n] * density ** b[n] * reduced_temp ** (-u[n]) * u[n] * (1 - u[n]) * 
            (D[n] + U[n] * exp(-c[n] * density ** k[n]))
            for n in range(58)
        )
    ]
        
    # (24)
    calc_reduced_pres = lambda density: density * reduced_temp * (1 + A[0](density))

    # (23-24)
    while abs((calc_reduced_pres(delta0) - reduced_pres) / reduced_pres) > 10**(-6):
        A0_delta0 = A[0](delta0)
        A1_delta0 = A[1](delta0)
        delta0 += ((reduced_pres / reduced_temp) - (1 + A0_delta0) * delta0) / (1 + A1_delta0)

    A0, A1, A2, A3 = A[0](delta0), A[1](delta0), A[2](delta0), A[3](delta0)

    # (27)
    z = 1 + A0

    # (28)-(29)
    theta = reduced_temp**(-1)
    cp0ri = lambda i: (
         B0i[i] + 
        (C0i[i] * ((D0i[i] * theta) / sinh(D0i[i] * theta)) ** 2 if sinh(D0i[i] * theta) != 0 else 0) + 
        (E0i[i] * ((F0i[i] * theta) / cosh(F0i[i] * theta)) ** 2 if cosh(F0i[i] * theta) != 0 else 0) + 
        (G0i[i] * ((H0i[i] * theta) / sinh(H0i[i] * theta)) ** 2 if sinh(H0i[i] * theta) != 0 else 0) + 
        (I0i[i] * ((J0i[i] * theta) / cosh(J0i[i] * theta)) ** 2 if cosh(J0i[i] * theta) != 0 else 0)
    )
    cp0r = sum(mix[i] * cp0ri(i) for i in range(len(ri)))

    # (30)
    adiab_ind = (1 + A1 + ((1 + A2) ** 2) / (cp0r - 1 + A3)) / z
    
    return adiab_ind

if __name__ == "__main__":
    ri = list(natural_gas.values())
    print(f'Показатель адиабаты:\nk = {main(ri)}')
