#!/usr/bin/env python
#-*- coding: utf-8 -*-

# imports

from polinomies import Polinomio

import numpy as np
from numpy import linalg as ln

# constantes

cl = 0.5
y = 1.5
l = 1.0
r = 0.02887 
R1 = 0.0
Rc = 0.0
R2 = 0.0
T1 = 0.0
Tc = 0.0
T2 = 1.0
ylr = y * (l / r) ** 2
rl = (r / l) ** 2

# Condiciones de borde

CONDICIONES_BORDE = {
                1: ([1], [0, 1], [1], [-1, 1]),
                2: ([1], [0, 1], [1], [1]),
                3: ([1], [1], [1], [1]),
                4: ([0, 1], [0, 1], [-1, 1], [-1, 1]),
                5: ([0, 1], [0, 1], [1], [-1, 1]),
                6: ([0, 1], [0, 1], [1], [1])
}


class Matriz():

    def __init__(self, size=5, value=4, borde=1):
        matrix_size = size ** 2
        self.matriz_k = ln.linalg.zeros((matrix_size,matrix_size)) # TODO : arreglar esto
        self.matriz_m = ln.linalg.zeros((matrix_size,matrix_size))
        self.value = value
        self.size = size
        self.genera_polinomios(borde, value)
        #self.fill_matriz()

    def genera_polinomios(self, borde, value):
        """Genera los polinomios de acuerdo a las condiciones de
        borde que se ingrese"""
        try:
            p1, q1, p2, q2 = CONDICIONES_BORDE[borde]
            self.polinomios_p1 = [Polinomio(p1).aumenta_grado(x)
                for x in xrange(value)]
            self.polinomios_q1 = [Polinomio(q1).aumenta_grado(x)
                for x in xrange(value)]
            self.polinomios_p2 = [Polinomio(p2).aumenta_grado(x)
                for x in xrange(value)]
            self.polinomios_q2 = [Polinomio(q2).aumenta_grado(x)
                for x in xrange(value)]

        except ValueError:
            print 'Condiciones de borde erroneas(1 - 6)'

    def fill_matriz(self):
        """Llena la matriz con los valores relevantes"""

        n = self.size
        # Matriz K de rigidez

        self.matriz_k[0:n,0:n] = self.sub_matrices(self.kaa1)
        self.matriz_k[n:2 * n,n:n] = self.matriz_k[n:n,n:2 * n] = self.sub_matrices(self.kab1)
        self.matriz_k[2 * n:n,2 * n:n] = self.sub_matrices(self.kbb1)
        self.matriz_k[n:2 * n,4 * n:n * n] = self.matriz_k[n * n:,1] = self.sub_matrices(self.lb1_lambda2)
        self.matriz_k[2,2] = self.sub_matrices(self.kaa2)
        self.matriz_k[2,3] = self.matriz_k[3,2] = self.sub_matrices(self.kab2)
        self.matriz_k[3,3] = self.sub_matrices(self.kbb2)
        self.matriz_k[3,4] = self.matriz_k[4,3] = self.sub_matrices(self.lb2_lambda2)

        # Matriz M de masa

        self.matriz_m[0,0] = self.sub_matrices(self.maa1)
        self.matriz_m[1,1] = self.sub_matrices(self.mbb1)
        self.matriz_m[2,2] = self.sub_matrices(self.maa2)
        self.matriz_m[3,3] = self.sub_matrices(self.mbb2)

    def sub_matrices(self, func):
        """Falta implementar"""
        matriz = []
        try:
            for i in range(self.value):
                filas = []
                for j in range(self.value):
                    filas.append(func(i,j))
                matriz.append(filas)
        except TypeError:
            for i in range(self.value):
                filas = [func(i)]
                matriz.append(filas)
        sub_matriz = np.matrix(matriz)
        return sub_matriz

    def kaa1(self, i, m):
        p1i = self.polinomios_p1[i-1]
        p1m = self.polinomios_p1[m-1]
        dp1i = p1i.get_derivada()
        dp1m = p1m.get_derivada()
        prod_der = dp1i * dp1m
        prod_pol = p1i * p1m
        ylr_prod = prod_pol * ylr
        integ =  prod_der + ylr_prod
        intdef = integ.integral_definida(0, cl)
        return intdef + (R1 * p1i(0) * p1m(0)) + (Rc * p1i(cl) * p1m(cl))

    def kab1(self, i, j):
        p1i = self.polinomios_p1[i-1]
        q1j = self.polinomios_q1[j-1]
        integ = p1i * q1j.get_derivada() * ylr
        print integ
        return -1 * integ.integral_definida(0, cl)

    def kbb1(self, j, n):
        q1j = self.polinomios_q1[i-1]
        q1n = self.polinomios_q1[m-1]
        dq1j = q1j.get_derivada()
        dq1n = q1n.get_derivada()
        prod_der = dq1j * dq1n
        integ =  prod_der * ylr
        intdef = integ.integral_definida(0, cl)
        return intdef + (Tl * q1j(0) * q1n(0)) + (Tc * q1j(cl) * q1n(cl))

    def la1_lambda1(self, i):
        pi = self.polinomios_p1[i-1]
        return pi(cl)

    def lb1_lambda2(self, j):
        qj = self.polinomios_q1[j-1]
        return qj(cl)

    def kaa2(self, i, m):
        p2i = self.polinomios_p2[i-1]
        p2m = self.polinomios_p2[m-1]
        dp2i = p2i.get_derivada()
        dp2m = p2m.get_derivada()
        prod_der = dp2i * dp2m
        prod_pol = p2i * p2m
        ylr_prod = prod_pol * ylr
        integ =  prod_der + ylr_prod
        intdef = integ.integral_definida(cl, 1)
        return intdef + (R2 * p2i(1) * p2m(1))

    def kab2(self, i, j):
        p2i = self.polinomios_p2[i-1]
        q2j = self.polinomios_q2[j-1]
        integ = p2i * q2j.get_derivada() * ylr
        return -1 * integ.integral_definida(cl, 1)

    def kbb2(self, j, n):
        q2j = self.polinomios_q2[i-1]
        q2n = self.polinomios_q2[m-1]
        dq2j = qj.get_derivada()
        dq2n = qn.get_derivada()
        prod_der = dq2j * dq2n
        integ =  prod_der * ylr
        intdef = integ.integral_definida(cl, 1)
        return intdef + (T2 * q2j(1) * q2n(1))

    def la2_lambda2(self, i):
        p2i = self.polinomios_p2[i-1]
        return -1 * p2i(cl)

    def lb2_lambda2(self, j):
        q2j = self.polinomios_q2[j-1]
        return -1 * q2j(cl)

    def maa1(self, i, m):
        p1i = self.polinomios_p1[i-1]
        p1m = self.polinomios_p1[m-1]
        integ = p1i * p1m * rl
        return integ.integral_definida(0, cl)

    def mbb1(self, j, n):
        q1j = self.polinomios_q1[i-1]
        q1n = self.polinomios_q1[m-1]
        integ = q1j * q1n
        return integ.integral_definida(0, cl)

    def maa2(self, i, m):
        p2i = self.polinomios_p2[i-1]
        p2m = self.polinomios_p2[m-1]
        integ = p2i * p2m * rl
        return integ.integral_definida(cl, 1)

    def mbb2(self, j, n):
        q2j = self.polinomios_q2[i-1]
        q2n = self.polinomios_q2[m-1]
        integ = q2j * q2n
        return integ.integral_definida(cl, 1)
