#!/usr/bin/env python
#-*- coding: utf-8 -*-

# imports

from polinomies import Polinomio

import numpy as np
from numpy import linalg as ln

# constantes

cl = y = l = r = R1 = Rc = R2 = T1 = Tc = T2 = 1.0
ylr = y * (l / r) ** 2
rl = (r / l) ** 2 
polinomiosp = [Polinomio(x+1) for x in range(4)]
polinomiosq = [Polinomio(x+1) for x in range(4)]

class Matriz():
    
    def kaa1(self, i, m):
        p1i = polinomiosp[i-1]
        p1m = polinomiosp[m-1]
        dp1i = p1i.get_derivada()
        dp1m = p1m.get_derivada()
        prod_der = dp1i * dp1m
        prod_pol = p1i * p1m
        ylr_prod = prod_pol * ylr
        integ =  prod_der + ylr_prod
        intdef = integ.integral_definida(0, cl)
        return intdef + (R1 * p1i(0) * p1m(0)) + (Rc * p1i(cl) * p1m(cl))

    def kab1(self, i, j):
        p1i = polinomiosp[i-1]
        q1j = polinomiosq[j-1]
        integ = p1i * q1j.get_derivada() * ylr
        print integ
        return -1 * integ.integral_definida(0, cl)

    def kbb1(self, j, n):
        q1j = polinomiosq[i-1]
        q1n = polinomiosq[m-1]
        dq1j = q1j.get_derivada()
        dq1n = q1n.get_derivada()
        prod_der = dq1j * dq1n
        integ =  prod_der * ylr
        intdef = integ.integral_definida(0, cl)
        return intdef + (Tl * q1j(0) * q1n(0)) + (Tc * q1j(cl) * q1n(cl))

    def la1_lambda1(self, i):
        pi = polinomiosp[i-1]
        return pi(cl)

    def lb1_lambda2(self, j):
        qj = polinomiosq[j-1]
        return qj(cl)

     def kaa2(self, i, m):
        p2i = polinomiosp[i-1]
        p2m = polinomiosp[m-1]
        dp2i = p2i.get_derivada()
        dp2m = p2m.get_derivada()
        prod_der = dp2i * dp2m
        prod_pol = p2i * p2m
        ylr_prod = prod_pol * ylr
        integ =  prod_der + ylr_prod
        intdef = integ.integral_definida(cl, 1)
        return intdef + (R2 * p2i(1) * p2m(1))

    def kab2(self, i, j):
        p2i = polinomiosp[i-1]
        q2j = polinomiosq[j-1]
        integ = p2i * q2j.get_derivada() * ylr
        return -1 * integ.integral_definida(cl, 1)

    def kbb2(self, j, n):
        q2j = polinomiosq[i-1]
        q2n = polinomiosq[m-1]
        dq2j = qj.get_derivada()
        dq2n = qn.get_derivada()
        prod_der = dq2j * dq2n
        integ =  prod_der * ylr
        intdef = integ.integral_definida(cl, 1)
        return intdef + (T2 * q2j(1) * q2n(1))

    def la2_lambda1(self, i):
        p2i = polinomiosp[i-1]
        return -1 * p2i(cl)

    def lb2_lambda2(self, j):
        q2j = polinomiosq[j-1]
        return -1 * q2j(cl)

    def maa1(self, i, m):
        p1i = polinomiosp[i-1]
        p1m = polinomiosp[m-1]
        integ = p1i * p1m * rl
        return integ.integral_definida(0, cl)

    def mbb1(self, j, n):
        q1j = polinomiosq[i-1]
        q1n = polinomiosq[m-1]
        integ = q1j * q1n
        return integ.integral_definida(0, cl)

    def maa2(self, i, m):
        p2i = polinomiosp[i-1]
        p2m = polinomiosp[m-1]
        integ = p2i * p2m * rl
        return integ.integral_definida(cl, 1)

    def mbb2(self, j, n):
        q2j = polinomiosq[i-1]
        q2n = polinomiosq[m-1]
        integ = q2j * q2n
        return integ.integral_definida(cl, 1)
