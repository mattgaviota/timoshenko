#!/usr/bin/env python
#-*- coding: utf-8 -*-

# imports

# Clases


class Polinomio():

    def __init__(self, args=[]):
        """Inicializa los coeficientes del polinomio si args es una lista,
        la lista debe representar un polinomio en orden ascendente en el que
        solo se pasa los coeficientes del mismo, para obtener el polinomio
        deseado este debe estar completo, es decir [A0, A1, A2,...,An]
        Si args es un número se genera un polinomio de la forma X^N(X-1)
        donde N es args y usa el metodo generar.
        Si no se le pasa parametros se genera un polinomio de la forma
        X(X-1)"""
        if type(args) == type([]):
            if args:
                self.coeficientes = args[:]
            else:
                self.coeficientes = []
                self.generar(1)
        else:
            try:
                self.coeficientes = []
                self.generar(args)
            except TypeError:
                print 'Argumentos invalidos(lista o número del grado)'

    def generar(self, grado):
        """Genera un polinomio de la forma X^N(X-1) donde N es un
        parametro opcional al generar el Polinomio"""
        self.coeficientes = [0] * grado + [-1, 1]

    def aumenta_grado(self):
        """Aumenta el grado del polinomio"""
        self.coeficientes = [0] + self.coeficientes

    def derivar(self):
        """Deriva el polinomio"""
        self.derivada = Polinomio([x * (i + 1) for i, x in
            enumerate(self.coeficientes[1:])])

    def get_derivada(self):
        """Retorna la derivada del polinomio como otro polinomio"""
        self.derivar()
        return self.derivada

    def integrar(self):
        """Integra el polinomio"""
        self.integral = Polinomio([0] + [float(x) / (i + 1)
            for i, x in enumerate(self.coeficientes)])

    def get_integral(self):
        """Retorna la integral del polinomio como otro polinomio"""
        self.integrar()
        return self.integral

    def evaluar(self, valor):
        """Metodo para evaluear el polinomio"""
        return sum([x * valor ** i for i, x in enumerate(self.coeficientes)])

    def __add__(self, other):
        """Sobrecarga de la suma.
        Se asume que other también es un polinomio"""
        r = [a + b for a, b in zip(self.coeficientes, other.coeficientes)]
        if len(self.coeficientes) < len(other.coeficientes):
            r += other.coeficientes[len(self.coeficientes):]
        else:
            r += self.coeficientes[len(other.coeficientes):]
        return Polinomio(r)

    def __sub__(self, other):
        """Sobrecarga de la resta.
        Se asume que other también es un polinomio"""
        r = [a - b for a, b in zip(self.coeficientes, other.coeficientes)]
        if len(self.coeficientes) < len(other.coeficientes):
            r += [-1 * v for v in other.coeficientes[len(self.coeficientes):]]
        else:
            r += self.coeficientes[len(other.coeficientes):]
        return Polinomio(r)

    def __mul__(self, other):
        """Sobrecarga del producto.
        Se asume que other también es un polinomio."""
        try:
            r = [0] * (len(self.coeficientes) + len(other.coeficientes) - 1)
            for sgrado, scoef in enumerate(self.coeficientes):
                for ogrado, ocoef in enumerate(other.coeficientes):
                    r[sgrado + ogrado] += scoef * ocoef
        except AttributeError:
            r = [other * coef for coef in self.coeficientes]
        return Polinomio(r)

    def __eq__(self, other):
        """Comparación de igualdad.
        Se asume que other también es un polinomio."""
        for a, b in zip(self.coeficientes, other.coeficientes):
            if a != b:
                return False
        for c in self.coeficientes[len(other.coeficientes):]:
            if c != 0:
                return False
        for c in other.coeficientes[len(self.coeficientes):]:
            if c != 0:
                return False
        return True

    def __ne__(self, other):
        """Comparación de desigualdad implementada
        como complementaria de la de igualdad."""
        return not self.__eq__(other)

    def __call__(self, x):
        """Evaluación en un punto con notación funcional."""
        return self.evaluar(x)

    def __str__(self):
        """Metodo para mostrar el polinomio"""
        poly = []
        polytoshow = ''
        coef = self.coeficientes
        for grado, coef in enumerate(coef):
            if coef and coef != 1:
                if grado == 1:
                    elem = '%sX' % (coef,)
                elif grado == 0:
                    elem = '%s' % (coef,)
                else:
                    elem = '%sX^%s' % (coef, grado)
                poly.append(elem)
            elif coef == 1:
                if grado == 1:
                    elem = 'X'
                elif grado == 0:
                    elem = '%s' % (coef,)
                else:
                    elem = 'X^%s' % (grado,)
                poly.append(elem)
        poly.reverse()
        for elem in poly:
            if elem.startswith('-'):
                if elem.startswith('-1X'):
                    polytoshow = polytoshow[:-1] + elem.replace('1', '') + '+'
                else:
                    polytoshow = polytoshow[:-1] + elem + '+'
            else:
                polytoshow = polytoshow + elem + '+'
        return polytoshow[:-1]
