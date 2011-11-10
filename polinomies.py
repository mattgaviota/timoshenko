#!/usr/bin/env python
#-*- coding: utf-8 -*-

# imports

# Clases


class Polinomio():

    def __init__(self, coeficientes=[]):
        """Inicializa los coeficientes del polinomio.
        Es en orden ascendente, es decir [A0, A1, A2,...,An]"""
        self.coeficientes = coeficientes[:]

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
        for i, x in enumerate(coef):
            if x and x != 1:
                if i == 1:
                    elem = '%sX' % (x,)
                elif i == 0:
                    elem = '%s' % (x,)
                else:
                    elem = '%sX^%s' % (x, i)
                poly.append(elem)
            elif x == 1:
                if i == 1:
                    elem = 'X'
                elif i == 0:
                    elem = '%s' % (x,)
                else:
                    elem = 'X^%s' % (i,)
                poly.append(elem)
        poly.reverse()
        for elem in poly:
            if elem.startswith('-'):
                polytoshow = polytoshow[:-2] + elem + '+ '
            else:
                polytoshow = polytoshow + elem + '+ '
        return polytoshow[:-2]
