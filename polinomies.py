#!/usr/bin/env python
#-*- coding: utf-8 -*-

# imports

# Functions

class Polynomies():

    def __init__(self, coeficientes):
        self.coeficientes = coeficientes
    
    def derivar(self):
        self.derivada = Polynomies([x*(i+1) for i,x in
            enumerate(self.coeficientes[1:])])

    def get_derivada(self):
        self.derivar()
        return self.derivada

    def integrar(self):
        self.integral = Polynomies([0] + [float(x)/(i+1)
            for i, x in enumerate(self.coeficientes)])

    def get_integral(self):
        self.integrar()
        return self.integral

    def evaluar(self,valor):
        return sum([x*valor**i for i, x in enumerate(self.coeficientes)])

    def __str__(self):
        poly = []
        polytoshow = ''
        coef = self.coeficientes
        for i,x in enumerate(coef):
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


def main():
    polinomio = Polynomies([2,-3,1])
    print 'polinomio: ', polinomio
    der = polinomio.get_derivada()
    print 'derivada: ', der
    inte = polinomio.get_integral()
    print 'integral: ', inte
    valor = raw_input('ingrese valor a evaluar: ')
    print 'polinomio evaluado en %s: %s' % (valor, polinomio.evaluar(float(valor)))

if __name__ == '__main__':
    main()
