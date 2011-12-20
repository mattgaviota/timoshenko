Interfaz
========

Notación
--------

    +---------+---------------------+----+
    |cl = c/l |Parametro geometrico |dato|
    +---------+---------------------+----+
    |E        |Modulo de Young      |dato|
    +---------+---------------------+----+

Condiciones de la viga
----------------------

Condiciones de borde clásicas con una rótula restringida elasticamente en
x = cl

    +--------------------+-----+-----+-----+-----+
    |Condiciones de borde|p1(1)|q1(1)|p1(2)|q1(2)|
    +====================+=====+=====+=====+=====+
    |        S-S         |  1  |  x  |  1  | x-1 |
    +--------------------+-----+-----+-----+-----+
    |        S-F         |  1  |  x  |  1  |  1  |
    +--------------------+-----+-----+-----+-----+
    |        F-F         |  1  |  1  |  1  |  1  |
    +--------------------+-----+-----+-----+-----+
    |        C-C         |  x  |  x  | x-1 | x-1 |
    +--------------------+-----+-----+-----+-----+
    |        C-S         |  x  |  x  |  1  | x-1 |
    +--------------------+-----+-----+-----+-----+
    |        C-f         |  x  |  x  |  1  |  1  |
    +--------------------+-----+-----+-----+-----+
    |Para describir las condiciones de contorno  |                                        
    |se introduce la siguiente terminologia:     |                    
    |S: Simplemente apoyado, C: empotrado,       |
    |F: Libre                                    |
    +--------------------------------------------+

Los polinomios de orden superior se obtienen como::

    pi(k) = p1(k)x^i-1, i=1,2,..,M k=1,2
    qi(k) = q1(k)x^j-1, j=1,2,..,N k=1,2
