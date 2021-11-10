# Evolución Tidal de Sistemas Planetarios

## Descripción
Este código fue desarrollado y utilizado para la materia *Evolución Tidal de Sistemas Planetarios*, de la FaMAF (UNC). 🇦🇷

Si principal función es la evolución de un sistema orbital de 3 cuerpos, incluyendo efectos **tidales**, **gravitatorios** y **relativistas**; utilizando un sistema de referencia de _Jacobi_.

### Parámetros
Al comienzo del archivo [main.F90](./main.F90#L9#L53) se encuentra un bloque para intorucir los parámetros de los 3 cuerpos, y de la integración.

Los parametros iniciales de los cuerpos utilizados son:
- m: Masa _m_. [M<sub>⊙</sub>]
- a: Semieje mayor _a_. [UA]
- e: Eccentricidad _e_.
- s: Spin Ω<sub>s</sub>. [Rad Día⁻¹]
- o: Oblicuidad ε. [Rad]
- vp: Longitud del pericentro . [Rad]
- r: Radio _R_. [UA]
- alp: Alpha de momento de inercia α (0.4 ⇒ _esfera_).
- Q: Parámetro tidal (_Tidal Quality Factor_).

Los parametros iniciales de la integración son:
- t0: Tiempo inicial _t_<sub>0</sub>. [Día]
- dt: Paso de tiempo mínimo e inicial _dt_. [Día]
- tf: Tiempo final _t_<sub>f</sub>. [Día]
- n_points: Cantiad aproximada de datos de salida.
- beta: Tasa de aprendizaje, en caso que se utilize un integrador de paso adaptativo β.
- e_tol: Error absoluto máximo, en caso que se utilize un integrador de paso adaptativo ϵ<sub>tol</sub> (≡ |<b>y</b><sub>real</sub> - <b>y</b><sub>pred</sub>|).
- filename: Nombre del archivo de salida.

### Implementaciones modificables por usuario
Se puede modificar el **integrador utilizado**, como también las **fuerzas involucradas**.
#### **Integrador**
Se debe utilar uno de los 4 métodos disponibles: *integ_caller*, *implicit_caller*, *rec_rk_adap*, o *embedded_caller*. 

Esto se realiza en la líneas _155-159_ del archivo [main.F90](./main.F90#L155#L159). 
```fortran
!!! Execute an integration method (uncomment/edit one of theese)
! call implicit_caller (t, y, dt, dydtidall, euler_centred, max_iter, e_tol, ynew)
! call integ_caller (t, y, dt, dydtidall, rungek6, ynew)
! call rk_half_step_caller (t, y, dt_adap, dydtidall, rungek6, 6, e_tol, beta, dt_min, dt, ynew)
call embedded_caller (t, y, dt_adap, dydtidall, Bulirsch_Stoer, e_tol, beta, dt_min, dt, ynew)
```

- *implicit_caller*: Se puede modificar el integrador (```euler_centred```) por cualquier otro que sea implícito (no embebido), del archivo [integrators.F90](./integrators.F90#L614#L644). El paso _dt_ será constante, el igual al valor mínimo _dt_<sub>min</sub>. 

- *integ_caller*: Igual a *implicit_caller*, pero solo se pueden introducir integradores _NDimensionales_ no implícitos (ni embebidos) (ver [integrators.F90](./integrators.F90#L646#L865)).

- *rk_half_step_caller*: Igual a *integ_caller*, pero al cambiar el integrador, también se debe introducir su orden (ej. ```rungek6``` -> O(6)). En este caso el integrador intentará utilizar un paso de tiempo adaptativo _dt_<sub>adap</sub> adecuado, según la toleranca de error ϵ<sub>tol</sub> introducida (ver [integrators.F90](./integrators.F90#L1167#L1209)). En este caso el error calculado será ϵ<sub>calc</sub> (≡ |<b>y</b><sub>real</sub> - <b>y</b><sub>pred</sub>|)/(2<sup>ord</sup> - 1).

- *embedded_caller* (**recomendado**): Similar a *rk_half_step_caller* (debido a que también calcula un paso de tiempo óptimo), pero solo se pueden introducir integradores _NDimensionales_ embebidos (ver [integrators.F90](./integrators.F90#L867#L1165)), incluyendo al integrador *Bulirsch_Stoer* (**recomendado** si se introduce un ϵ<sub>tol</sub> < 10<sup>-8</sup>, ver [bstoer.F90](./bstoer.F90#L18#44)).


#### **Fuerzas involucradas**
Esto se realiza en la líneas _356-360_ del archivo [tidall.F90](./tidall.F90#L356#L360).
``` fortran
call dydtidal (y0, y1, 1, y01t, y10t)
call dydtidal (y0, y2, 2, y02t, y20t)
call dydtgrav (y1, 1, y2, 2, y12g, y21g)
call dydtrela (y1, 1, y1r)
call dydtrela (y2, 2, y2r)
```
En caso de querer despreciar (o simplemente no calcular) alguna fuerza en especial (*tidal*, *gravitatoria* o *efectos relativistas*), se debe comentar la línea correspondiente.

## Requerimientos
Compilar y ejecutar:
- ```gfortran```

Plotear:
- ```python```
  - ```numpy```
  - ```matplotlib```
  - ```pandas```

Recompilar y ejecutar (en caso de realizar modificaciones al código que cambien las dependencias):
- ```pip``` 
- ```python```
  - [fortdepend](https://github.com/ZedThree/fort_depend.py) (paquete de python utilzado para generar nuevo archivo de dependencias)

## Modo de uso
Si se usa el código tal y como está, solo es necesario compilarlo y ejecutarlo.
En caso de realizarle modificaciones que alteren las dependencias, es necesario instalar [fortdepend](https://github.com/ZedThree/fort_depend.py).
### Crear y ejecutar
```console
$ make
$ ./tidal
```
En caso de que no exista el archivo [tidal.dep](./tidal.dep), es necesario tener instalado [fortdepend](https://github.com/ZedThree/fort_depend.py) para generarlo.
### Instalar _fortdepend_ y generar nuevo archivo de dependencias
```console
$ make install
```
### Borrar archivos de compilación para luego recompilar
```console
$ make clean
```
En este caso, también se borra el archivo de dependencias [tidal.dep](./tidal.dep), por lo que luego será necesario regenerarlo.

### Plotear
Se adjunta un pequeño código [plot.py](./plot.py) en ```python``` para realizar algunos plots resultantes.
```console
$ python plot.py
```

## Licencia
Distribucído bajo la Licencia MIT.  (Archivo [LICENSE](./LICENSE))

## Autor
Emmanuel Gianuzzi - egianuzzi@mi.unc.edu.ar
