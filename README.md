# Evolución Tidal de Sistemas Planetarios

## Descripción
Este código fue desarrollado y utilizado para la materia *Evolución Tidal de Sistemas Planetarios*, de la FaMAF (UNC). 🇦🇷

Si principal función es la evolución de un sistema orbital de 3 cuerpos, incluyendo efectos **tidales**, **gravitatorios** y **relativistas**; utilizando un sistema de referencia de _Jacobi_.

### Parámetros
Al comienzo del archivo [main.F90](./main.F90#L9#L54) se encuentra un bloque para intorucir los parámetros de los 3 cuerpos, y de la integración.

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
- dt_min: Paso de tiempo para integradores de paso fijo, o mínimo _dt_<sub>min</sub>. [Día]
- tf: Tiempo final _t_<sub>f</sub>. [Día]
- beta: Tasa de aprendizaje, en caso que se utilize un integrador de paso adaptativo β.
- e_tol: Error absoluto máximo, en caso que se utilize un integrador de paso adaptativo ϵ<sub>tol</sub> (≡ |<b>y</b><sub>real</sub> - <b>y</b><sub>pred</sub>|).
- n_points: Cantiad aproximada de datos de salida.
- logsp: Si la salida es log- (_True_) o equi- (_False_) espaciada.
- filename: Nombre del archivo de salida.

### Implementaciones modificables por usuario
Se puede modificar el **integrador utilizado**, como también las **fuerzas involucradas**.
#### **Integrador**
Se debe utilar uno de los 4 métodos disponibles: *integ_caller*, *rk_half_step_caller*, *embedded_caller*, o *BStoer_caller*. 

Esto se realiza en la líneas _149-153_ del archivo [main.F90](./main.F90#L149#L153). 
```fortran
    !!! Execute an integration method (uncomment/edit one of theese)
    ! call integ_caller (t, y, dt_adap, dydtidall, Runge_Kutta4, dt, ynew)
    ! call rk_half_step_caller (t, y, dt_adap, dydtidall, Runge_Kutta5, 5, e_tol, beta, dt_min, dt, ynew)
    ! call embedded_caller (t, y, dt_adap, dydtidall, Fehlberg4_5, e_tol, beta, dt_min, dt, ynew)
    call BStoer_caller (t, y, dt_adap, dydtidall, e_tol, dt_min, dt, ynew)
```

- *integ_caller*: Se puede modificar el integrador (```Runge_Kutta4```) por cualquier otro que no sea embebido, entre los listados en el archivo [integrators.F90](./integrators.F90#L260#L628). El paso _dt_<sub>adap</sub> será constante, el igual al valor mínimo _dt_<sub>min</sub>.

- *rk_half_step_caller*: Igual a *integ_caller*, pero al cambiar el integrador, también se debe introducir su orden (ej. ```Runge_Kutta5``` -> O(5)). En este caso el integrador intentará utilizar un paso de tiempo adaptativo _dt_<sub>adap</sub> adecuado, según la toleranca de error ϵ<sub>tol</sub> introducida (ver [integrators.F90](./integrators.F90#L934#L983)). En este caso el error calculado será ϵ<sub>calc</sub> (≡ |<b>y</b><sub>aux</sub> - <b>y</b><sub>pred</sub>|)/(2<sup>ord</sup> - 1).

- *embedded_caller*: Similar a *rk_half_step_caller* (debido a que también calcula un paso de tiempo óptimo), pero solo se pueden introducir integradores embebidos (ver [integrators.F90](./integrators.F90#L630#L932)). En este caso el error calculado será ϵ<sub>calc</sub> (≡ |<b>y</b><sub>aux</sub> - <b>y</b><sub>pred</sub>|).

- *BStoer_caller*:, Se utiliza el integrador *Bulirsch_Stoer* (**recomendado** si se introduce un ϵ<sub>tol</sub> < 10<sup>-8</sup>, ver [bstoer.F90](./bstoer.F90#L18#L62)).


#### **Fuerzas involucradas**
Esto se realiza en la líneas _358-362_ del archivo [tidall.F90](./tidall.F90#L358#L362).
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
