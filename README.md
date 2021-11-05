# Evolución Tidal de Sistemas Planetarios

## Descripción
Este código fue desarrollado y utilizado para la materia *Evolución Tidal de Sistemas Planetarios*, de la FaMAF (UNC). 🇦🇷

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
  - ```fortdepend``` (paquete de python utilzado para generar nuevo archivo de dependencias)

## Modo de uso
Si se usa el código tal y como está, solo es necesario compilarlo y ejecutarlo.
En caso de realizarle modificaciones que alteren las dependencias, es necesario instalar _fortdepend_
### Crear y ejecutar
```console
$ make
$ ./tidal
```
En caso de que no exista el archivo _tidal.dep_, es necesario tener instalado _fortdepend_ para generarlo.
### Instalar _fortdepend_ y generar nuevo archivo de dependencias
```console
$ make install
```
### Borrar archivos de compilación para luego recompilar
```console
$ make clean
```
En este caso, también se borra el archivo de dependencias _tidal.dep_, por lo que luego será necesario regenerarlo.

### Plotear
Se adjunta un pequeño código en ```python``` para realizar algunos plots resultantes.
```console
$ python plot.py
```

## Licencia
Distribucído bajo la Licencia MIT.  (Archivo `LICENSE`)

## Autor
Emmanuel Gianuzzi - egianuzzi@mi.unc.edu.ar
