# TER_DATABASES
Aqui se guardan las bases de datos generales del TER, Son Creadas con el programa SavvyCAN u otros programas que permitan la edición de .dbcs
**La idea es que este repositorio sirva como unica fuente de verdad sobre las señales del coche, por favor acuerdate de hacer un commit cada vez que hagas un cambio.**

## Automatización
Todos los dbcs que se suban a este repositorio serán convertidos a sus correspondientes .h y .c para su uso
en los programas, mediante github actions, se sanean todos los dbcs eliminando las opciones propias de cada editor y se procesan para generar los archivos mediante la librería cantools de python.

# Inclusion como Submodulo
En tus proyectos de código es recomendable incluir este repositorio como submodulo para tener siempre la ultima versión del dbc a la hora de compilar, ejecutando:
```{bash}
git submodule add https://github.com/Tecnun-eRacing/TeR_DATABASES.git
git submodule update
```
Así contarás con la ultima versión del dbc siempre en tu código. No olvides incluir la carpeta DBC en tu IDE para que pueda encontrar los .h/.c.

# Clonado de proyectos que incluyen el submódulo
Cuando bajamos un repositorio que incluye un submódulo necesitamos decirle a git que baje también el submodulo para esto podemos hacerlo desde el principio con:
```{bash}
git clone https://mirepo.org/sustituyeaquiturepo.git --recursive #Esto baja el repositorio y todos los submodulos
```
Si bien ya hemos bajado el repositorio, tendremos que inicializarlo y bajarlo con:

```{bash}
git submodule init #inicializa los submodulos
git submodule update #los descarga
```
