# help-desk-104

Dockerfile and support files
for reproducing https://github.com/csdms/help-desk/issues/104,
built on condaforge/miniforge3 (Debian GNU/Linux).

From the repo directory, build with:
```
$ docker build --tag help-desk-104 .
```

Launch a container a view its output with:
```
$ docker run --rm help-desk-104
The 2D Heat Equation
```
