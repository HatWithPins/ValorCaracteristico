# ValorCarateristico

:es: \
Implementación del método QR y de la potencia para obtener los autovalores de una matriz típica de EDP élípticas. Usa CMake para la compilación.

Para compilar en la línea de comandos:

`cmake CMakeList.txt` \
`cmake --build . --target ValorCarateristico`

Si se usa un IDE tipo Visual Studio, abre la carpeta como un proyecto y compila como cualquier otro proyecto. Asegúrate de que tienes la extensión para CMake.

Una vez se ha compilado, hay que pasarle el tamaño de la matriz N y el parámetro de control lambda.

`ValorCaracteristico N=10 lambda=0.20`

Si no recibe los argumentos que espera o si tienen algún problema, devuelve error.

## Estructura de los archivos.
- ValorCaracteristico contiene el main y las comprobaciones de los argumentos.
- MetodoQR contiene la clase que implementa el método QR.
- MetodoPotencias contiene la clase que implementa el método de la potencia.

:uk: \
Implementation of QR and power method to obtain eigenvalues for a typical elliptical PDE. Builds using cmake.

To build using command line:

`cmake CMakeList.txt` \
`cmake --build . --target ValorCarateristico`

In case you use and IDE like Visual Studio, open the folder as a project and build it like any other project. Ensure that your IDE has the extensions for CMake.

Once is built, it needs the matrix size N and control parameter lambda.

`ValorCaracteristico N=10 lambda=0.20`

If the program doesn't receive the expected arguments or if they have any issue, program will return an error.

## Folder structure.
- ValorCaracteristico contains the main function and is in charge of checking parsed arguments.
- MetodoQR contains the class that implements QR method.
- MetodoPotencias contains the class that implements power method.