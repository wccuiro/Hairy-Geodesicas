import os

try:
    os.mkdir("./hola")
    print("Directorio creado")
except FileExistsError:
    print("Directorio ya creado")
