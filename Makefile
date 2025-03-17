CC=mpicc

pruebas: pruebas.c
	$(CC) pruebas.c -o pruebas.exe

clean:
	rm *.exe
