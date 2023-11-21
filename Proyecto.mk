Documento.pdf: Documento.tex Presentacion.tex Matriz.dat Temperatura.jpg
	pdflatex Documento.tex
Presentacion.pdf: Presentacion.tex Temperatura.jpg
	pdflatex -synctex=1 -interaction=nonstopmode "Presentacion".tex
Matriz.dat: ProyectoHuevo.cpp
	g++ -std=c++11 ProyectoHuevo.cpp -o ProyectoHuevo.exe
	./ProyectoHuevo.exe
Temperatura.jpg: GraficarHuevo.py Matriz.dat
	python3 GraficarHuevo.py
