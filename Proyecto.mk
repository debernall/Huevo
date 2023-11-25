Documento.pdf: Documento.tex Presentacion.tex Matriz0.dat Matriz100.dat Temperatura.jpg TemperaturaMedia.jpg Volumen.jpg DistTemp.jpg DistTempFinal.jpg
	pdflatex Documento.tex
Presentacion.pdf: Presentacion.tex Temperatura.jpg
	pdflatex -synctex=1 -interaction=nonstopmode "Presentacion".tex
Matriz0.dat: ProyectoHuevo.cpp DatosIniciales.dat
	g++ -std=c++11 ProyectoHuevo.cpp -o ProyectoHuevo.exe
	./ProyectoHuevo.exe
Matriz100.dat: ProyectoHuevo.cpp DatosIniciales.dat
	g++ -std=c++11 ProyectoHuevo.cpp -o ProyectoHuevo.exe
	./ProyectoHuevo.exe
Temperatura.jpg: GraficarHuevo.py Matriz0.dat
	python3 GraficarHuevo.py
TemperaturaMedia.dat: ProyectoHuevo.cpp DatosIniciales.dat
	g++ -std=c++11 ProyectoHuevo.cpp -o ProyectoHuevo.exe
	./ProyectoHuevo.exe
TemperaturaMedia.jpg: TempMedia.py TemperaturaMedia.dat
	python3 TempMedia.py
Volumen.jpg: Volumen.py Volumen.dat
	python3 Volumen.py
Volumen.dat: ProyectoHuevo.cpp DatosIniciales.dat
	g++ -std=c++11 ProyectoHuevo.cpp -o ProyectoHuevo.exe
	./ProyectoHuevo.exe
DistTemp.jpg: DistTemp.py Matriz0.dat
	python3 DistTemp.py
DistTempFinal.jpg: DistTempFinal.py Matriz100.dat
	python3 DistTempFinal.py
