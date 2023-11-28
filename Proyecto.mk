Documento.pdf: Documento.tex Presentacion.tex Matriz0.dat Matriz100.dat Temperatura.jpg TemperaturaMedia.jpg Volumen.jpg DistTemp.jpg DistTempFinal.jpg Matriz0Copia.dat Matriz100Copia.dat TemperaturaCopia.jpg TemperaturaMediaCopia.jpg VolumenCopia.jpg DistTempCopia.jpg DistTempFinalCopia.jpg
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
Matriz0Copia.dat: ProyectoHuevoCopia.cpp DatosInicialesCopia.dat
	g++ -std=c++11 ProyectoHuevoCopia.cpp -o ProyectoHuevoCopia.exe
	./ProyectoHuevoCopia.exe
Matriz100Copia.dat: ProyectoHuevoCopia.cpp DatosInicialesCopia.dat
	g++ -std=c++11 ProyectoHuevoCopia.cpp -o ProyectoHuevoCopia.exe
	./ProyectoHuevoCopia.exe
TemperaturaCopia.jpg: GraficarHuevoCopia.py Matriz0Copia.dat
	python3 GraficarHuevoCopia.py
TemperaturaMediaCopia.dat: ProyectoHuevoCopia.cpp DatosInicialesCopia.dat
	g++ -std=c++11 ProyectoHuevoCopia.cpp -o ProyectoHuevoCopia.exe
	./ProyectoHuevoCopia.exe
TemperaturaMediaCopia.jpg: TempMediaCopia.py TemperaturaMediaCopia.dat
	python3 TempMediaCopia.py
VolumenCopia.jpg: VolumenCopia.py VolumenCopia.dat
	python3 VolumenCopia.py
VolumenCopia.dat: ProyectoHuevoCopia.cpp DatosInicialesCopia.dat
	g++ -std=c++11 ProyectoHuevoCopia.cpp -o ProyectoHuevoCopia.exe
	./ProyectoHuevoCopia.exe
DistTempCopia.jpg: DistTempCopia.py Matriz0Copia.dat
	python3 DistTempCopia.py
DistTempFinalCopia.jpg: DistTempFinalCopia.py Matriz100Copia.dat
	python3 DistTempFinalCopia.py
