G2DFEAS
Version 3
General two dimensional Finite Element Analysis Script
Requiere variables cargados
Usar unidades CONSISTENTES: [N,mm,Mpa,mm^2,mm^4] [N,m,Pa,m^2,m^4] etc.
  nod=[x1 y1;x2 y2... xn yn]
  elenod=[N1 N2;N2 N3 ] conexiones entre nodos (elementos). Importante
  tener bien definido el orden numerado de los elementos y los nodos!
  Despues los resultados se basan en el orden que se tiene en elenod.
  eletype=[1 2 2 1 3 4 1 2... 2 1 1];
      1=barra
      2=viga
      3=viga rotulada en su nodo de inicio
      4=viga rotulada en su nodo final o destino
      11,22,33,44= elemento simetrico (se divide su rigidez por dos)
  apoyos_simples=[N1 N5 ... N8] Nodos apoyados en x e y
  empotramientos=[] Nodos empotrados
  R=[fx1 fy1 mz1 fx2 fy2 mz2 ... ] vector de fuerzas sobre los nodos. Usar función fuerzapuntual() para aplicar cargas.

      MATERIALES:
  Ee=[200e9 210e9 ... 200e9] Modulos de young para cada elemento, se sigue el
  orden de elenod.
  Similarmente, se tiene que declarar Ae (area), Ie (2do momento de
  inercia, he (altura/radio de viga, puede ser arbitrario para
  barras), be (ancho/diametro para viga), Sye (tension de fluencia)
  safetyfactor=4; factor seguridad deseado

      OPCIONES AVANZADAS (opcional):
  CB=~~[1 1 0 1 1 ... 0  0  1]   Condiciones de borde (TRUE/1 anula grado de libertad)
  Todos los nodos tienen 3 grados de liberad. Tomar esto en cuenta al
  aplicar condiciones de borde. Las rotulas generan un grado de libertad
  extra al final de CB, siguiendo el orden de elenod.
  
  graficar=true;  grafica estructura modelada
  vigasinteresantes=[4:9] grafica esfuerzos sobre elementos 4 a 9 si
  graficar=true.
  showresults=false; No muestra resultados.
      RESULTADOS:   
  Muestra inmediatamente los resultados de tensiones sobre las barras.
  dangerzone:
  primera fila numera los elementos, segunda fila muestra el calculo:
  (tension/tensionadm)*factorseguridad. Cuanto mas alto, mas
  comprometido. En principio se desea que no exceda la unidad.
  Segunda fila muestra peligro de pandeo para barras articuladas en los extremos.
  
      OUTPUT :
  forzas: Devuelve estructura (fuerzas sobre cada elemento)
  loskrotados: estructura con matriz rigidez locales.
  
  Patricio Whittingslow CC-BY-NC  2018 
