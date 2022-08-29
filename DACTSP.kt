/**
 * Universidad Simón Bolívar. Proyecto 1
 *
 * Descripción: Programa que encuentra solución al
 * problema del TSP euleriano. El objetivo es el
 * poder encontrar un tour compuesto de lados (par de ciudades), 
 * que comience y finalice en la misma ciudad, con la menor 
 * distancia posible y sin pasar más de una vez por una ciudad.
 * 
 * @author Ka Fung <18-10492@usb.ve>
 */

import kotlin.math.*
import kotlin.Double.Companion.POSITIVE_INFINITY
import java.io.File

/******************* INICIO DE QUICKSORT (INTROSORT) ********************/
/**
 * Intercambia dos elementos de un arreglo A (arreglo de pares).
 * Precondición: A.size >= 0 && 0 <= i, j < A.size
 * Postcondición: A.size >= 0 && 0 <= i, j < A.size
 *
 * @param[A] Arreglo de pares de doubles
 * @param[i] Índice de un elemento del arreglo a intercambiar
 * @param[j] Índice de un elemento del arreglo a intercambiar
 */
fun exchange(A: Array<Pair<Double, Double>>, i: Int, j: Int) {
    val temp = A[i];
    A[i] = A[j];
    A[j] = temp;
}

/**
 * Verifica si cumple la relación de orden. Si el eje es X, se
 * compara el primer elemento del par1 con el del par2, si son iguales
 * se compara el segundo elemento del par1 con el del par2.
 * Precondición: par1.size == 2 && par2.size == 2 && (eje == "X" || eje == "Y")
 * Postcondición: (a < a0 || (a == a0 && b < b0))
 *
 * @param[par1] Par de double
 * @param[par2] Par de double
 * @param[eje] Eje de coordenada (X o Y)
 */
fun relacionDeOrden(par1: Pair<Double, Double>, par2: Pair<Double, Double>, eje: String): Boolean {
    val a = if (eje == "X") par1.first else par1.second
    val a0 = if (eje == "X") par2.first else par2.second
    val b = if (eje == "X") par1.second else par1.first
    val b0 = if (eje == "X") par2.second else par2.first
    return (a < a0 || (a == a0 && b < b0))
}

/**
 * Ordena un arreglo A con el algoritmo InsertionSort
 * Precondición: A.size >= 0 && 0 <= f, b <= A.size && (eje == "X" || eje == "Y")
 * Postcondición: Arreglo ordenado de A[f..(b - 1)]
 *
 * @param[A] Arreglo de pares de doubles a ordenar
 * @param[f] Índice de inicio a ordenar
 * @param[b] Índice de final a ordenar
 * @param[eje] Eje de coordenada con respecto se quiere ordenar (X o Y)
 */
fun insertionSort(A: Array<Pair<Double, Double>>, f: Int, b: Int, eje: String) {
    for (j in (f + 1) until b) {
        val key = A[j]
        var i = j - 1
        while (i >= 0 && relacionDeOrden(key, A[i], eje)) {
            A[i + 1] = A[i]
            i--
        }
        A[i + 1] = key
    }
}

/**
 * Retorna el hijo izquierdo o derecho de un nodo de índice i
 * Precondición: i >= 0
 * Postcondición: El valor a retornar debe ser mayor que i
 *
 * @param[i] Índice del nodo cual se quiere conseguir su hijo.
 *
 * @return Un entero que es el índice del hijo de un nodo i
 */
fun left(i: Int): Int = (2 * i) + 1
fun right(i: Int): Int = (2 * i) + 2

/**
 * Mantiene las propiedades de un Heap en un arreglo A, desde el nodo i.
 * Precondición: (0 <= i <= A.size / 2) && (eje == "X" || eje == "Y") 
 * Postcondición: Para todo nodo i, A[i] >= left(i) && A[i] >= right(i)
 *
 * @param[A] Arreglo de pares de doubles a ordenar
 * @param[i] Índice (entero) del nodo a aplicar el maxheap
 * @param[heapsize] (entero) Tamaño del heap
 * @param[dif] Offset (entero) de desde dónde ordenar
 * @param[eje] Eje de coordenada con respecto se quiere ordenar (X o Y)
 */
fun maxHeapify(A: Array<Pair<Double, Double>>, i: Int, heapsize: Int, dif: Int, eje: String) {
    val l = left(i) - dif
    val r = right(i) - dif
    var largest = i
    if (l < heapsize && relacionDeOrden(A[i], A[l], eje)) {
        largest = l
    }
    if (r < heapsize && relacionDeOrden(A[largest], A[r], eje)) {
        largest = r
    }
    if (largest != i) {
        exchange(A, i, largest)
        maxHeapify(A, largest, heapsize, dif, eje)
    }
}

/**
 * Convierte un subarreglo del arreglo de A en un Max-heap. 
 * Precondición: 0 <= begin, end < A.size && (eje == "X" || eje == "Y")
 * Postcondición: A[begin .. end] debe cumplir con las propiedades de un Max-Heap
 *
 * @param[A] Arreglo de pares que se quiere construir el Max-Heap
 * @param[begin] Índice (entero) del comienzo a construir el Max-Heap
 * @param[end] Índice (entero) del final a construir el Max-Heap
 * @param[eje] Eje de coordenada con respecto se quiere ordenar (X o Y)
 */
fun buildMaxHeap(A: Array<Pair<Double, Double>>, begin: Int, end: Int, eje: String) {
    for (i in (begin + ((end - begin) / 2) - 1) downTo begin) {
        maxHeapify(A, i, end, begin, eje)
    }
}

/** 
 * Ordena un subarreglo del arreglo de A con el algoritmo HeapSort.
 * Precondición: 0 <= begin, end < A.size && (eje == "X" || eje == "Y")
 * Postcondición: Para todo begin <= i < end: A[i] <= A[i + 1] (arreglo ordenado A de begin a end)
 *
 * @param[A] Arreglo de pares al que se quiere ordenar
 * @param[begin] Índice (entero) del comienzo de A a ordenar
 * @param[end] Índice (entero) del final de A a ordenar
 * @param[eje] Eje de coordenada con respecto se quiere ordenar (X o Y)
 */
fun heapSort(A: Array<Pair<Double, Double>>, begin: Int, end: Int, eje: String) {
    buildMaxHeap(A, begin, end, eje)
    var heapsize = end - begin - 1
    for (i in heapsize downTo 1) {
        exchange(A, begin, begin + i)
        maxHeapify(A, begin, begin + i, begin, eje)
    }
}

/**
 * Retorna un punto de partición del arreglo A.
 * Precondición: 0 <= p, r < A.size && (eje == "X" || eje == "Y") && (x in A)
 * Postcondición: 0 <= i < A.size
 *
 * @param[A] Arreglo de pares de doubles al que se quiere particionar
 * @param[p] Índice (entero) del comienzo de A a particionar
 * @param[q] Índice (entero) del final de A a particionar
 * @param[x] pivote (entero) del arreglo
 * @param[eje] Eje de coordenada con respecto se quiere ordenar (X o Y)
 *
 * @return índice del pivote, tras realizar la partición
 */
fun partitionHoare(A: Array<Pair<Double, Double>>, p: Int, r: Int, eje: String): Int {
    var i = p - 1
    var j = r
    val pivot = A[r - 1] 
    while (true) {
        do {
            j--
        } while (relacionDeOrden(pivot, A[j], eje));
        
        do {
            i++
        } while (relacionDeOrden(A[i], pivot, eje));

        if (i < j) {
            exchange(A, i, j)
        } else {
            return i
        }
    }
}

/**
 * Ordena un arreglo A con IntroSort y HeapSort.
 * Precondición: (Para todo f <= i < j < b : A[i] <= A[j] && i - j > size_threshold)
 * && (eje == "X" || eje == "Y")
 * Postcondición: Para todo f <= i <= b : A[i] <= A[i + 1] (arreglo ordenado de A de begin a end)
 *
 * @param[A] arreglo a ser ordenado, f indice del primer elemento del subarreglo de A
 * b indice del último elemento del subarreglo A, depth_limit valor límite de recurrencias.
 */
fun introsortLoop(A: Array<Pair<Double, Double>>, f: Int, b: Int, depth_limit: Int, eje: String) {
    var bb = b
    var dl = depth_limit
    val sizeThreshold = 32
    while (bb - f > sizeThreshold) {
        if (dl == 0) {
            heapSort(A, f, b, eje)
            return
        }
        dl--
        val p = partitionHoare(A, f, bb, eje)
        introsortLoop(A, p, bb, dl, eje)
        bb = p
    }
}

/**
 * Ordena un arreglo A de manera ascendente con el algoritmo QuickSort 
 * presentado por David Musser (IntroSort).
 * Precondición: (0 <= f, b <= A.size) && (eje == "X" || eje == "Y")
 * Postcondición: Para todo f <= i < b : A[i] <= A[i + 1] (arreglo ordenado A de f a (b - 1))
 *
 * @param[A] Arreglo de pares de doubles al que se quiere particionar
 * @param[f] Índice (entero) del comienzo de A a ordenar
 * @param[b] Índice (entero) del final de A a ordenar
 * @param[eje] Eje de coordenada con respecto se quiere ordenar (X o Y)
 *
 * @return índice del pivote, tras realizar la partición
 */
fun introsort(A: Array<Pair<Double, Double>>, f: Int, b: Int, eje: String) {
    introsortLoop(A, f, b, 2 * ( floor( log2( (b - f).toDouble() ) ).toInt() ), eje)
    insertionSort(A, f, b, eje)
}
/******************* FIN DE QUICKSORT (INTROSORT) ********************/

/**
 * Crea un ciclo de una ciudad.
 * Precondición: P.size == 1
 * Postcondición: Un arreglo de un par de dos pares iguales.
 *
 * @param[P] Arreglo de pares de doubles (ciudades)
 *
 * @return Arreglo de un par de pares de doubles (ciudad).
 */
fun cicloUnaCiudad(
    P: Array<Pair<Double, Double>>
): Array<Pair<Pair<Double, Double>, Pair<Double, Double>>> {
    return arrayOf(Pair(P[0], P[0]))
}

/**
 * Crea un ciclo de dos ciudadades.
 * Precondición: P.size == 2
 * Postcondición: Un arreglo de dos pares de dos pares.
 *
 * @param[P] Arreglo de pares de doubles (ciudades)
 *
 * @return Arreglo de dos pares de pares de doubles (ciudad).
 */
fun cicloDosCiudades(
    P: Array<Pair<Double, Double>>
): Array<Pair<Pair<Double, Double>, Pair<Double, Double>>> {
    return arrayOf(Pair(P[0], P[1]), Pair(P[0], P[1]))
}

/**
 * Crea un ciclo de tres ciudadades.
 * Precondición: P.size == 3
 * Postcondición: Un arreglo de tres pares de dos pares.
 *
 * @param[P] Arreglo de pares de doubles (ciudades)
 *
 * @return Arreglo de tres pares de pares de doubles (ciudad).
 */
fun cicloTresCiudades(
    P: Array<Pair<Double, Double>>
): Array<Pair<Pair<Double, Double>, Pair<Double, Double>>> {
    return arrayOf(Pair(P[0], P[1]), Pair(P[0], P[2]), Pair(P[1], P[2]))
}

/**
 * Calcula las coordenadas de un rectángulo dado un arreglo de ciudades
 * Precondición: P.size > 0
 * Postcondición: Cada ciudad (par) debe estar dentro del rango de las coords del rectángulo
 * 
 * @param[P] Arreglo de pares de doubles (ciudades).
 *
 * @return Arreglo de pares de doubles (coordenadas del rectángulo)
 */
fun rectanguloCoords(P: Array<Pair<Double, Double>>): Array<Pair<Double, Double>> {
    var (menorX, mayorX) = Pair(P[0].first, P[0].first)
    var (menorY, mayorY) = Pair(P[0].second, P[0].second)
    for (coord in P) {
        if (coord.first < menorX) menorX = coord.first
        if (coord.first > mayorX) mayorX = coord.first
        if (coord.second < menorY) menorY = coord.second
        if (coord.second > mayorY) mayorY = coord.second  
    }
    val coords = arrayOf(Pair(menorX, menorY), Pair(mayorX, menorY), Pair(menorX, mayorY), Pair(mayorX, mayorY))
    return coords
}

/**
 * Calcula las dimensiones de un rectángulo, dado sus coordenadas.
 * Precondición: rectangulo.size == 4
 * Postcondición: Pair(mayorX - menorX, mayorY - menorY)
 *
 * @param[rectangulo] Arreglo de pares de doubles (coordenadas)
 *
 * @return Par del tamaño de lado horizontal (X) y vertical (Y) del rectángulo
 */
fun rectanguloDim(rectangulo: Array<Pair<Double, Double>>): Pair<Double, Double> {
    var (menorX, mayorX) = Pair(rectangulo[0].first, rectangulo[1].first)
    var (menorY, mayorY) = Pair(rectangulo[0].second, rectangulo[2].second)
    return Pair(mayorX - menorX, mayorY - menorY)
}

/**
 * Ordena el arreglo de ciudades según la coordenada X o Y, y 
 * consigue el punto de corte (ciudad) del rectángulo.
 * Precondición: P.size >= 0
 * Postcondición: P[0] <= P[pos] <= P[P.size - 1]
 *
 * @param[P] Arreglo de pares de doubles (ciudades).
 * @param[ejeDeCorte] String que indica el eje de corte (X o Y).
 *
 * @return Par de doubles que es la coordenada de la ciudad que se encuentra más próxima al medio
 */
fun obtenerPuntoDeCorte(P: Array<Pair<Double, Double>>, ejeDeCorte: String): Pair<Double, Double> {
    val n = P.size
    val pos = ceil(n.toDouble() / 2.0).toInt() - 1
    introsort(P, 0, n, ejeDeCorte) // Ordena las ciudades de P con respecto a X o Y
    return P[pos]
}

/**
 * Consigue el punto de corte a la mitad del rectángulo según el eje X o Y
 * Precondición: (rectangulo.size == 4) && (ejeDeCorte == "X" || ejeDeCorte == "Y")
 * Postcondición: xmin <= puntoDeCorte.first && ymin <= puntoDeCorte.second
 *
 * @param[rectangulo] Arreglo de pares de doubles (coordenadas)
 * @param[ejeDeCorte] String del eje de corte (X o Y)
 *
 * @return Par de double que es la coordenada del punto de corte.
 */
fun obtenerPuntoDeCorteMitad(
    rectangulo: Array<Pair<Double, Double>>, 
    ejeDeCorte: String
): Pair<Double, Double> {
    val (xmin, ymin) = Pair(rectangulo[0].first, rectangulo[0].second)
    var puntoDeCorte: Pair<Double, Double>
    if (ejeDeCorte == "X") {
        // Obtener el valor del lado del rectángulo que es paralelo al eje X
        val xdim = rectangulo[1].first - rectangulo[0].first
        puntoDeCorte = Pair(xmin + (xdim / 2), ymin)
    } else {
        // Obtener el valor del lado del rectángulo que es paralelo al eje Y
        val ydim = rectangulo[2].second - rectangulo[0].second
        puntoDeCorte = Pair(xmin, ymin + (ydim / 2))
    }
    return puntoDeCorte
}

/**
 * Aplica el corte del rectangulo que contiene a todas las ciudades de P
 * trazando una recta perpendicular al ejeDeCorte (X o Y), retornando dos rectángulos.
 * Precondición: (ejeDeCorte == "X" || ejeDeCorte == "Y") && (rectangulo.size == 4)
 * Postcondición: Dos rectángulos tal que sus medidas no sobrepasen al rectangulo dado.
 *
 * @param[ejeDeCorte] String del eje de corte (X o Y)
 * @param[coord] Arreglo de pares de doubles que son el punto de corte
 * @param[rectangulo] Arreglo de pares de doubles (coordenadas)
 * @return Par de arreglos de pares de doubles que son las coordenadas de los rectángulos.
 */
fun aplicarCorte( 
    ejeDeCorte: String, 
    coord: Pair<Double, Double>, 
    rectangulo: Array<Pair<Double, Double>>
): Pair<Array<Pair<Double, Double>>, Array<Pair<Double, Double>>> {
    var rectanguloIzq: Array<Pair<Double, Double>>
    var rectanguloDer: Array<Pair<Double, Double>>
    val (xc, yc) = coord

    if (ejeDeCorte == "X") {
        // Recta perpendicular desde xc
        var medioInferior = Pair(xc, rectangulo[0].second)
        var medioSuperior = Pair(xc, rectangulo[2].second)
        rectanguloIzq = arrayOf(rectangulo[0], medioInferior, rectangulo[2], medioSuperior)
        medioInferior = Pair(xc + 0.000001, rectangulo[0].second)
        medioSuperior = Pair(xc + 0.000001, rectangulo[2].second)
        rectanguloDer = arrayOf(medioInferior, rectangulo[1], medioSuperior, rectangulo[3])
    } else {
        // Recta perpendicular desde yc
        var medioIzq = Pair(rectangulo[0].first, yc)
        var medioDer = Pair(rectangulo[1].first, yc)
        rectanguloIzq = arrayOf(rectangulo[0], rectangulo[1], medioIzq, medioDer)
        medioIzq = Pair(rectangulo[0].first, yc + 0.000001)
        medioDer = Pair(rectangulo[1].first, yc + 0.000001)
        rectanguloDer = arrayOf(medioIzq, medioDer, rectangulo[2], rectangulo[3])
    }
    return Pair(rectanguloIzq, rectanguloDer)
}

/**
 * Retorna las ciudades que se encuentran dentro de una rectángulo dado.
 * Precondición: (P.size >= 0) && (rectangulo.size == 4)
 * Postcondición: Todas las coordenadas de las ciudades a retornar deben estar dentro 
 * del rango de las coordenadas del rectángulo dado.
 *
 * @param[P] Arreglo de pares de doubles (ciudades).
 * @param[rectangulo] Arreglo de pares de doubles (coordenadas)
 *
 * @return Arreglo de pares de doubles que son las ciudades que se encuentran dentro del rectángulo
 */
fun obtenerPuntosRectangulo(
    P: Array<Pair<Double, Double>>, 
    rectangulo: Array<Pair<Double, Double>>
): Array<Pair<Double, Double>> {
    // Coordenadas x mínima, x máxima, y mínima, y máxima
    val (xMin, xMax) = Pair(rectangulo[0].first, rectangulo[1].first)
    val (yMin, yMax) = Pair(rectangulo[0].second, rectangulo[2].second)
    
    // Buscamos el index por donde comienza las ciudades del rectangulo
    val n = P.size
    var i = 0
    var j = n - 1
    while (i < n && !(xMin <= P[i].first && P[i].first <= xMax && 
        yMin <= P[i].second && P[i].second <= yMax)
    ) {
        i++
    }

    // Buscamos el index por donde finaliza las ciudades del rectangulo
    while (j >= 0 && !(xMin <= P[j].first && P[j].first <= xMax && 
        yMin <= P[j].second && P[j].second <= yMax)
    ) {
        j--
    }
    val ciudades = P.sliceArray(i..j)
    return ciudades
}

/** 
 * Calcula la distancia entre dos puntos i, j, redondeado a entero
 * con la fórmula de TSPLIB EUD_2D.
 * Precondición: i.size == 2 && j.size == 2
 * Postcondición: (sqrt((xd * xd) + (yd * yd)) + 0.5).toInt()
 *
 * @param[i] Par de doubles, coordenadas de un punto.
 * @param[j] Par de doubles, coordenadas de un punto.
 *
 * @return Distancia entre dos puntos i, j.
 */
fun distancia(i: Pair<Double, Double>, j: Pair<Double, Double>): Double {
    val xd = i.first - j.first
    val yd = i.second - j.second
    val dij = sqrt((xd * xd) + (yd * yd))
    return dij
}

/** 
 * Calcula la distancia ganada entre dos distancias viejas y dos nuevas distancias.
 * Precondición: dNEW1 >= 0 && dNEW2 >= 0 && dOLD1 >= 0 && dOLD2 >= 0
 * Postcondición: (dNEW1 + dNEW2) - (dOLD1 + dOLD2)
 *
 * @param[dOLD1] Número entero que representa una distancia.
 * @param[dOLD2] Número entero que representa una distancia.
 * @param[dNEW1] Número entero que representa una distancia.
 * @param[dNEW2] Número entero que representa una distancia.
 *
 * @return Distancia entre dos puntos i, j. 
 */
fun distanciaGanada(dOLD1: Double, dOLD2: Double, dNEW1: Double, dNEW2: Double): Double {
    return ((dNEW1 + dNEW2) - (dOLD1 + dOLD2))
}

/** 
 * Ordena un ciclo, tal que para cada lado y el siguiente lado compartan una ciudad.
 * Precondición: Ciclo.size > 0
 * Postcondición: Ciclo ordenado de manera que cumpla con la propiedad de un tour válido
 *
 * @param[Ciclo] Arreglo de pares de pares de doubles, arreglo de lados.
 */
fun selectionSortCiclo(Ciclo: Array<Pair<Pair<Double, Double>, Pair<Double, Double>>>) {
    val n = Ciclo.size
    for (i in 1 until (n - 1)) {
        var ladoCambiar = i      // Indice del actual lado a intercambiar
        for (j in i until n) {
            // Si el anterior lado comparte con un lado j, se intercambia el actual con ese lado
            if ((Ciclo[i - 1].first == Ciclo[j].first) || 
                (Ciclo[i - 1].first == Ciclo[j].second) ||
                (Ciclo[i - 1].second == Ciclo[j].first) || 
                (Ciclo[i - 1].second == Ciclo[j].second)
            ) {
                ladoCambiar = j
                break
            }
        }
        val tmp = Ciclo[i]
        Ciclo[i] = Ciclo[ladoCambiar]
        Ciclo[ladoCambiar] = tmp
    }
}

/**
 * Algoritmo 3: Combina dos ciclos, luego de revisar cada combinación de pares de lados 
 * de los ciclos y sustituirlos por dos nuevos lados que unan a los ciclos.
 * Precondición: Ciclo1.size >= 0 && Ciclo2.size >= 0
 * Postcondición: Ciclo3.size == Ciclo1.size + Ciclo2.size
 *
 * @param[Ciclo1] Arreglo de pares de pares de doubles, arreglo de lados.
 * @param[Ciclo2] Arreglo de pares de pares de doubles, arreglo de lados.
 *
 * @return Arreglo de pares de pares de doubles, que representan un ciclo, arreglo de lados. 
 */
fun combinarCiclos(
    Ciclo1: Array<Pair<Pair<Double, Double>, Pair<Double, Double>>>, 
    Ciclo2: Array<Pair<Pair<Double, Double>, Pair<Double, Double>>>
): Array<Pair<Pair<Double, Double>, Pair<Double, Double>>> {
    if (Ciclo1.size == 0) return Ciclo2
    if (Ciclo2.size == 0) return Ciclo1
    var minG = POSITIVE_INFINITY

    // Inicializar los lados a agregar y eliminar
    var ladosAgregarC1 = Ciclo1[0]
    var ladosAgregarC2 = Ciclo2[0]
    var ladosEliminarC1 = Ciclo1[0]
    var ladosEliminarC2 = Ciclo2[0]

    for (lado1 in Ciclo1) {
        val (a, b) = lado1
        val dOLD1 = distancia(a, b)         // Distancia actual del lado 1

        for (lado2 in Ciclo2) {
            val (c, d) = lado2
            val dOLD2 = distancia(c, d)     // Distancia actual del lado 2

            // Calcular las distancias de los posibles lados
            val dNEW1 = distancia(a, c) 
            val dNEW2 = distancia(b, d) 
            val dNEW3 = distancia(a, d) 
            val dNEW4 = distancia(b, c)

            // Calcular la distancia ganada
            val g1 = distanciaGanada(dOLD1, dOLD2, dNEW1, dNEW2)
            val g2 = distanciaGanada(dOLD1, dOLD2, dNEW3, dNEW4)

            // Si la ganancia es menor a la anterior, actualizar los lados a eliminar y agregar
            val ganancia = min(g1, g2)
            if (ganancia < minG) {
                minG = ganancia
                if (g1 < g2) {
                    ladosAgregarC1 = Pair(a, c)
                    ladosAgregarC2 = Pair(b, d)
                } else {
                    ladosAgregarC1 = Pair(a, d)
                    ladosAgregarC2 = Pair(b, c)
                }
                ladosEliminarC1 = Pair(a, b)
                ladosEliminarC2 = Pair(c, d)
            }
        }
    }
    // Agregar a Ciclo3 los lados de Ciclo1 y Ciclo2, sin tomar en cuenta ladosEliminarC1
    // y ladosEliminarC2, pero si agregando ladosAgregarC1 y ladosAgregarC2
    val (n1, n3) = Pair(Ciclo1.size, Ciclo1.size + Ciclo2.size)
    var (eliminado1, eliminado2) = Pair(false, false)
    val Ciclo3 = Array<Pair<Pair<Double, Double>, Pair<Double, Double>>>(n3, {
            if (it < n1) {
                if (!(eliminado1) && Ciclo1[it] == ladosEliminarC1) {
                    eliminado1 = true
                    ladosAgregarC1
                } else {
                    Ciclo1[it]
                }
            } else {
                if (!(eliminado2) && Ciclo2[it - n1] == ladosEliminarC2) {
                    eliminado2 = true
                    ladosAgregarC2
                } else {
                    Ciclo2[it - n1]
                }
            }
        }
    )
    // Ordenar Ciclo3 con una modificación de Selection Sort.
    selectionSortCiclo(Ciclo3)
    return Ciclo3
}


/**
 * Algoritmo 2: Divide un arreglo de ciudades en dos particiones. Se divide el espacio
 * en dos rectángulos, se obtienen las ciudades incluídas en cada partición.
 * Precondición: P.size >= 4
 * Postcondición: particionIzq.size > 0 && particionDer.size > 0
 *
 * @param[P] Arreglo de pares de doubles, arreglo de ciudades.
 *
 * @return Par de dos arreglos de pares de doubles, es decir, dos arreglos de ciudades. 
 */
fun obtenerParticiones(
    P: Array<Pair<Double, Double>>
): Pair<Array<Pair<Double, Double>>, Array<Pair<Double, Double>>> {
    val rectangulo = rectanguloCoords(P)                // Coords del rectángulo
    val (xDim, yDim) = rectanguloDim(rectangulo)        // Dimensiones del rectángulo
    var ejeDeCorte = if (xDim > yDim) "X" else "Y"      // Eje de corte según el lado más largo
    var (xc, yc) = obtenerPuntoDeCorte(P, ejeDeCorte)   // Buscar el punto de corte

    // Obtener el rectángulo izquierdo/derecho o superior/inferior
    var (rectanguloIzq, rectanguloDer) = aplicarCorte(ejeDeCorte, Pair(xc, yc), rectangulo)

    // Obtener las ciudades de los dos rectángulos
    var particionIzq = obtenerPuntosRectangulo(P, rectanguloIzq) 
    var particionDer = obtenerPuntosRectangulo(P, rectanguloDer)

    // Si una partición no tiene ciudades, realizar el corte de nuevo con el ejeDeCorte contrario
    if ((particionIzq.size == 0 && particionDer.size > 3) || 
        (particionIzq.size > 3 && particionDer.size == 0)
    ) {
        ejeDeCorte = if (ejeDeCorte == "X") "Y" else "X"
        var puntoDeCorteNuevo = obtenerPuntoDeCorte(P, ejeDeCorte)
        xc = puntoDeCorteNuevo.first
        yc = puntoDeCorteNuevo.second

        var rectangulosNuevo = aplicarCorte(ejeDeCorte, Pair(xc, yc), rectangulo)
        rectanguloIzq = rectangulosNuevo.first
        rectanguloDer = rectangulosNuevo.second

        particionIzq = obtenerPuntosRectangulo(P, rectanguloIzq)
        particionDer = obtenerPuntosRectangulo(P, rectanguloDer)

        // Si una partición no tiene ciudades, realizar el corte por la mitad
        if ((particionIzq.size == 0 && particionDer.size > 3) || 
            (particionIzq.size > 3 && particionDer.size == 0)
        ) {
            puntoDeCorteNuevo = obtenerPuntoDeCorteMitad(rectangulo, ejeDeCorte)
            xc = puntoDeCorteNuevo.first
            yc = puntoDeCorteNuevo.second

            rectangulosNuevo = aplicarCorte(ejeDeCorte, Pair(xc, yc), rectangulo)
            rectanguloIzq = rectangulosNuevo.first
            rectanguloDer = rectangulosNuevo.second
            
            particionIzq = obtenerPuntosRectangulo(P, rectanguloIzq)
            particionDer = obtenerPuntosRectangulo(P, rectanguloDer)
        } 
    }
    return Pair(particionIzq, particionDer)
}

/**
 * Algoritmo 1: Resuelve el TSP euleriano de manera recursiva.
 * Construye un ciclo de lados (par de ciudades).
 * Entrada: arreglo con ciudades. Salida: ciclo solución válida TSP
 * Precondición: P.size >= 0
 * Postcondición: P.size == tamaño del ciclo obtenido.
 * 
 * @param[P] Arreglo de pares de doubles, arreglo de ciudades.
 *
 * @return Arreglos de lados (pares de pares de doubles, par de dos ciudades), solución TSP.
 */
fun divideAndConquerTSP(
    P: Array<Pair<Double, Double>>
): Array<Pair<Pair<Double, Double>, Pair<Double, Double>>> {
    val n = P.size
    when (n) {
        0 -> return arrayOf()
        1 -> return cicloUnaCiudad(P)
        2 -> return cicloDosCiudades(P)
        3 -> return cicloTresCiudades(P)
        else -> {
            val (pizquierda, pderecha) = obtenerParticiones(P)
            val c1 = divideAndConquerTSP(pizquierda)
            val c2 = divideAndConquerTSP(pderecha)
            return combinarCiclos(c1, c2)
        }
    }
}

/**
 * Verifica que una solución TSP sea correcta.
 * Precondición: P.size >= 0 && sol.size >= 0
 * Postcondición: P.size == tamaño del ciclo obtenido.
 * 
 * @param[primeraCiudad] Par de doubles, representa una ciuda.
 * @param[P] Arreglo de pares de doubles, arreglo de ciudades.
 * @param[tour] Arreglo de pares de pares de doubles, arreglo de lados solución del TSP.
 *
 * @return Boolean. True si la solución TSP es válida, False si la solución no es válida.
 */
fun verificadorTSP(
    primeraCiudad: Pair<Double, Double>,
    P: Array<Pair<Double, Double>>, 
    tour: Array<Pair<Pair<Double, Double>, Pair<Double, Double>>>
): Boolean {
    val n = tour.size
    // Verificar que el tamaño de la solución es igual al tamaño de ciudades
    if (P.size != n) {
        println("Error: El tamaño del tour TSP debe ser igual al número de ciudades")
        return false
    }

    // Verificar que inicie y finalice en la primera ciudad
    if ((primeraCiudad != tour[0].first && primeraCiudad != tour[0].second) ||
        (primeraCiudad != tour[n - 1].first && primeraCiudad != tour[n - 1].second)
    ) {
        println("Error: La solución del tour TSP debe iniciar y terminar con la primera ciudad.")
        return false
    }

    // Verificar que cada ciudad sea visitado una vez y que estén todas las ciudades.
    for (ciudad in P) {
        var counter = 0
        for (lado in tour) {
            if (ciudad == lado.first || ciudad == lado.second) {
                counter++
            }
        }
        if (counter != 2) {
            println("Error: La ciudad ${ciudad} no está siendo visitada o está siendo visitada más de una vez.")
            return false
        }
    }

    // Verificar que cada lado conecte con el siguiente lado
    var ciudades = arrayOf(tour[0].first, tour[0].second)
    for (i in 1 until n) {
        if (!(tour[i].first in ciudades || tour[i].second in ciudades)) {
            println("Error: el lado ${tour[i]} no conecta con la anterior ciudad")
            return false
        } else {
            ciudades = arrayOf(tour[i].first, tour[i].second)
        }
    }
    return true
}

/**
 * Extrae las coordenadas de las ciudades dado un archivo path.
 * Precondición: filename sea un archivo existente.
 * Postcondición: Tamaño del arreglo de ciudades es igual a DIMENSION del archivo
 *
 * @param[filename] Dirección del archivo
 *
 * @return Arreglo de las coordenadas (pares de doubles) de ciudades
 */
fun extraerCoords(filename: String): Array<Pair<Double, Double>> {
    val datos: Array<String> = File(filename).readLines().toTypedArray()
    val dimension = datos[3].replace(Regex("[^0-9]"), "").toInt()
    val ciudades = Array<Pair<Double, Double>>(dimension, {
        val coord = (datos[it + 6].split(" ").filter {it != ""}).toTypedArray()
        Pair(coord[1].toDouble(), coord[2].toDouble())
    })
    return ciudades
}

/**
 * Retorna un arreglo con los índices de las ciudades en orden como fueron visitadas en el tour
 * Precondición: ciudades.size == tour.size
 * Postcondición: El arreglo a retornar debe estar constituído por números enteros del 
 * 1 al ciudades.size (sin repetir), que representan cada ciudad. 
 * El arreglo debe estar en orden del como fueron visitadas en el tour.
 *
 * @param[ciudades] Arreglo de pares de doubles, que son las coordenadas de las ciudades
 * @param[tour] Arreglo de pares de pares de double que es la solución tour TSP
 *
 * @return Arreglo de enteros, que son los índices de las ciudades según el tour TSP.
 */
fun tourSolucion(
    ciudades: Array<Pair<Double, Double>>, 
    tour: Array<Pair<Pair<Double, Double>, Pair<Double, Double>>>
): IntArray {
    val n = ciudades.size
    val tourSol = IntArray(n)
    tourSol[0] = 1
    var tourIndex = 1
    var ciudadAnterior = ciudades[0]

    // Recorremos el tour, determinando de que ciudad a que ciudad va
    for (t in tour) {
        if (tourIndex < n) {
            val ciudadActual = if (t.first != ciudadAnterior) t.first else t.second
            // Buscamos la posición de la ciudad
            var i = 0
            while (i < n && ciudades[i] != ciudadActual) {
                i++
            }
            tourSol[tourIndex] = i + 1
            tourIndex++
            ciudadAnterior = ciudadActual
        }
    }
    return tourSol
}

/**
 * Dado un archivo, la distancia total de un tour y un arreglo con el orden en que fueron 
 * visitadas las ciudades en el tour solución TSP, escribe dentro del archivo estos datos.
 * Precondición: instancias.size >= 0  && distanciaTotal > 0
 * Postcondición: Archivo out creado con la solución del tour.
 *
 * @param[filename] String del nombre del archivo a crear.
 * @param[distanciaTotal] Número entero de la suma de todas las distancias.
 * @param[instancias] Arreglo de los números índices de las ciudades según la solución TSP.
 */
fun escribirOutfile(filename: String, distanciaTotal: Int, instancias: IntArray) {
    // Crear un archivo con el filename dado
    var file = File(filename)
    file.createNewFile()
    
    // Agregar los datos de la solución TSP al archivo creado
    val dataPart = arrayOf(
        "NAME : ${filename}", "COMMENT : Length ${distanciaTotal}", "TYPE : TOUR",
        "DIMENSION : ${instancias.size}", "TOUR_SECTION"
    )
    file.printWriter().use { out ->
        dataPart.forEach {
            out.println("${it}")
        }
        instancias.forEach {
            out.println("${it}")
        }
        out.println("-1")
        out.println("EOF")
    }
}

/**
 * Main: Dado un archivo de instancias TSP y un archivo de salida, se obtiene la solución TSP de
 * dicha instancia y se escribe la solución en el archivo de salida.
 * Precondición: args.size == 2
 * Postcondición: Se imprime por consola el nombre de la instancia TSP y la distancia total.
 *
 * @param[args] Arreglo de String, que contiene la dirección del archivo de entrada y de salida.
 */
fun main(args: Array<String>) {
    val ciudades = extraerCoords(args[0])       // Coordenadas de las ciudades
    val ciudadesCopia = ciudades.copyOf()       // Copia de las ciudades originales
    var tour = divideAndConquerTSP(ciudades)    // Solución TSP (se modificará el orden de ciudades)

    // Asegurar que el tour debe comenzar con la primera ciudad del archivo
    val n = tour.size
    var ciudad1 = 0
    for (i in 0 until n) {
        if (i == 0 && (tour[i].first == ciudadesCopia[0] || tour[i].second == ciudadesCopia[0]) &&
            (i + 1 < n && tour[i + 1].first != ciudadesCopia[0] && tour[i + 1].second != ciudadesCopia[0])
        ) {
            break
        } else if (tour[i].first == ciudadesCopia[0] || tour[i].second == ciudadesCopia[0]) {
            ciudad1 = i + 1
            break
        }
    }
    tour = tour.sliceArray(ciudad1..(n - 1)) + tour.sliceArray(0..(ciudad1 - 1))

    //val verificador = verificadorTSP(ciudadesCopia[0], ciudades, tour)
    //if (!verificador) return
    
    // Calcular las distancias
    var distanciaTotal = 0
    for (t in tour) distanciaTotal += (distancia(t.first, t.second) + 0.5).toInt()
    println("TSP instance: ${args[0]}")
    println("Total distance: ${distanciaTotal}")

    // Buscamos las soluciones del tour (índices de la ciudad) y escribe en el archivo de salida.
    val solCiudades = tourSolucion(ciudadesCopia, tour)
    escribirOutfile(args[1], distanciaTotal, solCiudades) 
}