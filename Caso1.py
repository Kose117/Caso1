import numpy as np
from collections import deque
from pulp import *
import pandas as pd
import networkx as nx
from pulp import LpProblem, LpMinimize, LpVariable, lpSum, LpStatus
from collections import deque


class FordFulkersonWithDelay:
    def __init__(self, grafo, retardos):
        self.grafo = grafo
        # Creando una copia del grafo para el grafo residual
        self.grafo_residual = [fila[:] for fila in grafo]
        self.num_nodos = len(grafo)
        self.retardos = retardos
        self.total_retardo = 0  # Inicializar el retardo total

    def bfs(self, fuente, sumidero, padre):
        visitado = [False] * self.num_nodos
        cola = [fuente]
        visitado[fuente] = True

        while cola:
            u = cola.pop(0)

            for v, capacidad in enumerate(self.grafo_residual[u]):
                # Si no esta visitado y hay capacidad
                if not visitado[v] and capacidad > 0:
                    cola.append(v)
                    visitado[v] = True
                    padre[v] = u
                    if v == sumidero:
                        return True
        return False

    def ford_fulkerson(self, fuente, sumidero):
        padre = [-1] * self.num_nodos
        flujo_maximo = 0

        # Aumentar el flujo mientras haya un camino de la fuente al sumidero
        while self.bfs(fuente, sumidero, padre):
            flujo_camino = float('Inf')
            s = sumidero
            camino = []

            # Encontrar el flujo maximo a traves del camino encontrado
            while s != fuente:
                flujo_camino = min(
                    flujo_camino, self.grafo_residual[padre[s]][s])
                camino.insert(0, s)
                s = padre[s]

            camino.insert(0, fuente)  # Anadir la fuente al inicio del camino

            # Calcular el retardo total para este camino
            retardo_camino = sum(
                self.retardos[camino[i]][camino[i + 1]] for i in range(len(camino) - 1))
            self.total_retardo += retardo_camino

            # Actualizar las capacidades residuales de las aristas y las aristas inversas a lo largo del camino
            v = sumidero
            while v != fuente:
                u = padre[v]
                self.grafo_residual[u][v] -= flujo_camino
                self.grafo_residual[v][u] += flujo_camino
                v = padre[v]

            flujo_maximo += flujo_camino

            # Imprimir el camino tomado y el flujo a traves de ese camino
            print(f"Camino tomado para el nodo {
                  sumidero+1}: {camino}, Flujo: {flujo_camino}, Retardo: {retardo_camino}")

        return flujo_maximo, self.total_retardo

# --------------------- solucion punto 1 ---------------------


def punto1_algoritmo():
    print("\nPUNTO 1 - CYCLE CANCELLING")
    cost_df = pd.read_csv('retardos.csv')
    capacity_df = pd.read_csv('grafo.csv')

    G = nx.DiGraph()

    num_nodes = cost_df.shape[0]

    for i in range(num_nodes):
        for j in range(num_nodes):
            cost = cost_df.iloc[i, j]
            capacity = capacity_df.iloc[i, j]
            if cost != 999 and capacity > 0:  # Agregar arista solo si hay un camino valido con capacidad no cero
                G.add_edge(i+1, j+1, weight=cost, capacity=capacity)

    source = 1
    targets = [49, 50]
    demand = 90  # kb/s de trafico a enviar

    flow_dict = {edge: 0 for edge in G.edges}

    def find_min_cost_flow_corrected(G, source, targets, demand):
        residual_G = G.copy()

        for target in targets:
            path = nx.shortest_path(
                residual_G, source=source, target=target, weight='weight')
            path_edges = list(zip(path, path[1:]))
            min_capacity = min(residual_G[u][v]['capacity']
                               for u, v in path_edges)

            flow = min(min_capacity, demand / len(targets))
            for u, v in path_edges:
                flow_dict[(u, v)] += flow
                residual_G[u][v]['capacity'] -= flow
                if residual_G[u][v]['capacity'] == 0:
                    residual_G.remove_edge(u, v)

        while True:
            try:
                cycle_edges = nx.find_cycle(residual_G, orientation='ignore')
                cycle_cost = sum(residual_G[u][v]['weight']
                                 for u, v, _ in cycle_edges)

                if cycle_cost >= 0:
                    break

                min_cycle_capacity = min(
                    residual_G[u][v]['capacity'] for u, v, _ in cycle_edges)

                for u, v, direction in cycle_edges:
                    if direction == 'forward':
                        if (u, v) in flow_dict:
                            flow_dict[(u, v)] += min_cycle_capacity
                        else:
                            flow_dict[(v, u)] -= min_cycle_capacity
                    else:
                        if (v, u) in flow_dict:
                            flow_dict[(v, u)] += min_cycle_capacity
                        else:
                            flow_dict[(u, v)] -= min_cycle_capacity

                    residual_G[u][v]['capacity'] -= min_cycle_capacity
                    if residual_G[u][v]['capacity'] == 0:
                        residual_G.remove_edge(u, v)

            except nx.NetworkXNoCycle:
                break

        return flow_dict

    min_cost_flow_corrected = find_min_cost_flow_corrected(
        G, source, targets, demand)

    min_cost_paths_corrected = {}
    total_delay = {}

    for target in targets:
        path = nx.shortest_path(
            G, source=source, target=target, weight='weight')
        path_flow = min(G[u][v]['capacity'] for u, v in zip(path, path[1:]))
        delay = sum(G[u][v]['weight'] for u, v in zip(path, path[1:]))
        min_cost_paths_corrected[target] = (path, path_flow, delay)

    print("Minimum Cost Flow Paths with Total Delay:")
    for target, (path, flow, delay) in min_cost_paths_corrected.items():
        print(f"Target Node {target}:")
        print(f"  Path: {path}")
        print(f"  Flow: {flow} kb/s")
        print(f"  Total Delay: {delay}\n")


def solucion_punto_1():
    print("\nPUNTO 1 - MODELO MATEMÁTICO")
    # Lectura de las matrices de capacidades y retardos desde archivos Excel
    capacidades = pd.read_excel("Exceles/C2.xlsx").to_numpy()
    retardos = pd.read_excel("Exceles/R2.xlsx").to_numpy()
    num_nodos = len(capacidades[0])

    # Definición del problema de optimización para minimizar el costo
    problema = LpProblem('MinCost', LpMinimize)

    # Definir nodos y arcos en la red
    nodos = range(num_nodos)
    arcos = [(i, j) for i in nodos for j in nodos]

    # Variables de decisión para los flujos en cada arco con límites y tipo
    flujos = LpVariable.dicts("Flujo", (nodos, nodos),
                              lowBound=0, cat="Continuous")

    # Función objetivo: minimizar el costo total del flujo en la red
    problema += lpSum(flujos[i][j] * retardos[i][j] for (i, j) in arcos)

    # Restricciones del problema: flujo total en nodos de entrada y salida
    problema += lpSum(flujos[0][j] for j in nodos if (0, j) in arcos) == 90
    problema += lpSum(flujos[i][destino] for destino in [48, 49]
                      for i in nodos if (i, destino) in arcos) == 90

    # Restricciones de capacidad en los arcos
    for i, j in arcos:
        problema += flujos[i][j] <= capacidades[i][j]

    # Restricciones de conservación del flujo para nodos intermedios
    for k in nodos:
        if k != 0 and k not in [48, 49]:
            problema += (
                lpSum(flujos[i][k] for i in nodos if (i, k) in arcos) ==
                lpSum(flujos[k][j] for j in nodos if (k, j) in arcos)
            )

    # Resolver el problema de optimización
    problema.solve()

    # Imprimir el estado de la solución
    print(f"\nEstado de la solución: {LpStatus[problema.status]}\n")

    # Construcción de la red de flujos para análisis posterior
    red_flujo = {i: [] for i in nodos}
    for var in problema.variables():
        if var.varValue > 0:
            i, j = map(int, var.name.split('_')[1:])
            red_flujo[i].append(j)

    # Función para encontrar caminos mediante BFS
    def encontrar_camino(red, inicio, destinos):
        cola = deque([(inicio, [inicio])])
        visitados = set()

        while cola:
            nodo_actual, camino = cola.popleft()
            for siguiente in red[nodo_actual]:
                if siguiente not in visitados:
                    if siguiente in destinos:
                        return camino + [siguiente]
                    visitados.add(siguiente)
                    cola.append((siguiente, camino + [siguiente]))
        return None

    # Encontrar y mostrar caminos hacia los nodos destino 48 y 49
    print("Caminos desde el nodo 1 a los nodos destino:")
    print("-" * 40)
    for destino in [48, 49]:
        camino = encontrar_camino(red_flujo, 0, {destino})
        if camino:
            print(f"Camino hasta el nodo {
                  destino + 1}: {[n + 1 for n in camino]}")
        else:
            print(f"No se encontró un camino hasta el nodo {destino + 1}")
    print("-" * 40)

    # Mostrar el valor mínimo de la función objetivo
    print(f"\nEl costo mínimo total del flujo es: {value(problema.objective)}")

    
    punto1_algoritmo()


def punto2_algoritmo():
    print("\nPUNTO 2 - FORD FULKERSON")
    ruta_archivo_grafo = 'grafo.csv'
    ruta_archivo_retardos = 'retardos.csv'
    datos_grafo = pd.read_csv(ruta_archivo_grafo)
    datos_retardos = pd.read_csv(ruta_archivo_retardos)

    grafo = datos_grafo.values.tolist()
    retardos = datos_retardos.values.tolist()

    ff_delay = FordFulkersonWithDelay(grafo, retardos)

    flujo_max_1_a_49, retardo_total_1_a_49 = ff_delay.ford_fulkerson(0, 48)
    flujo_max_1_a_50, retardo_total_1_a_50 = ff_delay.ford_fulkerson(0, 49)

    flujo_max_total_ff = flujo_max_1_a_49 + flujo_max_1_a_50
    retardo_total_ff = retardo_total_1_a_49 + retardo_total_1_a_50

    print(f"Flujo maximo total: {flujo_max_total_ff}")
    print(f"Retardo total en la red: {retardo_total_ff}")

# --------------------- solucion punto 2 ---------------------


def solucion_punto_2():
    print("\nPUNTO 2 - MODELO MATEMÁTICO")
    # Leer los datos de capacidad y costo desde archivos Excel
    capacidad_red = pd.read_excel("./Exceles/C2.xlsx").to_numpy().tolist()
    costo_red = pd.read_excel("./Exceles/R2.xlsx").to_numpy().tolist()
    num_nodos = len(capacidad_red[0])

    # Definir el problema de optimización para maximizar el flujo
    problema_flujo_max = LpProblem("Maximo_Flujo", LpMaximize)

    # Crear nodos y arcos posibles
    nodos = range(num_nodos)
    arcos = [(origen, destino) for origen in nodos for destino in nodos]

    # Definir variables de flujo con límites y tipo
    flujo_variables = LpVariable.dicts(
        "Flujo", (nodos, nodos), lowBound=0, cat="Integer")

    # Función objetivo: maximizar el flujo hacia los nodos 48 y 49
    problema_flujo_max += (
        lpSum(flujo_variables[origen][48] for origen in nodos) +
        lpSum(flujo_variables[destino][49] for destino in nodos)
    )

    # Restricciones de capacidad por arco
    for origen in nodos:
        for destino in nodos:
            problema_flujo_max += flujo_variables[origen][destino] <= capacidad_red[origen][destino]

    # Restricciones de conservación de flujo para nodos intermedios (excepto nodos 0, 48, y 49)
    for nodo in nodos:
        if nodo > 0 and nodo not in (48, 49):
            problema_flujo_max += (
                lpSum(flujo_variables[origen][nodo] for origen in nodos) ==
                lpSum(flujo_variables[nodo][destino] for destino in nodos)
            )

    # Resolver el problema de optimización
    estado = problema_flujo_max.solve()

    # Imprimir el estado de la solución
    print(f"\nEstado de la solución: {LpStatus[problema_flujo_max.status]}\n")

    # Imprimir variables con valores positivos (flujos significativos) de manera organizada
    print("Flujos significativos:")
    print("-" * 40)
    print(f"{'Origen':<10} {'Destino':<10} {'Flujo':<10}")
    print("-" * 40)
    for variable in problema_flujo_max.variables():
        if variable.value() > 0:
            # Extraer origen y destino del nombre de la variable
            origen, destino = variable.name.split('_')[1:]
            print(f"{origen:<10} {destino:<10} {value(variable):<10}")
    print("-" * 40)

    # Mostrar el valor del flujo máximo encontrado
    print(f"\nEl flujo maximo posible enviado es: {
          value(problema_flujo_max.objective)}")

    # Calcular e imprimir el retardo total en la red
    retardo_total = sum(
        value(flujo_variables[origen][destino]) * costo_red[origen][destino]
        for origen in nodos for destino in nodos if value(flujo_variables[origen][destino]) > 0
    )
    print(f"\nEl retardo total en la red al enviar el flujo maximo es: {
          retardo_total}")

    punto2_algoritmo()

# --------------------- solucion punto 3, 4, 5, 6, 7 y 8 ---------------------


def modelo_matematico_dijkstra(csv, destination, num):
    # Lectura y conversión de la matriz de retardos desde el archivo CSV con el separador correcto
    retardos = pd.read_csv(csv+".csv", header=None, sep=',', skiprows=1).apply(
        pd.to_numeric, errors='coerce').replace(999, np.inf).to_numpy()
    num_nodos = len(retardos[0])

    # Definición del problema de optimización
    problema = LpProblem('RutaMasCorta', LpMinimize)

    # Nodos en la red
    nodos = range(num_nodos)

    # Condicional para arcos: si num == 7, filtrar valores mayores a 40
    if num == 7:
        arcos = [(i, j) for i in nodos for j in nodos if retardos[i]
                 [j] > 40 and retardos[i][j] < np.inf]
    else:
        arcos = [(i, j)
                 for i in nodos for j in nodos if retardos[i][j] < np.inf]

    # Variables de decisión: 1 si se usa el arco, 0 en caso contrario
    x = LpVariable.dicts("x", (nodos, nodos), cat="Binary")

    # Función objetivo: minimizar el retardo total
    problema += lpSum(x[i][j] * retardos[i][j] for (i, j) in arcos)

    # Restricción de flujo en el nodo fuente (nodo 0)
    problema += lpSum(x[0][j] for j in nodos if (0, j) in arcos) == 1

    # Restricción de flujo en el nodo destino
    problema += lpSum(x[i][destination]
                      for i in nodos if (i, destination) in arcos) == 1

    # Conservación del flujo en nodos intermedios
    for k in nodos:
        # Excluyendo nodo fuente (0) y nodo destino (48)
        if k != 0 and k != destination:
            problema += lpSum(x[i][k] for i in nodos if (i, k)
                              in arcos) == lpSum(x[k][j] for j in nodos if (k, j) in arcos)

    # Resolución del problema
    problema.solve()
    print("Estado solución:", LpStatus[problema.status])

    # Mostrar la ruta más corta encontrada
    ruta = []
    for i in nodos:
        for j in nodos:
            if x[i][j].varValue == 1:
                ruta.append((i, j))

    # Convertir nodos a formato 1-indexado para la salida
    ruta_mas_corta = [(i + 1, j + 1) for (i, j) in ruta]
    print("Ruta más corta desde el nodo 1 hasta el nodo " +
          str(destination) + ":", ruta_mas_corta)

    # Mostrar el valor de la función objetivo
    print("El retardo mínimo total es:", value(problema.objective))


def solucion_dijkstra(csv, source, destination, num):
    print(f"\nPUNTO {num} - DIJKSTRA")
    file_path = csv + '.csv'
    data = pd.read_csv(file_path)

    G = nx.Graph()

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            # Condicional para el punto 7, solo considera valores mayores a 40
            if num == 7 and data.iloc[i, j] <= 40:
                continue
            if data.iloc[i, j] != 999:  # 999 indica que no hay conexión
                G.add_edge(i, j, weight=data.iloc[i, j])

    try:
        path = nx.dijkstra_path(G, source, destination, weight='weight')
        delay = nx.dijkstra_path_length(
            G, source, destination, weight='weight')

        print(
            f"La ruta óptima desde el nodo {source + 1} al nodo {destination + 1} es:", path)
        print("El costo total de la ruta óptima es:", delay)

        # Llamar a la función de costos del servidor solo si num == 8
        if num == 8:
            calcular_costos_servidor(path, delay)

    except nx.NetworkXNoPath:
        print("No existe una ruta óptima entre los nodos especificados.")


def calcular_costos_servidor(path, delay):
    # Costos del servidor
    costoServidor = [57, 40, 34, 82, 61, 34, 65, 61, 52, 84, 68, 65, 49, 76, 50, 82, 84, 81,
                     61, 52, 60, 42, 39, 49, 57, 58, 82, 39, 61, 54, 67, 60, 65, 26, 83, 76,
                     83, 67, 81, 50, 66, 29, 36, 28, 31, 44, 29, 78, 39, 64]

    # Calcular la suma de los costos del servidor en la ruta óptima
    sum = 0
    for i in range(len(path)-1):
        sum += costoServidor[path[i]]

    print("La suma de los costos del servidor en la ruta optima es:", sum)
    print("Retardo total de la red es:", delay-sum)


# --------------------- Menu ---------------------
def main_menu():
    while True:
        print("\nMenu de Preguntas:")
        print("1. Cual es el flujo de red que minimiza el retardo entre el nodo 1 (fuente) y los nodos 49 y 50?")
        print("2. Cual es el trafico maximo que soporta la red entre el nodo 1 (fuente) y los nodos 49 y 50? Cual seria el retardo total?")
        print("3. Cual seria la ruta optima a minimo retardo entre el nodo 1 (fuente) y el nodo 49?")
        print("4. Cual seria la ruta optima a minimo retardo entre el nodo 1 (fuente) y el nodo 50?")
        print("5. Cual seria la ruta a minimo costo entre el nodo 1 (fuente) y el nodo 49?")
        print("6. Cual seria la ruta optima a minimo costo entre el nodo 1 (fuente) y el nodo 50?")
        print("7. Como cambian las soluciones de los puntos 3 y 4 si el retardo debe ser mayor a 40?")
        print("8. Cual es la ruta que tiene minimo costo y minimo retardo?")
        print("9. Salir")

        opcion = input("Seleccione una opcion (1-9): ")

        if opcion == '1':
            solucion_punto_1()

        elif opcion == '2':
            solucion_punto_2()

        elif opcion == '3':
            solucion_dijkstra("retardos", 0, 48, 3)
            print("\n/--------------------------------/")
            modelo_matematico_dijkstra("retardos", 48, 3)

        elif opcion == '4':
            solucion_dijkstra("retardos", 0, 49, 4)
            print("\n/--------------------------------/")
            modelo_matematico_dijkstra("retardos", 49, 4)

        elif opcion == '5':
            solucion_dijkstra("matriz_transformada", 0, 48, 5)
            print("\n/--------------------------------/")
            modelo_matematico_dijkstra("matriz_transformada", 48, 5)

        elif opcion == '6':
            solucion_dijkstra("matriz_transformada", 0, 49, 6)
            print("\n/--------------------------------/")
            modelo_matematico_dijkstra("matriz_transformada", 49, 6)

        elif opcion == '7':
            solucion_dijkstra("retardos", 0, 48, 7)
            modelo_matematico_dijkstra("retardos", 48, 7)
            print("\n/--------------------------------/")
            solucion_dijkstra("retardos", 0, 49, 7)
            modelo_matematico_dijkstra("retardos", 49, 7)

        elif opcion == '8':
            solucion_dijkstra("matriz_suma", 0, 48, 8)
            modelo_matematico_dijkstra("matriz_suma", 48, 8)
            print("\n/--------------------------------/")
            solucion_dijkstra("matriz_suma", 0, 49, 8)
            modelo_matematico_dijkstra("matriz_suma", 49, 8)

        elif opcion == '9':
            print("Saliendo del programa.")
            break

        else:
            print("Opcion no valida. Por favor, seleccione una opcion entre 1 y 9.")


if __name__ == "__main__":
    main_menu()
