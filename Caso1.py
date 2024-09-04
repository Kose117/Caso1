import pandas as pd
import networkx as nx

class FordFulkersonWithDelay:
    def __init__(self, grafo, retardos):
        self.grafo = grafo
        self.grafo_residual = [fila[:] for fila in grafo]  # Creando una copia del grafo para el grafo residual
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
                if not visitado[v] and capacidad > 0:  # Si no esta visitado y hay capacidad
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
                flujo_camino = min(flujo_camino, self.grafo_residual[padre[s]][s])
                camino.insert(0, s)
                s = padre[s]

            camino.insert(0, fuente)  # Anadir la fuente al inicio del camino

            # Calcular el retardo total para este camino
            retardo_camino = sum(self.retardos[camino[i]][camino[i + 1]] for i in range(len(camino) - 1))
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
            print(f"Camino tomado para el nodo {sumidero+1}: {camino}, Flujo: {flujo_camino}, Retardo: {retardo_camino}")

        return flujo_maximo, self.total_retardo

# --------------------- solucion punto 1 ---------------------
def solucion_punto_1():
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
            path = nx.shortest_path(residual_G, source=source, target=target, weight='weight')
            path_edges = list(zip(path, path[1:]))
            min_capacity = min(residual_G[u][v]['capacity'] for u, v in path_edges)
            
            flow = min(min_capacity, demand / len(targets))
            for u, v in path_edges:
                flow_dict[(u, v)] += flow
                residual_G[u][v]['capacity'] -= flow
                if residual_G[u][v]['capacity'] == 0:
                    residual_G.remove_edge(u, v)
                    
        while True:
            try:
                cycle_edges = nx.find_cycle(residual_G, orientation='ignore')
                cycle_cost = sum(residual_G[u][v]['weight'] for u, v, _ in cycle_edges)
                
                if cycle_cost >= 0:
                    break
                
                min_cycle_capacity = min(residual_G[u][v]['capacity'] for u, v, _ in cycle_edges)
                
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

    min_cost_flow_corrected = find_min_cost_flow_corrected(G, source, targets, demand)

    min_cost_paths_corrected = {}
    total_delay = {}

    for target in targets:
        path = nx.shortest_path(G, source=source, target=target, weight='weight')
        path_flow = min(G[u][v]['capacity'] for u, v in zip(path, path[1:]))
        delay = sum(G[u][v]['weight'] for u, v in zip(path, path[1:]))
        min_cost_paths_corrected[target] = (path, path_flow, delay)

    print("Minimum Cost Flow Paths with Total Delay:")
    for target, (path, flow, delay) in min_cost_paths_corrected.items():
        print(f"Target Node {target}:")
        print(f"  Path: {path}")
        print(f"  Flow: {flow} kb/s")
        print(f"  Total Delay: {delay}\n")

# --------------------- solucion punto 2 ---------------------
def solucion_punto_2():
    print("\nPUNTO 2 - FLUJO MAXIMO Y RETARDO TOTAL")
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

# --------------------- solucion punto 3 ---------------------
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
        delay = nx.dijkstra_path_length(G, source, destination, weight='weight')

        print(f"La ruta óptima desde el nodo {source + 1} al nodo {destination + 1} es:", path)
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
            solucion_dijkstra("retardos",0,48, 3)

        elif opcion == '4':
            solucion_dijkstra("retardos",0,49, 4)

        elif opcion == '5':
            solucion_dijkstra("matriz_transformada",0,48, 5)
            

        elif opcion == '6':
            solucion_dijkstra("matriz_transformada",0,49, 6)

        elif opcion == '7':
            solucion_dijkstra("retardos",0,48, 7)
            print("\n/--------------------------------/")
            solucion_dijkstra("retardos",0,49, 7)

        elif opcion == '8':
            solucion_dijkstra("matriz_suma",0,48, 8)
            print("\n/--------------------------------/")
            solucion_dijkstra("matriz_suma",0,49, 8)

        elif opcion == '9':
            print("Saliendo del programa.")
            break
        
        else:
            print("Opcion no valida. Por favor, seleccione una opcion entre 1 y 9.")

if __name__ == "__main__":
    main_menu()

