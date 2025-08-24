# Kantenkontraktionsalgorithmen für das Multicut-Problem

Dieses Repository enthält die Implementierungen der in der Bachelor-Arbeit verwendeten Algorithmen.

## Inhalt

- **`greedy_joining_extended.hxx`**: Implementierung des Greedy-Joining-Algorithmus sowie der Erweiterung zur parallelen Kantenkontraktion. Diese Erweiterung berechnet die Kontraktionsmenge mittels des Luby-Jones-Handshaking-Algorithmus zur Berechnung eines maximalen Matchings (mode != 'f') bzw. des Mutex-Kruskal-Algorithmus zur Konstruktion eines konfliktfreien Spannwaldes (mode = 'f') auf ungerichteten, gewichteten Graphen.  
- **`partition.hxx`**: Datenstruktur zur Partitionierung von Knoten, verwendet in den Algorithmen.  

## Nutzung

Die Implementierungen sind in C++ geschrieben und können in Projekte eingebunden werden, die C++17 oder höher unterstützen. 
