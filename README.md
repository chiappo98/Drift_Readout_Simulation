# Drift Readout Simulation

Simulazione del drift e del segnale indotto con Garfield++ e ROOT.

## Requisiti

- CMake
- Garfield++
- ROOT con `ROOTDataFrame`
- `nlohmann_json` (versione 3.2 o successiva)

## Build ed esecuzione

Con l'ambiente Garfield++ e ROOT configurato:

```sh
cmake -S . -B build
cmake --build build
cd build
./drift
```

L'eseguibile legge `config/config.json` e `config/simulation.json`. In alternativa, si possono indicare entrambi i percorsi come argomenti. I risultati del segnale vengono salvati in `results.root`.

## Configurazione

`config/config.json` contiene shaper, tensioni, geometria e percentuali del gas. `config/simulation.json` controlla le modalità di simulazione; `signalInduction` include numero di eventi, posizione e opzioni di stampa e grafici. `position` è una frazione tra 0 e 1.

I percorsi dei file `.gas` e delle mobilità ioniche in `drift.C` devono essere disponibili e, se necessario, adattati all'ambiente locale.