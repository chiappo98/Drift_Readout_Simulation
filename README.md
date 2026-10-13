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

L'eseguibile legge `config/config.json` e `config/simulation.json`. In alternativa, si possono indicare entrambi i percorsi come argomenti. I risultati del segnale vengono salvati in `results/results.root`.

Analisi degli istogrammi con ROOT: `root -l -q 'analyse_results.C'`. Il macro crea `results_histograms.root` e `results_histograms.png`; il file di input si può passare come primo argomento.

## Configurazione

`config/config.json` contiene shaper, tensioni, geometria e percentuali del gas. `config/simulation.json` controlla le modalità di simulazione; `signalInduction` include numero di eventi, limiti `y_min`/`y_max` e opzioni di stampa e grafici. I limiti sono frazioni tra 0 e 1 della semilarghezza della cella; valori uguali fissano una posizione costante.

I percorsi dei file `.gas` e delle mobilità ioniche in `drift.C` devono essere disponibili e, se necessario, adattati all'ambiente locale.