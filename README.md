# Drift Readout Simulation

Simulazione del drift e del segnale indotto con Garfield++ e ROOT.

## Requisiti

- CMake
- Garfield++
- ROOT con `ROOTDataFrame`
- `nlohmann_json` (versione 3.2 o successiva)

## Esecuzione

Da Linux o WSL, con Garfield++ e ROOT configurati:

```sh
bash setup.sh
bash run.sh 01 simulation
bash run.sh 01 analysis
```

`setup.sh` crea le directory necessarie; `run.sh` lo richiama automaticamente. Il numero del run è obbligatorio e viene usato nei nomi `results/runXX.*`. Se si rilancia un numero già presente, i file corrispondenti vengono spostati in una sottocartella datata di `results/backup` prima della simulazione. La modalità `simulation` configura CMake, compila `drift` e lo avvia. Per passare JSON personalizzati: `bash run.sh 01 simulation config/config.json config/simulation.json`.

La modalità `analysis` legge il ROOT del run e crea `results/runXX_histograms.root` e `results/runXX_histograms.png`.

## Configurazione

`config/config.json` contiene shaper, tensioni, geometria e percentuali del gas. `config/simulation.json` controlla le modalità di simulazione; `signalInduction` include numero di eventi, limiti `y_min`/`y_max` e opzioni di stampa e grafici. I limiti sono frazioni tra 0 e 1 della semilarghezza della cella; valori uguali fissano una posizione costante.

I percorsi dei file `.gas` e delle mobilità ioniche in `drift.C` devono essere disponibili e, se necessario, adattati all'ambiente locale.