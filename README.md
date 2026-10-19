# Drift Readout Simulation

Simulazione del drift e del segnale indotto con Garfield++ e ROOT.

## Requisiti

- CMake
- Garfield++
- ROOT con `ROOTDataFrame`
- PyROOT
- `nlohmann_json` (versione 3.2 o successiva)

## Esecuzione

Da Linux o WSL, con Garfield++ e ROOT configurati:

```sh
bash setup.sh
bash run.sh 01 simulation
bash run.sh 01 analysis
bash run.sh 01 waveform --event 0
```

`setup.sh` crea le directory necessarie; `run.sh` lo richiama automaticamente. Il numero del run è obbligatorio e viene usato nei nomi `results/runXX.*`. Se si rilancia un numero già presente, i file corrispondenti vengono spostati in una sottocartella datata di `results/backup` prima della simulazione. La modalità `simulation` configura CMake, compila `drift` e lo avvia. Per passare JSON personalizzati: `bash run.sh 01 simulation config/config.json config/simulation.json`.

Oltre a `results/runXX.root`, la simulazione scrive `results/runXX_waveform.root`: il TTree `waveforms` contiene una entry per evento, con indice `event` e vettore `float` di campioni `waveform`. La scrittura è incrementale, compressa e suddivisa in cluster da circa 30 MB. Si legge con `ROOT::RDataFrame("waveforms", "results/runXX_waveform.root")`.

Nel TTree `results` di `results/runXX.root`, ogni evento include `nPrimaries` (numero di cluster primari) e i vettori `clusterX`, `clusterY`, `clusterZ`, `clusterT`, `clusterE` e `clusterElectrons`. Ogni elemento corrisponde allo stesso cluster nell'ordine restituito da Garfield++; `clusterElectrons` contiene il numero di elettroni del cluster.

Quando sono attivi i grafici dei segnali, il primo pannello di ogni evento mostra anche un cerchio nella posizione x-y di ciascun cluster; il raggio è proporzionale al numero di elettroni rispetto al cluster più popolato dell'evento.

La modalità `analysis` legge il ROOT del run e crea `results/runXX_histograms.root` e `results/runXX_histograms.png`.
Gli istogrammi includono anche la distanza 3D e la differenza temporale assoluta tra cluster consecutivi, il numero di elettroni per cluster e l'energia per cluster.
La stessa analisi crea inoltre `results/runXX_correlations.root` e `results/runXX_correlations.png`, con i grafici di correlazione tra carica totale, ampiezza di picco, durate temporali e posizione y della traccia.

La modalità `waveform` esegue `analyse_waveform.py` con `python3` sul file `results/runXX_waveform.root`. Per visualizzare un evento specifico:

```sh
bash run.sh 01 waveform --event 0
bash run.sh 01 waveform --overlay
```

Senza opzioni lo script salva il grafico dell'evento 0; `--event` seleziona un altro indice da zero. Con `--overlay` seleziona tutte le waveform con un solo massimo locale di ampiezza almeno il 10% del massimo assoluto, allinea il picco a `t = 0` e normalizza l'ampiezza assoluta a 1, quindi salva tutte le curve sovrapposte in un PDF. `--output percorso.pdf` permette di scegliere il nome del file in entrambe le modalità; le opzioni si possono passare attraverso `run.sh`.

## Configurazione

`config/config.json` contiene shaper, tensioni, geometria e percentuali del gas. `config/simulation.json` controlla le modalità di simulazione; `signalInduction` include numero di eventi, limiti `y_min`/`y_max` e opzioni di stampa e grafici. I limiti sono frazioni tra 0 e 1 della semilarghezza della cella; valori uguali fissano una posizione costante.

Il parametro `shaper` accetta `0` (VTX), `1` o `3`. Lo shaper VTX usa `transferVTX`: per `t > t0` calcola `gain * ((t - t0) / (n * tau))^n * exp(n - (t - t0) / tau)`, altrimenti restituisce zero. I parametri sono `gain = -39.6/40 mV/fC`, `n = 3.34`, `tau = 3.06 ns` e `t0 = 1.463 ns`. L'uscita è in mV: imposta `thr_rise` e `thr_fall` in mV. `shaping_time_ns` si applica solo agli shaper 1 e 3, le cui soglie restano in fC.

I percorsi dei file `.gas` e delle mobilità ioniche in `drift.C` devono essere disponibili e, se necessario, adattati all'ambiente locale.