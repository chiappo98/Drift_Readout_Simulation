# Drift Readout Simulation

Simulazione del drift e del segnale indotto con Garfield++ e ROOT.

## Requisiti

- CMake
- Garfield++
- ROOT con `ROOTDataFrame`
- PyROOT
- NumPy
- `nlohmann_json` (versione 3.2 o successiva)

## Esecuzione

Da Linux o WSL, con Garfield++ e ROOT configurati:

```sh
bash setup.sh
bash run.sh 01 simulation
bash run.sh 01 analysis
bash run.sh 01 waveform --event 0
bash run.sh 01 raw-signal --event 0
```

`setup.sh` crea le directory necessarie; `run.sh` lo richiama automaticamente. Il numero del run è obbligatorio e viene usato nei nomi `results/runXX.*`. Se si rilancia un numero già presente, i file corrispondenti vengono spostati in una sottocartella datata di `results/backup` prima della simulazione. La modalità `simulation` configura CMake, compila `drift` e lo avvia. Per passare JSON personalizzati: `bash run.sh 01 simulation config/config.json config/simulation.json`.

Ogni simulazione salva in `results/runXX_configurations.txt` il contenuto dei due JSON effettivamente utilizzati (`config.json` e `simulation.json`), con intestazioni che li distinguono.

Oltre a `results/runXX.root`, per impostazione predefinita la simulazione scrive `results/runXX_waveform.root` e `results/runXX_raw_signal.root`. Il TTree `waveforms` contiene una entry per evento, con indice `event` e vettore `float` di campioni `waveform` dopo lo shaper. Il TTree `raw_signals` nel file raw contiene una entry per evento, con indice `event` e vettore `float` di campioni `raw_signal` della corrente prima dello shaper. Entrambi i file sono scritti in modo incrementale, compresso e suddiviso in cluster da circa 30 MB. Si leggono con `ROOT::RDataFrame("waveforms", "results/runXX_waveform.root")` e `ROOT::RDataFrame("raw_signals", "results/runXX_raw_signal.root")`.

Nel TTree `results` di `results/runXX.root`, ogni evento include `totalCharge` come integrale del valore assoluto della corrente campionata, `nPrimaries` (numero di cluster primari) e i vettori `clusterX`, `clusterY`, `clusterZ`, `clusterT`, `clusterE` e `clusterElectrons`. Ogni elemento corrisponde allo stesso cluster nell'ordine restituito da Garfield++; `clusterElectrons` contiene il numero di elettroni del cluster.

Quando sono attivi i grafici dei segnali, il primo pannello di ogni evento mostra anche un cerchio nella posizione x-y di ciascun cluster; il raggio è proporzionale al numero di elettroni rispetto al cluster più popolato dell'evento.

La modalità `analysis` legge il ROOT del run e crea `results/runXX_histograms.root` e `results/runXX_histograms.png`.
Gli istogrammi includono anche la distanza 3D e la differenza temporale assoluta tra cluster consecutivi, il numero di elettroni per cluster e l'energia per cluster.
La stessa analisi crea inoltre `results/runXX_correlations.root` e `results/runXX_correlations.png`, con i grafici di correlazione tra carica totale, ampiezza di picco, rapporto tra ampiezza di picco e carica totale, ToT, tempo di salita e posizione y della traccia.

La modalità `waveform` esegue `analyse_waveform.py` con `python3` sul file `results/runXX_waveform.root`. Per visualizzare un evento specifico:

```sh
bash run.sh 01 waveform --event 0
bash run.sh 01 waveform --overlay
```

Senza opzioni lo script salva il grafico dell'evento 0; `--event` seleziona un altro indice da zero. Con `--overlay` seleziona tutte le waveform con un solo massimo locale di ampiezza almeno il 10% del massimo assoluto, allinea il picco a `t = 0` e normalizza l'ampiezza assoluta a 1, quindi salva tutte le curve sovrapposte in un PDF. `--output percorso.pdf` permette di scegliere il nome del file in entrambe le modalità; le opzioni si possono passare attraverso `run.sh`.

La modalità `raw-signal` esegue `analyse_raw_signal.py` sul file `results/runXX_raw_signal.root` e salva il grafico della corrente raw dell'evento 0 in PDF. `--event` seleziona un evento diverso e `--output percorso.pdf` permette di scegliere il file di output. Con `--fourier` calcola lo spettro di ampiezza monolatero di ogni evento, li somma e salva un unico grafico in `results/runXX_raw_signal_fourier_spectrum.pdf`:

```sh
bash run.sh 01 raw-signal --fourier
```

La frequenza è espressa in GHz (campionamento a 1 ns); l'ampiezza rappresentata è la somma degli spettri dei singoli eventi.

## Configurazione

`config/config.json` contiene shaper, gain, saturazione ASIC, tensioni, geometria e percentuali del gas. `asic_saturation_mV` imposta la saturazione ASIC (600 mV di default): limita il modulo del picco salvato come `q_peak` a `asic_saturation_mV / gain`, mantenendone il segno; non modifica i campioni della waveform. `config/simulation.json` controlla le modalità di simulazione; `signalInduction` include numero di eventi, limiti `y_min`/`y_max`, opzioni di stampa e grafici e i booleani `save_waveforms` e `save_raw_signals` per abilitare o disabilitare separatamente la scrittura dei rispettivi file ROOT. I limiti sono frazioni tra 0 e 1 della semilarghezza della cella; valori uguali fissano una posizione costante.

Il parametro `shaper` accetta `0` (VTX), `1` o `3`. Lo shaper VTX usa `transferVTX`: per `t > t0` calcola `gain * ((t - t0) / (n * tau))^n * exp(n - (t - t0) / tau)`, altrimenti restituisce zero. I parametri sono `gain = -39.6/40 mV/fC`, `n = 3.34`, `tau = 3.06 ns` e `t0 = 1.463 ns`. L'uscita è in mV: imposta `thr_rise` e `thr_fall` in mV. `shaping_time_ns` si applica solo agli shaper 1 e 3, le cui soglie restano in fC.

La macro ROOT `utils/plot_transfer_functions.C` sovrappone le risposte VTX e degli shaper di primo e terzo ordine, normalizzate al rispettivo picco per confrontarne la forma. Il tempo di shaping predefinito è 25 ns; per modificarlo e generare `plot/transfer_functions.png`:

```sh
root -l -q 'utils/plot_transfer_functions.C(50.0)'
```

La macro `utils/simulate_shaper3_square_pulse.C` simula 10.000 eventi, ciascuno con un treno di 400 impulsi quadri da 0,4 µA e 1 ns. Per ogni evento, il ritardo `delayNs` è uniforme tra 0 e 100 ns e i tempi di arrivo dei singoli impulsi sono uniformi tra `delayNs` e `delayNs + TimeSpreadNs`. Convolge il treno con lo shaper di terzo ordine (tempo di shaping predefinito 25 ns e guadagno 2 mV/fC), misura il primo attraversamento della soglia di 1 mV e lo rappresenta in funzione della posizione `delayNs × 50 µm/ns`. `TimeSpreadNs` è il terzo argomento della macro e vale 100 ns per default. Salva il grafico PNG e un file ROOT con nomi `plot/DTrelation_tauXX_spreadYY` (dove `XX` e `YY` sono tempo di shaping e TimeSpread in ns, con `p` per eventuali decimali); il ROOT contiene il grafico, la waveform del primo evento e i dati dei 10.000 eventi:

```sh
root -l -q 'utils/simulate_shaper3_square_pulse.C()'
root -l -q 'utils/simulate_shaper3_square_pulse.C(25.0, 2.0, 50.0)'
```

I percorsi dei file `.gas` e delle mobilità ioniche in `drift.C` devono essere disponibili e, se necessario, adattati all'ambiente locale.