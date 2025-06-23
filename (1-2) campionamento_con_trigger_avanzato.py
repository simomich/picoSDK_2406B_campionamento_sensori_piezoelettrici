#!/usr/bin/env python
# coding


# author:                  Simone Micelli
# email:        simonemicelli47@gmail.com
# date:                        05/05/2025
# description:
#     Campionamento con trigger avanzato per oscilloscopio PicoScope 2406B
#     Utilizza la libreria picosdk per il controllo dei driver.
#     La libreria picosdk è disponibile su GitHub:

# Copyright (CC) 2025 Simone Micelli



import ctypes
from picosdk.ps2000a import ps2000a as ps
from picosdk.functions import mV2adc, adc2mV, assert_pico_ok
import numpy as np

#! DEFINIZIONE DELLE COSTANTI
SEMI_INTERVALLO_TENSIONE_INT = 20          #[V] (±20V)
SOGLIA_TRIGGER = 2000                      #[mV] Soglia del trigger
TOTAL_SAMPLES = 195000                     #* Numero totale di campioni minimo
PRE_TRIGGER_PERCENT = 0.3                  #* Percentuale di campioni pre-trigger (30% del totale)
DELTA_T_CAMPIONAMENTO = 50E6               #[ns] ∆t di campionamento voluto: 50 [ms] = 50E6 [ns]
SOGLIA_ISTERESI_TRIGGER = 0.05             #* Percentuale di isteresi del trigger (5% del totale)


#* Creaiamo un handle per il dispositivo
#* chandle sarà utilizzato nei metodi dove necessario handle come parametro
#* status invece è un dict che riceverà le info sullo stato dell'oscilloscopio
status = {}
chandle = ctypes.c_int16()

#* Nessun downsampling
DOWNSAMPLING_MODE = ps.PS2000A_RATIO_MODE["PS2000A_RATIO_MODE_NONE"]

#? Avvia il dispositivo + 1 flash del led
status["openunit"] = ps.ps2000aOpenUnit(ctypes.byref(chandle), None)

#? Blocco di codice per la gestione dell'errore
try:
    assert_pico_ok(status["openunit"])
except:
    # powerstate becomes the status number of openunit
    powerstate = status["openunit"]

    # If powerstate is the same as 282 then it will run this if statement
    if powerstate == 282:
        # Changes the power input to "PICO_POWER_SUPPLY_NOT_CONNECTED"
        status["ChangePowerSource"] = ps.ps2000aChangePowerSource(chandle, 282)
    # If the powerstate is the same as 286 then it will run this if statement
    elif powerstate == 286:
        # Changes the power input to "PICO_USB3_0_DEVICE_NON_USB3_0_PORT"
        status["ChangePowerSource"] = ps.ps2000aChangePowerSource(chandle, 286)
    else:
        raise

    assert_pico_ok(status["ChangePowerSource"])



INTERVALLO_TENSIONE = ps.PS2000A_RANGE["PS2000A_" + str(SEMI_INTERVALLO_TENSIONE_INT) + "V"] #? Corrispondente a ±20[V]

attivo = 1
disattivo = 0
offset_analogico = 0.0

#* SetChannel dei canali A, B, C e D
#? ps.PS2000A_COUPLING["PS2000A_DC"] Accoppiamento DC:
#? 1 MOhm di impedenza; il canale accetta tutte le frequenze in input
#? da zero (DC) fino al massimo analogico della larghezza di banda

for chan in ["A", "B", "C", "D"]:
    status["setCh" + chan] = ps.ps2000aSetChannel(chandle,
                                                  ps.PS2000A_CHANNEL["PS2000A_CHANNEL_" + chan],
                                                  attivo,
                                                  ps.PS2000A_COUPLING["PS2000A_DC"],
                                                  INTERVALLO_TENSIONE,
                                                  offset_analogico)

    assert_pico_ok(status["setCh" + chan])



#* Troviamo il valore massimo di ADC COUNTS (Vedi guida alla programmazione dell'oscilloscopio {Voltage Range})
#? Puntatore ctypes.byref(maxADC) alla variabile maxADC di tipo c_int16
#? maxADC è un valore che rappresenta il massimo valore ADC COUNTS
maxADC = ctypes.c_int16()
status["maximumValue"] = ps.ps2000aMaximumValue(chandle, ctypes.byref(maxADC))
print("Massimo valore ADC COUNTS (Codifica della tensione):", maxADC.value)
assert_pico_ok(status["maximumValue"])


soglia = SOGLIA_TRIGGER
soglia_adc = mV2adc(SOGLIA_TRIGGER, INTERVALLO_TENSIONE, maxADC)
print("Soglia di trigger in counts adc: ", soglia_adc)
autotrigger_ms = 0

#? Modalità che evita un nuovo trigger se il valore non scende sotto la soglia di isteresi.
#? In questo caso non è utile, ma in altri casi in cui il campionamento è continuo allora è molto importante.
percIsteresi = SOGLIA_ISTERESI_TRIGGER

condizione_attiva = ps.PS2000A_TRIGGER_STATE["PS2000A_CONDITION_TRUE"]
condizione_ignorata = ps.PS2000A_TRIGGER_STATE["PS2000A_CONDITION_DONT_CARE"]

condizioni_ch = [0, 1, 2, 3]
num_condizioni = len(condizioni_ch) #* Lunghezza dell'array che sarà passato alla funzione ps2000aSetTriggerChannelConditions

for ch in condizioni_ch:
    #! Ogni struct, condizione del trigger, effettua un AND logico
    lista_condizioni = [condizione_ignorata] * 8
    lista_condizioni[ch] = condizione_attiva        #? Posizione 0: canale A, 1: canale B, 2: canale C, 3: canale D

    #? La funzione PS2000A_TRIGGER_CONDITIONS accetta le condizioni del trigger come argomento (AND logico).
    #? e restituisce un oggetto PS2000A_TRIGGER_CONDITIONS che rappresenta le condizioni del trigger per il canale specificato.
    #? In particolare noi decomponiamo la lista delle condizioni negli 8 parametri richiesti dalla funzione.
    condizioni_ch[ch] = ps.PS2000A_TRIGGER_CONDITIONS(*lista_condizioni)

#* Creazione dell'array degli oggetti PS2000A_TRIGGER_CONDITIONS
arr_condizioni_ch = (ps.PS2000A_TRIGGER_CONDITIONS * len(condizioni_ch))(*condizioni_ch)

status["setConditions"] = ps.ps2000aSetTriggerChannelConditions(chandle, ctypes.byref(arr_condizioni_ch), num_condizioni)
assert_pico_ok(status["setConditions"])


d_rising = ps.PS2000A_THRESHOLD_DIRECTION["PS2000A_RISING"]
non_set = ps.PS2000A_THRESHOLD_DIRECTION["PS2000A_NONE"]

status["setDirections"] = ps.ps2000aSetTriggerChannelDirections(chandle, d_rising, #* Canale A
                                                                         d_rising, #* Canale B
                                                                         d_rising, #* Canale C
                                                                         d_rising, #* Canale D
                                                                         non_set, 
                                                                         non_set)
assert_pico_ok(status["setDirections"])


prop_ch = []

#* Impostiamo il trigger per ogni canale con le relative soglie [mV in counts adc]
for chan in ["A", "B", "C", "D"]:
    properties = ps.PS2000A_TRIGGER_CHANNEL_PROPERTIES(soglia_adc,                                     #* Soglia del trigger in counts adc
                                                      int(soglia_adc * percIsteresi),                  #* Soglia di isteresi in counts adc
                                                      (soglia_adc * -1),                               #* Soglia di trigger negativa in counts adc
                                                      int(soglia_adc * percIsteresi),                  #* Soglia di isteresi negativa in counts adc
                                                      ps.PS2000A_CHANNEL["PS2000A_CHANNEL_" + chan],   #* Canale
                                                      ps.PS2000A_THRESHOLD_MODE["PS2000A_LEVEL"])      #* Modalità di soglia del trigger

    prop_ch.append(properties)

num_prop = len(prop_ch)

#* Creazione dell'array degli oggetti PS2000A_TRIGGER_CHANNEL_PROPERTIES
arr_prop_ch = (ps.PS2000A_TRIGGER_CHANNEL_PROPERTIES * len(prop_ch))(*prop_ch)

status["setProperties"] = ps.ps2000aSetTriggerChannelProperties(chandle, ctypes.byref(arr_prop_ch), num_prop, 0, autotrigger_ms)
assert_pico_ok(status["setProperties"])



#* Impostiamo il tempo di campionamento
totalSamples = TOTAL_SAMPLES

delta_t_camp = DELTA_T_CAMPIONAMENTO
#? int() tronca il valore a intero verso lo zero, buono per il nostro caso dove il tempo
#? di campionamento deve essere un numero intero e diventa più piccolo dell'intervallo teorico necessario calcolato
timeInterval_ns_obbiettivo = delta_t_camp / (totalSamples - 1)

#* Ora utilizziamo la funzione che restituisce il timeInterval_ns per la data timebase
#* ed il massimo valore del numero totale dei campioni, in questo caso il 2406B con 4 canali attivi

#? Ricordiamo che la timebase (o base dei tempi) è il nome del circuito che genera
#? il clock di campionamento dell'oscilloscopio.
timebase = 0
if timeInterval_ns_obbiettivo >= 8:
    timebase = int(0.125 * timeInterval_ns_obbiettivo + 2)
else:
    timebase = int(np.log2(timeInterval_ns_obbiettivo))

timeInterval_ns = ctypes.c_float()
returnedMaxSamples = ctypes.c_int32()

#? oversample è utilizzato dalla tecnica Resolution enhancement che aumenta la risoluzione
#? verticale del campionamento a discapito della velocità massima di campionamento disponibile.
#? Se lavoriamo con segnali a sufficiente bassa frequenza possiamo utilizzare questa tecnica
#? per aumentare la risoluzione verticale del campionamento.
#? In questo caso non la utilizziamo (anche se potremmo), quindi impostiamo a zero.
oversample = ctypes.c_int16(0)

segmentIndex = 0 #* Non usiamo la memoria segmentata


#! Nelle funzioni della libreria, in particolare nella seguente, dove vediamo *
#! indicante un puntatore ad una variabile, dobbiamo passare il riferimento alla variabile
#! tramite ctypes.byref() altrimenti non funziona.
#! In molti casi, come questo, il puntatore passato modificherà la variabile originale
#! (in questo caso timeInterval_ns e returnedMaxSamples)


status["getTimebase2"] = ps.ps2000aGetTimebase2(chandle,
                                                timebase,
                                                totalSamples,
                                                ctypes.byref(timeInterval_ns),
                                                oversample,
                                                ctypes.byref(returnedMaxSamples),
                                                segmentIndex)

assert_pico_ok(status["getTimebase2"])


if timeInterval_ns.value > timeInterval_ns_obbiettivo:
    print("Attenzione: il tempo di campionamento non è quello richiesto!")
    print("Richiesta:", timeInterval_ns_obbiettivo, "ns")
    print("Restituito:", timeInterval_ns.value, "ns")

print("Timebase:", timebase)
print("Intervallo di tempo [ns] tra ogni campione:", int(timeInterval_ns.value))
print("Massimo valore totale di campioni con questa configurazione:", returnedMaxSamples.value)



totalSamples = np.ceil(delta_t_camp / timeInterval_ns.value) + 1   #* Val aggiustato
preTriggerPercent = PRE_TRIGGER_PERCENT                            #* 30% del totale
preTriggerSamples = int(totalSamples * preTriggerPercent)          #* Calcolo del numero di campioni pre-trigger
postTriggerSamples = totalSamples - preTriggerSamples              #* 70% del totale


#? Variabile che traccia se il campionamento è terminato o meno
callback_eseguito = False

def callback_campionamento_terminato(handle, statusCode, pParam):

    print("Campionamento terminato")

    global callback_eseguito
    callback_eseguito = True

cFuncPtr = ps.BlockReadyType(callback_campionamento_terminato)

#* Startiamo il campionamento
#? oversample = 0 = oversample non utilizzato
#? time indisposed ms = None = Non necessario
#? segment index = 0 (Non utilizziamo la memoria segmentata)
#? lpReady = cFuncPtr = Funzione di callback
#? pParameter = None = Non necessario attualmente, è utilizzabile per scrivere
#?                     qualunque dato nel callback accessibile poi all'esterno
status["runBlock"] = ps.ps2000aRunBlock(chandle,
                                        preTriggerSamples,
                                        postTriggerSamples,
                                        timebase,
                                        oversample,
                                        None,
                                        segmentIndex,
                                        cFuncPtr,
                                        None)
assert_pico_ok(status["runBlock"])


print("Campionamento in corso (In attesa del trigger)...")


#TODO: Implementare un comportamento asincrono per il check del callback eseguito
#? Check, per l'acquisizione ed il termine del callback, ogni 0.01 secondi
#! ATTUALMENTE BLOCKING THE THREAD
import time
while callback_eseguito == False:
    time.sleep(0.01)



#* dict contenente i buffer per i canali A, B, C e D
#? Ogni buffer è un array di int16 (2 byte VEDI ADC COUNT AC2DC CONVERTER) di lunghezza totalSamples
buffers = {
    "bufferA": (ctypes.c_int16 * totalSamples)(),
    "bufferB": (ctypes.c_int16 * totalSamples)(),
    "bufferC": (ctypes.c_int16 * totalSamples)(),
    "bufferD": (ctypes.c_int16 * totalSamples)()
}

#* Settaggio delle variabili buffer
#? La funzione richiede un puntatore al buffer by reference
#? buffer length = totalSamples (numero di campioni totale)
#? segment index = 0 (Non utilizziamo la memoria segmentata)
#? ratio mode = PS2000A_RATIO_MODE_NONE = 0
for chan in ["A", "B", "C", "D"]:
    status["setDataBuffers" + chan] = ps.ps2000aSetDataBuffer(chandle,
                                                              ps.PS2000A_CHANNEL["PS2000A_CHANNEL_" + chan],
                                                              ctypes.byref(buffers["buffer" + chan]),
                                                              totalSamples,
                                                              segmentIndex,
                                                              DOWNSAMPLING_MODE)
    assert_pico_ok(status["setDataBuffers" + chan])

#* Rappresenta un bit field che indica in quali canali è avvenuto l'OverRange
#? 0000000000000001 = canale A
#? 0000000000000010 = canale B
#? 0000000000000100 = canale C
#? 0000000000001000 = canale D
#? 0000000000000000 = nessun canale
#? 0000000000000101 = canale A e C
#? 0000000000001111 = tutti e quattro i canali
#? E le altre combinazioni possibili
overflow = ctypes.c_int16()



#* Conversione del tipo di totalSamples perché dopo è richiesto by reference
cTotalSamples = ctypes.c_int32(totalSamples)

#? startIndex è l'indice di partenza per il campionamento, in questo caso il primo campione
start_index = 0
downsampling_ratio = 0 #* Non usiamo il downsampling, quindi 0

#* Invio dei dati ai buffer
status["getValues"] = ps.ps2000aGetValues(chandle,
                                          start_index,
                                          ctypes.byref(cTotalSamples),
                                          downsampling_ratio,
                                          DOWNSAMPLING_MODE,
                                          segmentIndex,
                                          ctypes.byref(overflow))
assert_pico_ok(status["getValues"])



#* Convertiamo i valori ADC COUNTS in [V]
#? Espressione di adc2mV modificata (VEDI LINK SOTTOSTANTE):
#? https://github.com/picotech/picosdk-python-wrappers/blob/master/picosdk/functions.py

def adc2V(nome_buffer: str) -> np.array:
    return np.array([(np.int64(x) * SEMI_INTERVALLO_TENSIONE_INT) / maxADC.value for x in buffers[nome_buffer]])

chanA_V = adc2V("bufferA")
chanB_V = adc2V("bufferB")
chanC_V = adc2V("bufferC")
chanD_V = adc2V("bufferD")

#* Creazione del numpy array del tempo di campionamento totale [ms]
tempo_pre_trigger_ns = delta_t_camp * preTriggerPercent
tempo_post_trigger_ns = (totalSamples - 1) * timeInterval_ns.value - tempo_pre_trigger_ns

tempo_camp_ms = np.linspace( -1 * (tempo_pre_trigger_ns) / 1E6, tempo_post_trigger_ns / 1E6, totalSamples)



#* Plot dei dati dei canali A, B, C e D
import matplotlib.pyplot as plt

plt.figure(dpi=150, layout="tight")
plt.plot(tempo_camp_ms, chanA_V, linewidth = "1")
plt.plot(tempo_camp_ms, chanB_V, linewidth = "1")
plt.plot(tempo_camp_ms, chanC_V, linewidth = "1")
plt.plot(tempo_camp_ms, chanD_V, linewidth = "1")
plt.title("Rilevamento Oscilloscopio")
plt.legend(["Canale A", "Canale B", "Canale C", "Canale D"])
plt.grid()
plt.ticklabel_format(axis="both", style="plain")
plt.xlabel("Time [ms]")
plt.ylabel("Voltage [V]")
plt.show()

# Stop the scope
status["stop"] = ps.ps2000aStop(chandle)
assert_pico_ok(status["stop"])

# Close unit / Disconnect the scope
status["close"] = ps.ps2000aCloseUnit(chandle)
assert_pico_ok(status["close"])

# display status returns
# print(status)



#* Informiamo l'utente di un eventuale overrange
lista_canali_overrange = []
if overflow.value != 0:
    #* Convertiamo il valore in binario (16 bit, complemento a due)
    binario = format(overflow.value, "016b")
    print("Overrange:", binario)

    if binario[-1] == "1":
        lista_canali_overrange.append("A")
    if binario[-2] == "1":
        lista_canali_overrange.append("B")
    if binario[-3] == "1":
        lista_canali_overrange.append("C")
    if binario[-4] == "1":
        lista_canali_overrange.append("D")

    print("Attenzione: Overrange di tensione")
    if len(lista_canali_overrange) == 1:
        print("Impatto registrato dal canale:", lista_canali_overrange[0])
    else:
        print("OverRange sui canali", lista_canali_overrange, "\nImpatto molto forte!")


#* Chiediamo all'utente se intende salvare i dati in un file CSV e grafico in PNG
if input("Vuoi salvare i dati in un file CSV? (s/n): ").lower() == "s":

    import pandas as pd

    #* Creazione di un DataFrame con i dati dei canali A, B, C e D
    df = pd.DataFrame({
        "Tempo [ms]":    tempo_camp_ms,
        "Canale A [V]":        chanA_V,
        "Canale B [V]":        chanB_V,
        "Canale C [V]":        chanC_V,
        "Canale D [V]":        chanD_V
    })

    #* Export del DataFrame in CSV
    #? Chiediamo all'utente il nome del file
    file_csv_name = input("Nome del file CSV (senza estensione): ")

    #* Creiamo una cartella per i file CSV se non esiste
    import os
    newpath = "shots"
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    #* Salviamo il DataFrame in CSV nella cartella newpath
    #? Se il file esiste già, lo sovrascriviamo
    df.to_csv(newpath + "/" + file_csv_name + ".csv", index=False, header=True)

    print("File salvato in:", newpath + "/" + file_csv_name + ".csv")
