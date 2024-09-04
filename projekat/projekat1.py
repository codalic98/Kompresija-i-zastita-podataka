import heapq
import collections
import math
import struct

# Pomoćne funkcije

def izracunaj_entropiju(podaci):
    brojac_byteova = collections.Counter(podaci)
    ukupno_byteova = len(podaci)
    entropija = 0
    for broj in brojac_byteova.values():
        verovatnoća = broj / ukupno_byteova
        entropija -= verovatnoća * math.log2(verovatnoća)
    return entropija

def bitovi_u_byteove(bitovi):
    return int(bitovi, 2).to_bytes((len(bitovi) + 7) // 8, byteorder='big')

def byteovi_u_bitove(podaci):
    return ''.join(f'{byte:08b}' for byte in podaci)

def sacuvaj_kod_tabelu(ime_fajla, kodna_tabela):
    with open(ime_fajla, 'wb') as f:
        for simbol, kod in kodna_tabela.items():
            f.write(struct.pack('B', simbol))
            f.write(struct.pack('B', len(kod)))
            f.write(kod.encode())

def ucitaj_kod_tabelu(ime_fajla):
    kodna_tabela = {}
    with open(ime_fajla, 'rb') as f:
        while True:
            simbol = f.read(1)
            if not simbol:
                break
            simbol = struct.unpack('B', simbol)[0]
            duzina = struct.unpack('B', f.read(1))[0]
            kod = f.read(duzina).decode()
            kodna_tabela[simbol] = kod
    return kodna_tabela

def ispis_kod_tabele(kodna_tabela):
    for simbol, kod in sorted(kodna_tabela.items()):
        print(f'Kod {simbol}: {kod}')

# Shannon-Fano Kodiranje

class Simbol:
    def __init__(self, simbol, verovatnoca):
        self.simbol = simbol
        self.verovatnoca = verovatnoca
        self.kod = ''

def SF_code(simboli):
    if len(simboli) <= 1:
        return
    
    trenutno = 0
    razlika = []
    suma = sum(simbol.verovatnoca for simbol in simboli)
    
    simboli = sorted(simboli, key=lambda x: x.verovatnoca, reverse=True)
    
    for i, simbol in enumerate(simboli):
        trenutno += simbol.verovatnoca
        ostatak = suma - trenutno
        razlika.append(abs(trenutno - ostatak))
    
    targetIndex = razlika.index(min(razlika))
    
    for i, simbol in enumerate(simboli):
        simbol.kod += "0" if i <= targetIndex else "1"

    SF_code(simboli[:targetIndex + 1])
    SF_code(simboli[targetIndex + 1:])


def findElement(lista):
    return lista.index(min(lista))

def shanon_fano_kodiranje(podaci):
    frekvencije = collections.Counter(podaci)
    ukupno = len(podaci)
    simboli = [Simbol(byte, count / ukupno) for byte, count in frekvencije.items()]
    SF_code(simboli)
    return {simbol.simbol: simbol.kod for simbol in simboli}

def shanon_fano_dekodiraj(kodirani_podaci, kodna_tabela):
    bitovi = byteovi_u_bitove(kodirani_podaci)
    dekodirani = bytearray()
    bafer = ''
    kodovi = {v: k for k, v in kodna_tabela.items()}
    for bit in bitovi:
        bafer += bit
        if bafer in kodovi:
            dekodirani.append(kodovi[bafer])
            bafer = ''
    return dekodirani

# Huffman Kodiranje

class Cvor:
    def __init__(self, simbol, frekvencija):
        self.simbol = simbol
        self.frekvencija = frekvencija
        self.levo = None
        self.desno = None
    
    def __lt__(self, drugi):
        return self.frekvencija < drugi.frekvencija

def huffman_kodiranje(podaci):
    frekvencije = collections.Counter(podaci)
    prioritetna_reda = [Cvor(simbol, frekvencija) for simbol, frekvencija in frekvencije.items()]
    heapq.heapify(prioritetna_reda)
    
    while len(prioritetna_reda) > 1:
        levo = heapq.heappop(prioritetna_reda)
        desno = heapq.heappop(prioritetna_reda)
        spojeni = Cvor(None, levo.frekvencija + desno.frekvencija)
        spojeni.levo = levo
        spojeni.desno = desno
        heapq.heappush(prioritetna_reda, spojeni)
    
    def napravi_tabelu_kodova(cvor, kod=''):
        if cvor is None:
            return {}
        if cvor.simbol is not None:
            return {cvor.simbol: kod}
        levo_tabela = napravi_tabelu_kodova(cvor.levo, kod + '0')
        desno_tabela = napravi_tabelu_kodova(cvor.desno, kod + '1')
        return {**levo_tabela, **desno_tabela}
    
    huffman_stablo = prioritetna_reda[0]
    return napravi_tabelu_kodova(huffman_stablo)

def huffman_dekodiraj(kodirani_podaci, kodna_tabela):
    bitovi = byteovi_u_bitove(kodirani_podaci)
    dekodirani = bytearray()
    bafer = ''
    kodovi = {v: k for k, v in kodna_tabela.items()}
    for bit in bitovi:
        bafer += bit
        if bafer in kodovi:
            dekodirani.append(kodovi[bafer])
            bafer = ''
    return dekodirani

# LZ77 Kompresija

def lz77_kompresuj(podaci, velicina_prozora=256, velicina_pregleda=16):
    i = 0
    kompresovani = []
    while i < len(podaci):
        duzina_podudaranja = 0
        pomak_podudaranja = 0

        pocetak = max(0, i - velicina_prozora)
        kraj_pregleda = min(len(podaci), i + velicina_pregleda)

        for j in range(pocetak, i):
            duzina = 0
            while (i + duzina < kraj_pregleda) and (podaci[j + duzina] == podaci[i + duzina]):
                duzina += 1
                if duzina > velicina_pregleda:
                    break

            if duzina > duzina_podudaranja:
                duzina_podudaranja = duzina
                pomak_podudaranja = i - j

        if duzina_podudaranja > 0:
            kompresovani.append((pomak_podudaranja, duzina_podudaranja, podaci[i + duzina_podudaranja] if (i + duzina_podudaranja) < len(podaci) else 0))
            i += duzina_podudaranja + 1
        else:
            kompresovani.append((0, 0, podaci[i]))
            i += 1
    
    return kompresovani

def lz77_dekompresuj(kompresovani, velicina_prozora=256):
    dekompresovani = bytearray()
    for pomak, duzina, sledeci_byte in kompresovani:
        start = len(dekompresovani) - pomak
        dekompresovani.extend(dekompresovani[start:start + duzina])
        dekompresovani.append(sledeci_byte)
    return dekompresovani

# LZW Kompresija

def lzw_kompresuj(podaci):
    rečnik = {bytes([i]): i for i in range(256)}
    w = bytes()
    kompresovani = []
    velicina_recnika = 256
    for byte in podaci:
        wc = w + bytes([byte])
        if wc in rečnik:
            w = wc
        else:
            kompresovani.append(rečnik[w])
            rečnik[wc] = velicina_recnika
            velicina_recnika += 1
            w = bytes([byte])
    if w:
        kompresovani.append(rečnik[w])
    return kompresovani

def lzw_dekompresuj(kompresovani):
    rečnik = {i: bytes([i]) for i in range(256)}
    w = bytes([kompresovani.pop(0)])
    dekompresovani = bytearray(w)
    velicina_recnika = 256
    for kôd in kompresovani:
        if kôd in rečnik:
            entry = rečnik[kôd]
        elif kôd == velicina_recnika:
            entry = w + w[0:1]
        else:
            raise ValueError('Kôd nije pronađen u rečniku.')
        dekompresovani.extend(entry)
        rečnik[velicina_recnika] = w + entry[0:1]
        velicina_recnika += 1
        w = entry
    return dekompresovani

# Glavna Funkcija za Testiranje

def glavni():
    ime_fajla = 'primer.bin'
    
    with open(ime_fajla, 'rb') as f:
        podaci = f.read()
    
    print("Byte Entropija:", izracunaj_entropiju(podaci))
    
    # Shannon-Fano Kodiranje
    shanon_fano_stablo = shanon_fano_kodiranje(podaci)
    shanon_fano_kodirani = ''.join(shanon_fano_stablo[byte] for byte in podaci)
    shanon_fano_kodirani_byteovi = bitovi_u_byteove(shanon_fano_kodirani)
    
    sacuvaj_kod_tabelu('shannon_fano_tabela.bin', shanon_fano_stablo)
    with open('shannon_fano_kodirani.bin', 'wb') as f:
        f.write(shanon_fano_kodirani_byteovi)
    
    # Ispis Shannon-Fano kodne tabele
    print("Shannon-Fano kodna tabela:")
    ispis_kod_tabele(shanon_fano_stablo)
    
    # Shannon-Fano Dekodiranje
    kod_tabela = ucitaj_kod_tabelu('shannon_fano_tabela.bin')
    with open('shannon_fano_kodirani.bin', 'rb') as f:
        shanon_fano_kodirani_podaci = bytes(f.read())
        shanon_fano_dekodirani = shanon_fano_dekodiraj(shanon_fano_kodirani_podaci, kod_tabela)
        with open('shannon_fano_dekodirani.bin', 'wb') as f:
            f.write(shanon_fano_dekodirani)
    
    if shanon_fano_dekodirani == podaci:
        print("Shannon-Fano kodiranje i dekodiranje uspešno.")
    else:
        print("Shannon-Fano kodiranje i dekodiranje nije uspelo.")
    
    # Huffman Kodiranje
    huffman_stablo = huffman_kodiranje(podaci)
    huffman_kodirani = ''.join(huffman_stablo[byte] for byte in podaci)
    huffman_kodirani_byteovi = bitovi_u_byteove(huffman_kodirani)
    
    sacuvaj_kod_tabelu('huffman_tabela.bin', huffman_stablo)
    with open('huffman_kodirani.bin', 'wb') as f:
        f.write(huffman_kodirani_byteovi)
    
    # Ispis Huffman kodne tabele
    print("Huffman kodna tabela:")
    ispis_kod_tabele(huffman_stablo)
    
    # Huffman Dekodiranje
    kod_tabela = ucitaj_kod_tabelu('huffman_tabela.bin')
    with open('huffman_kodirani.bin', 'rb') as f:
        huffman_kodirani_podaci = bytes(f.read())
        huffman_dekodirani = huffman_dekodiraj(huffman_kodirani_podaci, kod_tabela)
        with open('huffman_dekodirani.bin', 'wb') as f:
            f.write(huffman_dekodirani)
    
    if huffman_dekodirani == podaci:
        print("Huffman kodiranje i dekodiranje uspešno.")
    else:
        print("Huffman kodiranje i dekodiranje nije uspelo.")
    
    # LZ77 Kompresija
    lz77_kompresovani = lz77_kompresuj(podaci)
    with open('lz77_kompresovani.bin', 'wb') as f:
        for pomak, duzina, sledeci_byte in lz77_kompresovani:
            f.write(struct.pack('>H', pomak))
            f.write(struct.pack('>H', duzina))
            f.write(struct.pack('B', sledeci_byte))
    
    # LZ77 Dekompresija
    with open('lz77_kompresovani.bin', 'rb') as f:
        lz77_kompresovani = []
        while True:
            pomak = f.read(2)
            duzina = f.read(2)
            sledeci_byte = f.read(1)
            if not pomak or not duzina or not sledeci_byte:
                break
            pomak = struct.unpack('>H', pomak)[0]
            duzina = struct.unpack('>H', duzina)[0]
            sledeci_byte = struct.unpack('B', sledeci_byte)[0]
            lz77_kompresovani.append((pomak, duzina, sledeci_byte))
        lz77_dekodirani = lz77_dekompresuj(lz77_kompresovani)
        with open('lz77_dekodirani.bin', 'wb') as f:
            f.write(lz77_dekodirani)
    
    if lz77_dekodirani == podaci:
        print("LZ77 kompresija i dekompresija uspešna.")
    else:
        print("LZ77 kompresija i dekompresija nije uspela.")
    
    # LZW Kompresija
    lzw_kompresovani = lzw_kompresuj(podaci)
    with open('lzw_kompresovani.bin', 'wb') as f:
        for kôd in lzw_kompresovani:
            f.write(struct.pack('>H', kôd))
    
    # LZW Dekompresija
    with open('lzw_kompresovani.bin', 'rb') as f:
        lzw_kompresovani = []
        while True:
            kôd = f.read(2)
            if not kôd:
                break
            lzw_kompresovani.append(struct.unpack('>H', kôd)[0])
        lzw_dekodirani = lzw_dekompresuj(lzw_kompresovani)
        with open('lzw_dekodirani.bin', 'wb') as f:
            f.write(lzw_dekodirani)
    
    if lzw_dekodirani == podaci:
        print("LZW kompresija i dekompresija uspešna.")
    else:
        print("LZW kompresija i dekompresija nije uspela.")

if __name__ == "__main__":
    glavni()
