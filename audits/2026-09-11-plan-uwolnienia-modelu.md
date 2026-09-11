# Plan uwolnienia modelu od narzuconego n=1

Stan na 2026-09-11, po sekcjach 26-29 audytu para/orto.

Cel postawiony przez autora: zdjąć z modelu narzucone założenia, w
szczególności wpisane `n=1`, przy intuicji, że **konfiguracje momentów
magnetycznych powinny same narzucić poziomy energetyczne** — i zrobić to tak,
żeby ocalały wyniki nietautologiczne, zwłaszcza te zgodne z fizyką kwantową
albo tłumaczące efekty kwantowe.

---

## 0. Co musi przeżyć, a co jest tautologią

To jest warunek wstępny całej reszty: bez tej listy „ocalić najcenniejsze"
nie ma treści.

### 0a. Wyniki nietautologiczne — do ochrony regresją

| # | wynik | dlaczego nie jest tautologią |
|---|---|---|
| N1 | \(\mu/(ec)\) **jest** promieniem bariery Comptona, co do 15 cyfr | bariera była wpisana jako „gdzie kończy się elektrodynamika punktowa", a wychodzi jako własna długość konfiguracji momentów; te dwa uzasadnienia są niezależne |
| N2 | kształt reguły Landaua-Yanga z klasycznej magnetostatyki | koherentny kanał M1 liczy się z \(\mathbf m=\boldsymbol\mu_1+\boldsymbol\mu_2\) i znika **dokładnie** przy antyrównoległości; dla pozytonium \(\mathbf m\propto\mathbf S_1-\mathbf S_2\), czyli waga kanału jest dokładnie operatorem mieszania singlet-tryplet |
| N3 | parzystość wymienna pola wzajemnego (sekcja 27) | człon motionalny nieparzysty, dipolowy parzysty dla para i nieparzysty dla orto — stąd \(\omega_1=\omega_2\) ściśle w orto; wyprowadzone, nie wpisane |
| N4 | libracja \(|S|\) para wracająca do **dokładnie** zera (sekcja 27b) | klasyczny cień zachowania \(S\); nikt tego nie zakładał |
| N5 | skalowanie \(\alpha^{-5}\lambda_C/c\) czasu kolapsu | ta sama potęga co \(\Gamma_{2\gamma}\), a nie została wpisana |
| N6 | klasyczna energia dipol-dipol jako ułamek rozszczepu nadsubtelnego, z resztą poprawnie zidentyfikowaną jako człon kontaktowy Fermiego i anihilacja wirtualna | liczba i jej deficyt są sprawdzalne wobec 203,3941 GHz |

### 0b. Tautologie — nie wolno ich cytować jako wyników

| # | rzecz | dlaczego tautologia |
|---|---|---|
| T1 | \(n=1\) ma promień i energię Bohra | \(a_{\rm pair}\equiv\hbar^2/(\mu k)\), więc \(-k/(2a_{\rm pair})\) **jest** algebraicznie \(-\mu k^2/(2\hbar^2)\) (sekcja 28a) |
| T2 | \(L=\hbar\) przy przygotowaniu | orbita kołowa przy \(n^2a_{\rm Ps}\) ma \(L=n\hbar\) z definicji |
| T3 | orto wychodzi poprawnie 12/12 | ścisła symetria konfiguracji, nie test, który model mógł oblać (sekcja 27a) |
| T4 | energia terminalna pod podłogą stanu podstawowego | to wartość, którą podłoga dostała |

---

## 1. Co momenty magnetyczne naprawdę dostarczają — zmierzone

Intuicja autora jest **w połowie trafna i w połowie trafiona w złą skalę**.

**Momenty niosą własną długość, i jest nią skala Comptona.**

```
mu/(e c)              = 1.933035388e-13 m
comptonBarrierRadius  = 1.933035388e-13 m     iloraz = 1.000000000000000
mu/(ec) / a_pair      = 1.826453741e-03  =  alpha * (g/2) / 4
```

**Ale nie ma rezonansu, więc nie ma czym kwantować.** Skan stosunku częstości
precesji spinu do częstości obiegu, para, momenty wzdłuż normalnej orbity:

| \(a/a_{\rm pair}\) | \(\beta\) | \(\omega_s/\omega_o\) | \(r/\lambda_C\) |
|---|---|---|---|
| \(1{,}0\cdot10^{-1}\) | 0,0231 | \(4{,}11\cdot10^{-4}\) | 27,41 |
| \(1{,}0\cdot10^{-2}\) | 0,0730 | \(8{,}68\cdot10^{-3}\) | 2,74 |
| \(3{,}16\cdot10^{-3}\) | 0,1298 | \(4{,}16\cdot10^{-2}\) | 0,867 |
| \(1{,}0\cdot10^{-3}\) | 0,2308 | \(2{,}93\cdot10^{-2}\) | 0,274 |
| \(1{,}0\cdot10^{-4}\) | 0,7297 | \(1{,}79\cdot10^{-3}\) | 0,027 |

Maksimum: \(\omega_s/\omega_o = 0{,}0440\) przy \(r=2{,}802\cdot10^{-13}\) m
\(=0{,}726\,\lambda_C\). **Do rezonansu brakuje 22,7 razy**, a krzywa zawraca,
więc żadne zamknięcie częstości nie nastąpi nigdzie.

**Wniosek.** Konfiguracja momentów wyznacza skalę Comptona i tylko ją. Do
skali Bohra brakuje **dokładnie jednej potęgi \(\alpha\)** — tej samej, której
brakowało w pracy nad tempem anihilacji (sekcja 24, „poszukiwanie brakującego
wierzchołka"). To jest ta sama luka widziana z drugiej strony, i to jest
najmocniejszy argument, że jest prawdziwa.

**Dlaczego poziomy nie wychodzą dziś.** Reguła emisji brzmi
\(E_\gamma=\hbar\omega_{\rm orb}\) i **nie zawiera momentów w ogóle**. Momenty
nigdy nie wchodzą do księgowania energii, więc nie mają jak niczego wyznaczyć.
Zmierzony ciąg lądowań kaskady (sekcja 29b) to
\(2\to1{,}414\to0{,}910\to0{,}509\to0{,}228\) — bez śladu preferencji dla
liczb całkowitych.

---

## 2. Kolejność uwalniania

Uporządkowana kosztem i tym, co każdy krok naraża.

### Krok 1 — darmowy: bariera przestaje być importem

`comptonBarrierRadius` jest dziś wpisana jako osobna stała z uzasadnieniem
„gdzie kończy się elektrodynamika punktowa". Skoro równa się \(\mu/(ec)\)
co do 15 cyfr, można ją **wyprowadzić** z konfiguracji momentów zamiast
wpisywać. Liczba się nie zmienia, więc żaden wynik nie jest zagrożony, a lista
importów skraca się o jeden. To jedyny krok bez ceny.

### Krok 2 — tani: promień startowy jako parametr ciągły

Dziś `--level n` pozwala tylko na \(n^2a_{\rm pair}\). Wprowadzić
`--start-radius` w metrach i puścić wachlarz promieni **niebędących**
promieniami Bohra. Pytanie: czy cokolwiek wyróżnia całkowite \(n\)?

Przewidywanie na podstawie sekcji 29b: **nie**. Ten krok najpewniej
potwierdzi brak emergentnych poziomów — i to jest wynik, nie porażka, bo
dziś nie jest zmierzony, tylko wnioskowany z jednego startu.

### Krok 3 — właściwy test intuicji: reguła emisji z momentami

Zamienić \(E_\gamma=\hbar\omega_{\rm orb}\) na regułę, w której momenty
występują. Kandydat nieuznaniowy: **foton na obrót precesji spinu**, zamiast
na obieg orbity. Wtedy skala emisji jest \(\omega_s\), a nie \(\omega_o\).

Zmierzona krzywa z części 1 przewiduje, co to zrobi: przy szczycie foton
padałby raz na \(\approx23\) obiegi, a poza szczytem znacznie rzadziej, więc
kaskada dostałaby **wyróżniony promień przy \(0{,}726\,\lambda_C\)**. To jest
falsyfikowalne i tanie.

Ostrzeżenie: ta reguła wyróżni skalę Comptona, nie Bohra. Jeśli intuicja ma
dać poziomy Bohra, sam krok 3 nie wystarczy.

### Krok 4 — jedyny kandydat na drugą skalę: bilans ZPF

Pole zerowe jest jedyną rzeczą w modelu, która niesie skalę niezależną od
\(\mu\), \(e\), \(m\), \(c\). Droga stochastycznej elektrodynamiki jest
dokładnie taka: promień, przy którym **moc wpompowana przez ZPF równa się
mocy wypromieniowanej**, jest samodzielnie wybraną skalą, a nie wpisaną.

Wstępna wskazówka z liczb, które model już drukuje (\(E_{\rm rms}\) ZPF rośnie
jak \(\omega^2\sim r^{-3}\), pole wiązania jak \(r^{-2}\), stosunek
\(9{,}93\cdot10^{-4}\) przy \(a_{\rm pair}\)): zrównanie **pól** wypada przy
\(1{,}05\cdot10^{-13}\) m \(=0{,}272\,\lambda_C\). Znowu skala Comptona.

**To nie jest pomiar bilansu.** Bilans SED dotyczy MOCY, nie pól, i ma inne
wykładniki. Trzeba go zmierzyć osobno — i to jest najważniejszy pojedynczy
pomiar w tym planie, bo tylko on może dać skalę Bohra.

### Krok 5 — WYKONANY, i obalił własne sformułowanie

Plan mówił: skali Bohra brakuje \(\alpha\), a \(\alpha\) w QED pochodzi
z wierzchołka. **Pierwsza połowa jest fałszywa, druga pyta o złą rzecz.**

*\(\alpha\) nie brakuje.* Jest w modelu jako samo sprzężenie,
\(k=\alpha\hbar c\) z dokładnością ULP, i każda długość to
\(\lambda_C\) razy jej potęga:

| długość | w jednostkach \(\lambda_C\) | zgodność |
|---|---|---|
| \(a_{\rm pair}\) | \(2/\alpha=274{,}0720\) | 3,7e−16 |
| \(r^*\) | \((g/2)/2=0{,}5006\) | dokładna |
| \(r_e\) | \(\alpha=7{,}297\cdot10^{-3}\) | dokładna |

`fineStructureConstant` i `reducedComptonWavelength` są teraz stałymi
pochodnymi, a `particle_species.hpp` przypina wszystkie cztery relacje
asercjami. Model ma **jedno** wejście wymiarowe, \(\lambda_C\), które
konfiguracja momentów już niesie, i **jedno** bezwymiarowe, \(\alpha\).
Obie są zwykłymi stałymi mierzonymi, a \(\alpha\) jest wejściem także
w QED.

*Brakuje powodu, żeby para SIEDZIAŁA przy \(a_{\rm pair}\).* Klasycznie,
przy ustalonym \(L\), \(E(r)=L^2/(2\mu r^2)-k/r\) ma minimum przy
\(L^2/(\mu k)\), a dla \(L=\hbar\) to jest \(a_{\rm pair}\) co do
\(1{,}1\cdot10^{-13}\). Skala bierze się więc **wyłącznie z założenia
\(L=\hbar\)**, a nie z \(\alpha\). To argument z modelu Bohra, który
mechanika kwantowa uchyliła: stan 1S ma \(L=0\), więc skala nie może
pochodzić z kwantowania momentu pędu. Kwantowo to samo minimum daje
\(\hbar^2/(2\mu r^2)\), czyli energia **lokalizacji**.

*Uczciwe przeformułowanie.* Treść skalowa modelu to nie sześć importów i nie
jedna \(\alpha\), tylko **jedno założenie \(L=\hbar\)**, zastępujące
energię lokalizacji, której model nie ma. Pozostałe pięć importów to reguły
o dyskretności, a reguły nie ustalają skali. Dlatego kroki 1–3 mogły zdjąć
barierę z listy importów, uwolnić promień startowy i wstawić momenty do
reguły emisji, nie powodując pojawienia się żadnego poziomu: żaden z nich nie
dotknął tego jednego założenia.

*I dlatego krok 4 jest jedyną pozostałą drogą.* Pole zerowe dostarcza
dokładnie tego, co zastępuje \(L=\hbar\): fluktuacji opierającej się
zapadnięciu, bez ustalania momentu pędu. To jest mechanizm, nie skala —
skala jest już dobra.

---

## 3. Jak chronić to, co cenne

Każdy z wyników N1-N6 dostaje **test regresyjny** w `positronium_validation`,
zanim ruszy krok 2. Dziś żaden z nich nie jest tam zamknięty, więc każda
zmiana reguły emisji może je po cichu zepsuć.

Priorytet: N2 i N3, bo są strukturalne i najłatwiej je złamać zmianą
transportu spinu; N1 jako jednolinijkowa tożsamość; N5 jako test skalowania.

---

## 4. Czego ten plan NIE obiecuje

- Nie da tempa anihilacji. Brak kanału kontaktowego jest osobną luką
  (sekcje 24-26) i krok 4 jej nie dotyka.
- Nie usunie \(\hbar\) z modelu. \(\hbar\) siedzi w \(\mu\), więc wchodzi
  razem z momentem; „uwolnienie od kwantu" nie jest tu na stole, na stole
  jest uwolnienie od **narzuconego poziomu**.
- Nie uratuje \(n=1\) jako wyniku. \(n=1\) jest i zostanie warunkiem
  początkowym, dopóki krok 4 nie da skali. Można jedynie przestać go
  cytować jako odtworzenie.
