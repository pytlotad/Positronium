# Klasyczny model układu elektron–pozyton

## Status i zakres modelu

Program jest numeryczną wizualizacją **klasycznej, fenomenologicznej** dynamiki
elektronu i pozytonu w próżni. Łączy elektrostatykę, magnetyzm poruszających się
ładunków, klasyczne oddziaływanie dwóch dipoli magnetycznych, relatywistyczną
kinematykę oraz przybliżony opis emisji i reakcji promieniowania.

Nazwy `Para-positronium` i `Ortho-positronium` oznaczają w programie dwie klasy
związanych trajektorii rozróżniane orientacją klasycznych momentów dipolowych.
Nie są one klasycznym wyprowadzeniem rzeczywistych stanów para- i
ortopozytonium. Spin, statystyka stanów i anihilacja pozytonium są zjawiskami
kwantowymi. Dlatego wyniki programu mogą służyć do badania konsekwencji
przyjętego modelu klasycznego, ale nie do ilościowego przewidywania własności
fizycznego pozytonium.

Wszystkie obliczenia wykonywane są w jednostkach SI. Elektron i pozyton mają tę
samą masę \(m_e\), przeciwne ładunki \(-e\) i \(+e\) oraz momenty magnetyczne o
stałej wartości magnetonu Bohra \(\mu_B\). Początkowa odległość wynosi jeden
promień Bohra \(a_0\), a środek masy układu spoczywa.

## Uwzględniane efekty fizyczne

### 1. Relatywistyczna kinematyka cząstek

Stan dynamiczny jest przechowywany za pomocą położenia i pędu. Zależność
między pędem i prędkością ma postać

\[
\mathbf p=\gamma m\mathbf v,\qquad
\gamma=\frac{1}{\sqrt{1-v^2/c^2}}.
\]

Energia kinetyczna każdej cząstki jest liczona jako
\(K=(\gamma-1)mc^2\). Ograniczenie \(v<c\) wynika więc z użytej relacji
pęd–prędkość, a nie z dodatkowego obcinania prędkości.

### 2. Przyciąganie Coulomba

Podstawową siłą wiążącą jest chwilowa siła Coulomba

\[
\mathbf F_C=\frac{1}{4\pi\varepsilon_0}
\frac{q_1q_2}{r^3}\mathbf r,
\qquad
U_C=\frac{1}{4\pi\varepsilon_0}\frac{q_1q_2}{r},
\]

gdzie \(\mathbf r=\mathbf r_1-\mathbf r_2\). Dla przeciwnych znaków ładunku
jest to siła przyciągająca. Program nie dodaje sztucznej siły odpychającej,
która miałaby stabilizować orbitę; niezerowy orbitalny moment pędu może jednak
powodować klasyczną barierę odśrodkową wynikającą z samego ruchu.

### 3. Magnetyzm orbitalny i magnetyczna część siły Lorentza

Poruszający się ładunek jest źródłem pola w przybliżeniu niskich prędkości

\[
\mathbf B=\frac{\mu_0}{4\pi}
q\frac{\mathbf v\times\mathbf R}{R^3}.
\]

Druga cząstka podlega magnetycznej części siły Lorentza
\(\mathbf F_B=q\,\mathbf v\times\mathbf B\). Jej działanie jest realizowane
jako obrót pędu (krok typu Boris/Rodrigues), dzięki czemu sama część
magnetyczna nie zmienia modułu pędu i nie wykonuje pracy. Jest to poprawna
własność siły \(q\mathbf v\times\mathbf B\).

Jest to model chwilowy i niskoprędkościowy, a nie pełne rozwiązanie równań
Maxwella z czasem opóźnionym. Nie stanowi też kompletnego przybliżenia Darwina,
ponieważ nie uwzględnia wszystkich poprawek zależnych od prędkości tego samego
rzędu.

### 4. Oddziaływanie klasycznych dipoli magnetycznych

Każdej cząstce przypisany jest punktowy moment magnetyczny \(\boldsymbol\mu\)
o wartości \(\mu_B\). Dla dwóch dipoli używana jest standardowa energia

\[
U_{dd}=\frac{\mu_0}{4\pi r^3}
\left[\boldsymbol\mu_1\!\cdot\!\boldsymbol\mu_2
-3(\boldsymbol\mu_1\!\cdot\!\hat{\mathbf r})
(\boldsymbol\mu_2\!\cdot\!\hat{\mathbf r})\right],
\]

a siła dipol–dipol jest wyznaczana z gradientu tej energii. Orientacje dipoli
są losowane na początku i pozostają stałe podczas symulacji. Program nie
całkuje precesji ani momentów sił obracających dipole.

Punktowe wzory magnetyczne są osobliwe dla \(r\to0\), dlatego pole orbitalne,
energia i siła dipolowa są płynnie tłumione funkcją

\[
w(r)=\frac{1}{1+(0{,}70a_0/r)^6}.
\]

W sile dipolowej uwzględniany jest również gradient \(w(r)\), dzięki czemu
siła pozostaje zgodna z regularizowaną energią potencjalną \(wU_{dd}\).
Regularizacja jest założeniem numeryczno-modelowym, a nie nowym prawem
fizycznym.

### 5. Koherentne promieniowanie elektromagnetyczne pary

Elektron i pozyton są traktowane jako jedno, koherentne źródło promieniowania
w przybliżeniu dipola elektrycznego. Moment dipolowy całego układu i jego druga
pochodna wynoszą

\[
\mathbf p=\sum_i q_i\mathbf r_i=e(\mathbf r_p-\mathbf r_e),
\qquad
\ddot{\mathbf p}=e(\mathbf a_p-\mathbf a_e).
\]

Moc promieniowania jest liczona z kwadratu sumy amplitud pól obu ładunków:

\[
P_{dip}=\frac{|\ddot{\mathbf p}|^2}{6\pi\varepsilon_0c^3}.
\]

Zawiera to automatycznie człony interferencyjne, których brakowałoby w sumie
dwóch niezależnych mocy Larmora lub Liénarda. Moc jest całkowana numerycznie w
czasie i prezentowana jako skumulowana energia promieniowania. Przybliżenie
wymaga, aby rozmiar układu był mały względem charakterystycznej długości fali;
nie obejmuje jeszcze promieniowania dipola magnetycznego ani kwadrupola
elektrycznego.

### 6. Kolektywna reakcja promieniowania

Wpływ emisji na ruch jest wprowadzony bezpośrednio do równania dynamiki przez
pole reakcji odpowiadające temu samemu przybliżeniu dipolowemu:

\[
\mathbf E_{RR}=\frac{\dddot{\mathbf p}}
{6\pi\varepsilon_0c^3},
\qquad
\mathbf F_{RR,i}=q_i\mathbf E_{RR}.
\]

Trzecia pochodna momentu dipolowego jest obliczana przez centralną różnicę
przyspieszeń wyznaczonych z siły Coulomba, regularizowanej siły dipol–dipol i
orbitalnej siły magnetycznej. Jest to redukcja rzędu: nie wprowadza niezależnej
zmiennej przyspieszenia i nie dopuszcza rozwiązań samoprzyspieszających.
Reakcja wpływa na kolejne położenia i pędy; energia promieniowania nie służy do
późniejszego korygowania trajektorii. Zastosowany wzór pozostaje
nierelatywistycznym przybliżeniem długofalowym, a nie pełnym kowariantnym
równaniem reakcji promieniowania.

### 7. Energia bliskiego pola Schotta

Do diagnostycznego bilansu energii dołączana jest energia Schotta

\[
E_S=-\frac{\ddot{\mathbf p}\cdot\dot{\mathbf p}}
{6\pi\varepsilon_0c^3},
\]

dla wspólnego momentu dipolowego pary. Reprezentuje ona odwracalną wymianę
energii z bliskim polem układu. Raportowana wielkość kontrolna to

\[
E_{diag}=K_e+K_p+U_C+wU_{dd}+E_{rad}+E_S.
\]

Bilans ten jest diagnostyką konsekwencji równań ruchu, a nie narzuconą zasadą.
Nie obejmuje energii pola magnetycznego wytworzonego przez ruch ładunków, więc
nie należy oczekiwać jego dokładnej stałości. Pozostaje też błąd dyskretyzacji
numerycznej.

## Warunki początkowe i klasyfikacja zjawiska

Program losuje kierunki dipoli oraz radialną i styczną składową względnej
prędkości. Wybrana pozycja menu ogranicza zakres losowania tak, aby wszystkie
cztery rodzaje trajektorii były osiągalne. Następnie wynik określa jedno
klasyczne drzewo decyzyjne, używające:

- energii ruchu względnego \(E_{rel}=\tfrac12\mu v_{rel}^2-ke^2/r\), gdzie
  \(\mu=m_e/2\) jest masą zredukowaną;
- orbitalnego momentu pędu \(L_{orb}=\mu|\mathbf r\times\mathbf v_{rel}|\);
- przewidywanej odległości największego zbliżenia \(r_{min}\) z problemu
  Keplera–Coulomba;
- współczynnika ustawienia dipoli
  \(A=(\boldsymbol\mu_e\cdot\boldsymbol\mu_p)/\mu_B^2\).

Klasyfikacja przebiega następująco:

1. Ruch do wewnątrz i \(r_{min}<10^{-14}\,\mathrm m\) daje
   `Direct collision`.
2. W przeciwnym przypadku \(E_{rel}\ge 0\) daje `Scattering`.
3. Dla \(E_{rel}<0\), gdy \(A\ge0{,}5\), wynikiem jest
   `Para-positronium`.
4. Pozostały związany przypadek daje `Ortho-positronium`.

Próg orientacji dipoli jest fenomenologiczną regułą programu. Nie jest
klasycznym odpowiednikiem kwantowego singletu i trypletu.

`Lifetime` oznacza faktyczny czas symulacji do osiągnięcia umownego progu
anihilacji \(r\le10^{-14}\,\mathrm m\). Jeśli próg nie zostanie osiągnięty w
oknie czasowym, wyświetlana jest nieskończoność; dla rozproszenia oczekiwanym
wynikiem jest `∞`. Program nie modeluje samego procesu anihilacji ani produktów
końcowych.

## Metoda numeryczna

Równania są całkowane relatywistycznym schematem predyktor–korektor. Siły
zależne od położenia wykonują dwa półkroki pędu, część magnetyczna obraca pęd,
a pomiędzy nimi wykonywany jest krok położenia. Krok czasu jest adaptacyjny:
jest ograniczony do \(2\times10^{-18}\,\mathrm s\) i dodatkowo zmniejszany na
podstawie chwilowej częstości Coulombowskiej. Symulacja zapisuje najwyżej 1200
klatek i kończy się wcześniej po przekroczeniu progu zderzenia.

Równe masy i przeciwne siły wewnętrzne pozwalają utrzymywać środek masy w
spoczynku z dokładnością numeryczną. Obie cząstki są rysowane jako kulki tej
samej wielkości.

## Wielkości prezentowane na ekranie

Widok ROOT pokazuje trajektorie i wektory momentów magnetycznych. Strzałki
`Spin ↑/↓` w tabeli są jedynie skrótem graficznym: wskazują znak składowej
\(z\) klasycznego momentu magnetycznego, a nie pomiar kwantowego spinu.

Nad tabelą widoczne są: wynik klasyfikacji, \(E_{rel}\), \(L_{orb}\),
\(r_{min}\) i `Lifetime`. Tabela bieżąca zawiera czas, odległość, energię
mechaniczną oraz skumulowaną energię promieniowania. Wykres położenia używa
jednostek \(a_0\), odległość jest podawana w pm, czas w ps, a energie w eV.

## Uruchomienie i odtwarzalność

Wymagany jest CERN ROOT z programem `root-config` dostępnym w `PATH`.

```bash
make build
./positronium
```

Samo `make` kompiluje program i od razu go uruchamia.

Menu początkowe:

```text
1 -> Para-positronium
2 -> Ortho-positronium
3 -> Direct collision
4 -> Scattering
```

Wybór i ziarno generatora można podać bez interakcji:

```bash
./positronium --phenomenon 2 --seed 42
```

Tryb diagnostyczny nie otwiera okna i wypisuje zakres odległości oraz bilans
energii:

```bash
./positronium --diagnose --phenomenon 4 --seed 42
```

Przyciski `STOP`/`START` sterują animacją, a `EXIT` zamyka program.

## Ograniczenia — efekty nieuwzględniane

- mechanika kwantowa, funkcja falowa, zasada Pauliego, splątanie i kwantowa
  struktura stanów para-/ortopozytonium;
- kwantowa anihilacja oraz emisja dwóch albo trzech fotonów;
- pełne, wzajemnie opóźnione pola Liénarda–Wiecherta i kompletna dynamika pola
  Maxwella;
- kompletne relatywistyczne oddziaływanie dwóch ładunków, w tym pełny człon
  Darwina;
- dynamiczna precesja dipoli, sprzężenie spin–orbita i momenty sił działające
  na dipole;
- skończony rozmiar cząstek i struktura krótkiego zasięgu — zastępuje je próg
  \(10^{-14}\,\mathrm m\);
- pola zewnętrzne, zderzenia z materią, grawitacja i wpływ ośrodka;
- energia orbitalnego pola magnetycznego w prezentowanym bilansie energii;
- wyższe multipole promieniowania: dipol magnetyczny, kwadrupol elektryczny i
  kolejne człony rozwinięcia.

## Literatura pomocnicza

1. NIST, *2022 CODATA Recommended Values of the Fundamental Physical
   Constants*, [JPCRD 53, 030201 (2024)](https://physics.nist.gov/cuu/pdf/JPCRD2022CODATA.pdf).
2. H. Spohn, *The critical manifold of the Lorentz–Dirac equation*,
   [Europhys. Lett. 50, 287 (2000)](https://doi.org/10.1209/epl/i2000-00233-3).
3. C. Bild, D.-A. Deckert i H. Ruhl, *Radiation reaction in classical
   electrodynamics*, [Phys. Rev. D 99, 096001 (2019)](https://doi.org/10.1103/PhysRevD.99.096001).
4. C. G. Darwin, *The Dynamical Motions of Charged Particles*,
   [Philosophical Magazine 39, 537–551 (1920)](https://doi.org/10.1080/14786440508636066).

Literatura uzasadnia użyte elementy elektrodynamiki klasycznej, ale nie
waliduje fenomenologicznego utożsamienia orientacji dipoli ze stanami
para-/ortopozytonium.
