# Klasyczny model układu elektron–pozyton

## Użycie i wskazanie autorstwa modelu CREM

Program może być używany przez każdego. Publikacja wyników uzyskanych za jego
pomocą — również wyników opartych na zmodyfikowanym kodzie albo zmodyfikowanym
modelu CREM — powinna wskazywać położenie katalogu projektu Positronium:

```text
https://github.com/pytlotad/Positronium
```

oraz autora modelu CREM: **Tadeusz Sławomir Pytlos**, dostępnego pod adresem
`tadeusz.slawomir.pytlos@gmail.com`.

## Status i zakres modelu

Program jest numeryczną wizualizacją **klasycznej, fenomenologicznej** dynamiki
elektronu i pozytonu w próżni. Łączy elektrostatykę, magnetyzm poruszających się
ładunków, klasyczne oddziaływanie dwóch dipoli magnetycznych, relatywistyczną
kinematykę oraz przybliżony opis emisji i reakcji promieniowania.

W trybie wizualnym nazwy `Para-positronium` i `Ortho-positronium` oznaczają
dwie klasy związanych trajektorii rozróżniane orientacją klasycznych momentów
dipolowych. Nie są one klasycznym wyprowadzeniem rzeczywistych stanów para- i
ortopozytonium. Spin, statystyka stanów i anihilacja są zjawiskami kwantowymi.
Dlatego wyniki dynamiki orbitalnej służą do badania konsekwencji przyjętego
modelu klasycznego, a nie do ilościowego przewidywania własności fizycznego
pozytonium. Tryb statystyczny stanów związanych korzysta z osobnego, jawnie
oznaczonego generatora zaniku w idealnej próżni; nie udaje, że anihilacja
wynika z dojścia klasycznej trajektorii do małej odległości.

### Promień regularyzacji dipola magnetycznego

Pole dipola magnetycznego jest wygładzane wagą \(w(r)=1/(1+(a/r)^p)\) z
\(p=6\) i

\[
a=\sqrt[3]{
  \frac{\mu_0}{4\pi}\frac{\mu_1\mu_2}{E/2}\,K_6
}=68{,}47\,\mathrm{fm},
\qquad K_6=1{,}5244264856,
\]

gdzie \(E\) jest energią spoczynkową **lżejszego** składnika pary. Wchodzą oba momenty,
\(\mu_1\mu_2\), a nie kwadrat jednego: to ta sama liczba wyłącznie dla cząstki
i jej antycząstki, bo moment protonu jest 658 razy mniejszy od elektronowego,
więc p+e⁻ potrzebuje 7,87 fm tam, gdzie e⁺e⁻ potrzebuje 68,47. Promień jest
wyprowadzany z wybranej pary, a nie wpisany jako stała.

Momenty są przy tym tymi, które niesie dynamika, czyli \((g/2)\mu_B\) dla
elektronu, a nie gołym magnetonem Bohra. Wzór z \(\mu_B^2\) i masą elektronu
pomijałby zależność od wybranej pary. Współczynnik \(K_6\) nie pochodzi z
samego \(w/r^3\), lecz z rzeczywistego pola
\(\mathbf B_{\rm reg}=\boldsymbol\nabla\times\mathbf A_{\rm reg}\). Największa
wartość poprzeczna wynosi
\(x^3(5-x^6)/(1+x^6)^2=K_6\) przy
\(x=r/a=(9-2\sqrt{19})^{1/6}\). Dzięki temu klasyczna energia
\(-\boldsymbol\mu_1\cdot\mathbf B_{2,\rm reg}\), a nie jej uproszczony
substytut, jest ograniczona do \(E/2\) w całej dziedzinie (dla e⁺e⁻:
\(m_ec^2/2\)).
Limit dotyczy tablicowych wartości momentów w ich wspólnym układzie
spoczynkowym; nie jest obietnicą niezmienniczego ograniczenia laboratoryjnej
energii po dowolnie dużym boostcie.

Poprzednia wartość \(0{,}5\cdot\)`nuclearCutoff` \(=5\,\mathrm{fm}\) była
dobrana wyłącznie tak, aby wygładzanie pozostało wewnątrz raportowanej
dziedziny. Uproszczony skalar \(w/r^3\) wskazywał
\(2{,}1\cdot10^8\,\mathrm{eV}\), lecz pełne `curl(A)` dopuszczało około
\(6{,}6\cdot10^8\,\mathrm{eV}\), czyli 1284 masy spoczynkowe elektronu.
Wartość 5 fm czyniła też sektor magnetyczny 106 razy
bardziej punktowym niż sektor ładunkowy, który już zakłada źródło o skończonym
promieniu `chargeCloudRestRadius`. Promień 68 fm nadal leży daleko poniżej
skali atomowej; przy dziesiątkach pikometrów korekta wagi jest mniejsza niż
kilka części na \(10^{12}\).

Dla e⁺e⁻ nie usuwa to krótkozasięgowej bariery dipolowej. Bariera leży tam, gdzie
nieregularyzowane energie dipolowa i kulombowska się zrównują,
\(r^\*=\sqrt{(\mu_0/4\pi)\mu_1\mu_2/(ke^2)}=193\,\mathrm{fm}\), czyli w
połowie zredukowanej długości Comptona i poza promieniem wygładzania.
Poniżej tej skali klasyczny
opis punktowego momentu magnetycznego nie ma podstaw fizycznych.

Wszystkie obliczenia wykonywane są w jednostkach SI. Elektron i pozyton mają tę
samą masę \(m_e\), przeciwne ładunki \(-e\) i \(+e\) oraz momenty magnetyczne o
stałej wartości \((g/2)\mu_B\). Początkowa odległość wynosi promień Bohra
**pozytonium** \(a_{Ps}=2a_0=1{,}058\,\mathrm{\AA}\), a środek masy układu
spoczywa.

\(a_0\) jest promieniem Bohra wodoru, zbudowanym z masy elektronu; pozytonium
ma masę zredukowaną \(\mu=m_e/2\), więc jego skala to \(\hbar^2/(\mu ke^2)=2a_0\).
Wcześniej para była przygotowywana w \(a_0\) z podkołową prędkością styczną, co
dawało orbity o \(L=0{,}51\!-\!0{,}69\,\hbar\), półosi \(\approx0{,}8a_0\) i
energii wiązania \(\approx17\,\mathrm{eV}\) — czyli **wewnątrz** stanu
podstawowego, który model ma reprezentować. Ponieważ klasyczny czas inspirali
skaluje się jak \(a^3\), samo to skracało kolaps ponad dwudziestokrotnie.
Pasmo prędkości stycznej jest teraz wyśrodkowane na orbicie kołowej w
\(a_{Ps}\), czyli na stanie \(L=\hbar\).

\(L=\hbar\) to wartość Bohra/SED, nie kwantowa wartość stanu 1s, która wynosi
\(L=0\). Klasyczna orbita o \(L=0\) jest radialnym spadkiem, więc nie może
zastąpić stanu podstawowego; \(L=\hbar\) jest orbitą klasyczną odtwarzającą
poprawną energię wiązania i średni promień.

### Granica zderzenia

Separacja, przy której trajektoria liczy się jako zderzenie, też skaluje się
z parą:

\[
r_{\rm coll}=0{,}005\,a_{\rm pary},
\]

co daje 0,529177 pm dla e⁺e⁻, 0,00255928 pm dla µ⁺µ⁻ i 0,000288199 pm dla
p+p̄. Wcześniej była to stała 0,01·a₀ = 529 fm, zbudowana z promienia Bohra
**wodoru**, podczas gdy separacja startowa od dawna brała promień Bohra
**wybranej pary**. Dla µ⁺µ⁻ start wypada na 512 fm, a dla p+p̄ na 57,6 fm,
czyli już wewnątrz tej granicy: pętla uznawała trajektorię za zderzoną przed
pierwszym krokiem i zwracała awarię numeryczną. Po przeskalowaniu eksperyment 5
dla mionium schodzi z 2/12 zderzeń do 0/12. Granica jest wyprowadzana
w `particle_species.hpp` jako `collisionBoundaryOf(pair)`, obok promienia
regularyzacji — obie opuściły `physical_constants.hpp` z tego samego powodu:
niosły skalę gatunku, a nie czystą stałą fizyczną.

### Siatka historii retardowanej

Węzeł historii nie jest już zapisywany na każdy przyjęty krok. Odstęp węzłów
jest przypięty do okna retardowanego,

\[
\Delta t_{\rm hist}=\frac{\max(10^{-20}\,\mathrm{s},\;4r/c)}{128},
\]

a nowszy stan, który nie sięga tego odstępu, tylko odświeża czoło zamiast
zagęszczać siatkę.

Powód jest mierzalny. Pochodne retardowane są różnicowane po \(h=2\Delta
t_{\rm hist}\), a sektor Schotta potrzebuje **trzeciej** pochodnej, której
stencil dzieli przez \(h^3\). Gdy odstęp szedł za krokiem, połowienie kroku
ośmiokrotnie wzmacniało zaokrąglenia w tym członie. Zagęszczanie kroku
psuło więc kanał prędkości: na trajektorii kolapsu reszta **pozycji** opadała
jak \(dt^2\) z \(5{,}0\cdot10^{-9}\) do \(4{,}6\cdot10^{-13}\),
podczas gdy reszta **prędkości** w tych samych trzech połowieniach *rosła*
z \(1{,}25\cdot10^{-7}\) do \(4{,}31\cdot10^{-7}\). Silnik adaptacyjny
wyczerpywał wtedy głębokość podziału i zgłaszał awarię numeryczną przy 1,1–2,2
promienia zderzenia — dlatego tolerancji nie dawało się zacieśnić poniżej
\(10^{-6}\) w ogóle.

Po zmianie tolerancja schodzi niżej, a czas kolapsu zbiega monotonicznie:
3,252639e-05, 3,2526313e-05, 3,252627e-05 i 3,252626e-05 ps dla tolerancji
\(10^{-5}\) do \(10^{-8}\), gdzie wcześniej dwie ostatnie kończyły się
awarią. Sektor koherentnej reakcji wyraźnie zyskał: stosunek strumienia
przeszedł z \(1{,}1139\) na \(0{,}9965\) przy ideale 1, a reszta kroku
z \(8{,}40\cdot10^{-6}\) na \(1{,}70\cdot10^{-8}\), czyli 493 razy.
Koszt jest ujemny — walidacja skróciła się z 19,77 s do 16,56 s, bo węzłów
jest mniej. Pogorszyła się jedna pozycja: skumulowana reszta czteropędu po
boostcie z \(2{,}08\cdot10^{-6}\) na \(2{,}33\cdot10^{-6}\).

Gęstość 128 węzłów nie jest dowolna. Przy 64 reszta kroku przestawała zbiegać
z tolerancją (9,2/6,8/6,2\(\cdot10^{-6}\) zamiast 5,4/1,3/0,30), a przy 256
znikał zysk w sektorze koherentnym — siatka musi być dość gęsta, by rozdzielić
okno, i dość rzadka, by odczepić się od kroku. Limit `maximumHistoryNodes`
jest teraz jej dwukrotnością i w zwykłym biegu nie jest osiągany.

Cztery inne wyjaśnienia zostały po drodze **obalone pomiarem**: limit 128
węzłów (wynik bit-w-bit identyczny), rozejście się historii między próbą
zgrubną a dokładną, podniesienie `maximumDepth` do 16 (445-krotny czas i dalej
awaria) oraz sztywna podłoga \(h\), która przy drobnym kroku dała
\(2{,}3\cdot10^{8}\) eV energii wypromieniowanej.

Zmiana kupuje **odporność, nie domknięty bilans radiacyjny**. Po usunięciu
ściany zacieśnienie tolerancji do \(10^{-7}\) wciąż przesuwa energię
wypromieniowaną z 0,412 na 0,345 eV, a niedopasowanie reakcji zostaje przy
2,4 eV wobec 0,41 eV promieniowania.

Warto zapisać, **czym ta reszta nie jest**, bo nasuwa się to samo błędne
wyjaśnienie. Nie jest załamaniem rozwinięcia zredukowanego: na granicy
zderzenia e⁺e⁻ parametr Landaua–Lifshitza wynosi \(\tau\omega=3{,}7\cdot
10^{-4}\) przy \(\tau=2r_e/3c=6{,}26\cdot10^{-24}\) s i
\(\omega=5{,}9\cdot10^{19}\) rad/s, a \(\beta=0{,}10\). Rozwinięcie
siły własnej jest tam znakomicie zbieżne.

Nie jest też **brakującym członem interferencyjnym**.
`individualLandauLifshitzSelfForces` zwraca `firstSelf+firstMutual`, gdzie
część wzajemna to dokładnie wyraz \(2q_1q_2\,\mathbf a_1\!\cdot\!\mathbf
a_2\); ta droga jest więc algebraicznie tożsama z reakcją dipola koherentnego
\(q_i\dddot{\mathbf d}/\mathrm{norm}\). Człon w kodzie jest.

Rozkład wzdłuż kolapsu (e⁺e⁻, ziarno 12345, w eV): strumień retardowany 0,422,
E1 koherentny 0,527, suma indywidualnych 0,268, **interferencja 0,259**, praca
reakcji **+2,071**. Kontrola instrumentacji: E1/indywidualny = 1,97 wobec
dokładnie 2 wymaganych dla pary symetrycznej. Interferencja jest więc
dziewięć razy za mała, by wyjaśnić lukę 2,43 eV, a praca wychodzi **dodatnia**,
czyli zdominowana przez odwracalny człon Schotta, nie przez tłumienie.
\(E_S\) na starcie zmierzono, a nie założono: −1,2·10⁻⁶ eV.

Strukturalnie luka idzie za \(-\Delta E_S^{\rm koh}\):
\(\int\Phi-\Delta E_S^{\rm koh}-\int P_{E1}=0{,}422+2{,}667-0{,}527=2{,}562\)
wobec zmierzonych 2,432. Człon Schotta — odwracalna energia pola bliskiego, a
nie złamanie zachowania — jest więc główną jej częścią, a właściwa część
promienista (strumień − \(P_{E1}=-0{,}105\) eV) to mniej więcej poprawka
retardacyjna przy \(\beta=0{,}1\).

Sama `individualLandauLifshitzSelfForces` jest **poprawna**. Sprawdzona
lokalnie, krok po kroku, jej tożsamość \(W=-\Delta E_S^{\rm koh}-P_{E1}dt\)
domyka się do 2,5–10%, a e⁺e⁻ i µ⁺µ⁻ zgadzają się na cztery cyfry przy równym
\(r/r_{\rm coll}\) (0,0255 wobec 0,0254) — sektor reakcji jest niezmienniczy
względem skali, jak wymaga wspólne \(\tau\omega=3{,}7\cdot10^{-4}\).

Dwa wcześniejsze odczyty tej wielkości były **artefaktami pomiaru** i są tu
zapisane, żeby nikt ich nie odkrywał ponownie. Porównywanie **całki** pracy
z **punktową** wartością \(E_S\) na końcu każe tożsamości chybiać o 75–97%
dla µ⁺µ⁻ — bo \(E_S\) zmienia się najszybciej właśnie na końcu. Liczenie
pracy z lewego końca kroku przy dokładnym \(\Delta E_S\) zostawia pozorną
resztę \(-0{,}5\,P\,dt\). Żaden z nich nie przeżywa wyrównania kwadratur.
Obniżenie podłogi `derivativeStep` \(10^{-24}\) tysiąckrotnie przesuwa wynik
mionowy o 0,6%, więc i ona jest bez związku.

Przy sprawdzaniu tego wyszła jedna pułapka pomiarowa warta zapisania: na
**zwykłej orbicie związanej** okno retardowane \(\max(10^{-20},4r/c)\) jest
**krótsze niż krok całkowania**, więc historia trzyma około trzech węzłów — za
mało na stencil trzeciej pochodnej, który wtedy poprawnie zwraca zero. Droga
koherentna jest tam z konstrukcji niedostępna, a porównywanie jej z drogą LL
oznacza porównywanie z wektorem bliskim zeru (stosunek rzędu \(10^8\),
kosinus skaczący od −0,95 do +0,77). To artefakt pomiaru, nie druga usterka.

Cztery pary dają 33/33. Wyniki produkcyjne przy tolerancji \(10^{-5}\) są
dla e⁺e⁻ bit-w-bit te same; dla µ⁺µ⁻ czas kolapsu różni się na ósmej cyfrze
(1,5733404e-07 wobec 1,5733401e-07 ps), czyli poniżej wykazanego poziomu
zbieżności.

### Kwadratura księgowości radiacyjnej

Szybkości promieniowania są próbkowane na **początku** kroku, więc ważenie
każdej próbki jej własnym \(dt\) było sumą Riemanna z lewego końca — regułą
pierwszego rzędu. Blisko kolapsu moc podwaja się mniej więcej co krok, więc
ten błąd sięga połowy samej mocy.

Zmierzone na trajektorii kolapsu wobec tożsamości
\(W=-\Delta E_S^{\rm koh}-P_{E1}dt\), co piętnasty krok, przy \(r/r_{\rm
coll}\) biegnącym 104 → 4,9:

| reguła | reszta / \(P\,dt\) |
|---|---|
| lewy koniec | −0,298, −0,435, −0,492, −0,503 |
| trapez | +0,026, +0,041, +0,052, +0,064 |

e⁺e⁻ i µ⁺µ⁻ zgadzają się przy tym na cztery cyfry, więc rzecz dotyczy
kwadratury, nie pary.

Akumulator ma teraz postać

\[
\Delta A_k=f_k\,dt_k+\tfrac12\,(f_k-f_{k-1})\,dt_{k-1},
\]

czyli każdy **zamknięty** przedział jest dokładnym trapezem, a bieżący nosi
tymczasową wartość z lewego końca, którą poprawia następny krok. Kosztuje to
jeden `bool` i piętnaście liczb w `State` — żadnej dodatkowej ewaluacji siły.
Wcześniejsza próba trapezu brała prawy koniec ze **stanu próbnego** i zaniżała
o 34%, bo tamta rekonstrukcja jest niespójna; ta postać nie sięga poza
trajektorię zatwierdzoną.

Odrzucono też wariant ważący próbkę przez \((dt_{k-1}+dt_k)/2\). Jest to ta
sama reguła przegrupowana, ale gubi ostatnie pół przedziału, co kosztowało
0,34% na kontroli `Larmor accumulation` (0,99998 → 0,99656) — deficyt jest tam
połową **wartości** mocy. Postać powyżej zostawia połowę jej **przyrostu**,
który znika, gdy moc zmienia się wolno.

Skutek jest taki, że energia wypromieniowana **zaczęła zbiegać**:

| tolerancja | lewy koniec | trapez |
|---|---|---|
| \(10^{-5}\) | 0,41196 | 0,42603 |
| \(10^{-6}\) | 0,41882 (+1,67%) | 0,42611 (**+0,019%**) |
| \(10^{-7}\) | 0,34528 (−17,6%) | 0,42820 (**+0,49%**) |

Rozrzut spadł z 17,6% do pół procenta, a wartość produkcyjna przy kolapsie
przesunęła się o +3,4%, bo lewy koniec ją zaniżał. W audycie `pair-field
IDENTITY` poprawiło się pięciokrotnie (1,03·10⁻²⁵ → 2,05·10⁻²⁶), skumulowana
reszta czteropędu po boostcie o 2,5% (2,33 → 2,27·10⁻⁶), a `Larmor
accumulation` zostało nietknięte (0,99998 → 0,99997).

Czego to **nie** zmienia: stosunek niedopasowania reakcji przy kolapsie stoi
(5,905 → 5,912). I tak być powinno — jak opisano wyżej, tworzy go człon
Schotta, a nie błąd kwadratury. Poprawka naprawia całkowanie, nie
interpretację.

Cztery pary dają 33/33.

### Podłoga zasięgu CREM: bariera Comptona zamiast granicy zderzenia

Poprzednia podsekcja usunęła sprzężenie siatki historii retardowanej z krokiem
integratora, które wcześniej powodowało awarię numeryczną przy dowolnym
żądaniu peryapsis poniżej z grubsza \(529\) fm (patrz niżej, "Klasyczna
gęstość kontaktowa..."). Ta poprawka nie została tam zmierzona pod kątem
głębokości — sprawdzono tylko zbieżność czasu kolapsu i sektora reakcji.
Zmierzono to osobno, po fakcie.

Metoda: `estimateCremCollapse` podmieniono tymczasowo tak, by cel pomiaru
mechanicznego (dotąd `collisionBoundaryRadius` = \(529\) fm, rozmiar chmury
ładunku) wskazywał \(1\) fm — celowo nieosiągalnie głęboko, żeby wymusić serię
prawdziwych pomiarów mechanicznych aż do faktycznej awarii, i zlokalizować, GDZIE
ona teraz leży. Powtórzone na sześciu ziarnach (42, 1, 7, 23, 99, 2024),
eksperyment 1, po \(6\) trajektorii na ziarno.

Wynik, \(21\) zarejestrowanych awarii numerycznych (peryapsis w chwili awarii,
w fm):

| ziarno | awarie (fm) |
|---|---|
| 42 | 164,6 / 191,8 / 133,5 / 149,3 / 211,4 / 231,1 |
| 1 | 287,0 / 155,3 / 68,6 |
| 7 | 163,9 / 46,6 / 131,1 |
| 23 | 139,5 / 188,8 |
| 99 | 62,4 / 137,5 / 108,8 / 495,1 / 133,3 |
| 2024 | 210,7 / 227,9 |

Mediana \(155\) fm, średnia \(173\) fm, odchylenie standardowe \(92\) fm,
zakres \(47\)–\(495\) fm. **Żadne** ziarno nie odtworzyło starej skali
\(529\) fm poza jednym odstającym przypadkiem (ziarno 99, \(495\) fm) — awarie
skupiają się o rząd wielkości głębiej, w paśmie obejmującym barierę Comptona
\(r^*=193{,}30\) fm z poprzedniej podsekcji: \(71\%\) awarii wypada poniżej
niej, \(43\%\) w promieniu \(\pm20\%\) od niej. To nie jest ostra ściana —
zgodnie z wcześniejszym pomiarem na starym silniku ("nie ma ściany, tylko
malejący z głębokością odsetek dochodzących") — ale środek ciężkości
niepowodzeń przesunął się z \(529\) fm na coś bliskiego \(r^*\).

Wniosek praktyczny: stara wartość `collisionBoundaryRadius` nigdy nie była
właściwym celem dla tego pomiaru — była pod ręką, bo `chargeCloudRestRadius`
(rozdzielczość przestrzenna modelu) to jedyna skala "jak blisko jest za blisko"
dostępna w chwili, gdy CREM powstawał. Fizycznie właściwym celem jest bariera
Comptona: miejsce, w którym lokalizacja elektronu wymaga pędów rzędu \(mc\) i
klasyczna elektrodynamika punktowej cząstki przestaje obowiązywać — nie
rozmiar chmury ładunku, który jest osobnym, niepowiązanym wyborem
regularyzacji.

Zmiana produkcyjna: `comptonBarrierRadius` (`physical_constants.hpp`,
\(=g\hbar/4m_ec=193{,}30\) fm) zastąpił `collisionBoundaryRadius` jako cel
pomiaru mechanicznego w `crem_collapse.hpp` — wyłącznie tam; ogólna granica
zderzenia (`collisionBoundaryRadius`, eksperyment 5) zostaje nietknięta, bo to
osobne, poprawnie działające pojęcie. `finalApproachMultiple=10` zostawiono
bez zmian: próg zatrzymania pętli sekularnej wynosi teraz \(10\times
r^*=1933\) fm — wciąż bezpiecznie powyżej całego zmierzonego pasma awarii
(\(47\)–\(495\) fm), więc pomiar mechaniczny w normalnej pracy nigdy nie
zbliża się do strefy ryzyka. Potwierdzone: \(5\) ziaren \(\times\) \(6\)
trajektorii = \(30/30\) bez żadnej awarii po zmianie, walidacja pól nadal
\(33/33\) (sektor CREM jej nie dotyczy).

Czego to nie zmienia: obcięcie ostatniego odcinka nadal jest tylko \(0{,}1\%\)
czasu kolapsu (ten sam argument \(t\sim a^3\), teraz liczony od \(r^*\) zamiast
od \(529\) fm), więc raportowany czas kolapsu praktycznie się nie zmienia.
Zmienia się natomiast GDZIE model deklaruje koniec swojej ważności — z
przypadkowej skali chmury ładunku na skalę wyprowadzoną z fizyki.

Ta sama poprawka dotyczy też Eksperymentu 5 ("Interactions"): schwytana para
tam też bywa klasyfikowana jako para- lub orto-pozytonium
(`InteractionOutcome::ParaPositronium`/`OrthoPositronium`), a "Bound-state CREM
collapse time" ekstrapoluje dokładnie ten sam sekularny model co powyżej. Cel
tej ekstrapolacji zmieniony analogicznie, warunkowo: `comptonBarrierRadius`
dla schwytanej pary e⁺e⁻, `collisionBoundaryRadius` bez zmian dla każdej innej
pary (`comptonBarrierRadius` jest stałą specyficzną dla elektronu, więc
podmiana na sztywno byłaby fizycznie błędna dla mionium czy protonium).
Granica, przy której trajektoria jest deklarowana jako dosłowne "Collision"
(zatrzymanie całkowania mechanicznego w czasie rzeczywistym, niezależne od
tego, czy para się związała) — **pierwotnie zostawiona nietknięta jako osobny,
ogólny mechanizm numeryczny — okazała się podlegać tej samej chorobie, i to
mocniej niż ekstrapolacja czasu kolapsu.** Historyczne uzasadnienie 529 fm
("529 fm" == `chargeCloudRestRadius`-owa `collisionBoundaryRadius`) było
praktyczne: trajektorie zbliżające się do bariery 193 fm grzęzły w sztywnym
rejonie tuż nad nią, potrzebując kroków tysiąc razy krótszych, aż budżet
zegara je cenzurował — więc 529 fm wybrano, żeby kończyć PRZED tym rejonem,
nie żeby go rozstrzygać. To uzasadnienie opierało się na tym samym błędzie
sprzężenia siatki historii ze krokiem, naprawionym w `d164f69`.

Zmierzone bezpośrednio: to samo ziarno, 200 zdarzeń, przed/po podmianie celu
na `comptonBarrierRadius` (warunkowo dla e⁺e⁻, `collisionBoundaryRadius` bez
zmian dla innych par):

| ziarno | granica | Collision | Scattering | związane (para+orto) |
|---|---|---|---|---|
| 7 | 529 fm | 72 (36,0%) | 115 (57,5%) | 13 |
| 7 | 193,3 fm | 42 (21,0%) | 141 (70,5%) | 17 |
| 42 | 529 fm | 81 (40,5%) | 107 (53,5%) | 12 |
| 42 | 193,3 fm | 43 (21,5%) | 141 (70,5%) | 16 |

Zero `NumericalFailure`/`Unresolved` w obu wariantach na obu ziarnach, koszt
zegara wyższy tylko o ~8,5% — nie "godziny", przed czym ostrzegał stary
komentarz. Prawie **połowa** zdarzeń dotąd klasyfikowanych jako "Collision"
to w rzeczywistości trajektorie, które naprawiony integrator potrafi
rozstrzygnąć dalej: większość wraca jako `Scattering`, część jako prawdziwie
związane `Para`/`Ortho-Positronium`. To zmiana głównych statystyk wyjściowych
eksperymentu dla e⁺e⁻, nie kosmetyka — dużo większa niż w CREM, gdzie
raportowany czas kolapsu prawie się nie zmienił (~0,01%). Poprawka wdrożona
identycznie: `comptonBarrierRadius` warunkowo dla e⁺e⁻, `collisionBoundaryRadius`
bez zmian dla każdej innej pary.

Produkcyjny przebieg N=1000 (ziarno 42, e⁺e⁻, bez pola zewnętrznego, 4,9
minuty) potwierdza kierunek na pełnej próbie: `Collision` 191 (19,1%),
`Scattering` 693 (69,3%), `Para-Positronium` 36 (3,6%), `Ortho-Positronium` 80
(8,0%), razem związanych 116 — mniej `Collision` niż w tabeli powyżej (193 fm
zamiast starych 529 fm), co jest spójne z kierunkiem, ale nie liczbowo
porównywalne 1:1, bo N=1000 to inny rozmiar próby niż N=200 wyżej. Stosunek
para:orto wynosi 36:80 (\(\approx1{:}2{,}2\)) wobec izotropowego 1:3 — przy
116 zdarzeniach związanych to wciąż szum statystyczny, nie odchylenie
systemowe. Zero `NumericalFailure`/`Unresolved` na całej próbie.

Ponownie przeliczone po serii ośmiu poprawek warstwy sekularnej CREM
(`ecf7380`…`0ebc7d2`: zaciśnięcie prawa siły/energii pod barierą Comptona,
regularyzowane periapsis/apoapsis/okres, zamknięty wzór na wykładnik momentu
pędu \(k(e)\) i strumień \(E_1\) jako sygnał energii zamiast przybliżenia
Darwina) — e⁺e⁻, bez pola zewnętrznego, bez ustalonego ziarna, 4,0 minuty:
`Collision` 204 (20,4%), `Scattering` 708 (70,8%), `Para-Positronium` 27
(2,7%), `Ortho-Positronium` 61 (6,1%), razem związanych 88, para:orto 27:61
(\(\approx1{:}2{,}26\)) — w granicach tego samego szumu statystycznego co
poprzedni przebieg, zero `NumericalFailure`/`Unresolved`. Ta seria poprawek
dotyczyła wyłącznie warstwy sekularnej CREM (eksperymenty 1/2), której
eksperyment 5 nie używa, więc zgodność w granicach błędu jest oczekiwana, nie
przypadkowa.

Kolejna, czwarta z rzędu poprawka tej samej warstwy (`f477866`…`ff48f08`)
usunęła ostatnie źródło niepewności w tym pomiarze: kink w pierwszej pochodnej
siły dokładnie na granicy zacisku (`clampedSeparationVector`, twardy
`max(r,floor)`), który blokował adaptacyjny integrator dokładnie tam, gdzie
orbita pierwszy raz przekraczała podłogę. Zastąpiony gładkim zmiękczeniem
Plummera (`r_eff=sqrt(r²+floor²)`, ciągłym w pierwszej pochodnej wszędzie) —
zweryfikowane bezpośrednio: ta sama trajektoria (ziarno 5), która wcześniej
grzęzła na checkpoincie 27, po zmianie dochodzi do checkpointu 32, trzy razy
głębiej (\(a=63\) fm).

Usunięcie kinku otworzyło pytanie zadane wprost: skoro integrator już nie
grzęźnie na granicy, jak głęboko naprawdę bezpiecznie sięga próg
`finalApproachMultiple` (dotąd \(10\), patrz wyżej), i czy można go
zunifikować z granicą Collision (`comptonBarrierRadius`) używaną gdzie
indziej? Zmierzone bezpośrednio na tych samych sześciu ziarnach:
`finalApproachMultiple=1.0` (zatrzymanie dokładnie na barierze) wychodzi
czysto na wszystkich sześciu, `0.5` (już za barierą) zawodzi na 2 z 6. `1.0`
wyglądał więc na zmierzoną granicę bezpieczeństwa, zgodną co do rzędu
wielkości z granicą Collision — i tak też został zacommitowany (`898dc95`).

Ta sama unifikacja zastosowana do granicy Collision eksperymentów 3/4 (dotąd
`nuclearCutoff`=10 fm, gola stała bez wyprowadzenia fizycznego, inny rząd
wielkości niż `comptonBarrierRadius`=193,3 fm) ujawniła po drodze osobny
błąd: pierwsza próba sprzęgła siatkę próbkowania eksperymentu 3
(`defaultImpactMaximum`) z tą samą granicą, co poszerzyło pole próbkowania
~361x i rozcieńczyło wszystkie kategorie wyniku na stałym N (`captured`
376→6). Naprawione rozdzieleniem: siatka próbkowania zostaje
`nuclearCutoff`-based (bez zmian, własna decyzja strojenia statystyki),
granica Collision służy wyłącznie klasyfikacji i referencji analitycznej
(`340b67b`). Sprawdzone też, czy istnieje bezpieczny punkt pośredni między
`nuclearCutoff` a `comptonBarrierRadius`: średnia geometryczna (\(44\) fm)
daje w eksperymencie 3 tylko 2% `collision`, ale ta sama głębokość w CREM
zawodzi na 1 z 3 przetestowanych ziaren — dokładnie ta sama strefa załamania
przybliżenia Darwina opisana niżej. Żaden punkt pośredni nie okazał się
bezpieczny; `comptonBarrierRadius` pozostał jedyną zmierzoną bezpieczną
wartością.

Pełny bieg produkcyjny N=1000 dla `finalApproachMultiple=1.0` ujawnił jednak,
że sześć ręcznie wybranych ziaren nie było reprezentatywną próbą: 236 błędów
numerycznych + 223 cenzurowanych budżetem zegara na 1000 (46% bez czystego
wyniku) — prawdopodobieństwo trafienia 6/6 czystych przy prawdziwym ~54%
wskaźniku sukcesu to tylko ~2,6%. Zdiagnozowane wprost (instrumentacja w
miejscu awarii): wszystkie zmierzone `NumericalFailure` miały stosunek
okres/czas-przelotu-światła między \(37{,}8\) a \(71{,}6\) — znana strefa
załamania przybliżenia Darwina — ale około połowa tych awarii wystąpiła
WYRAŹNIE POWYŻEJ `comptonBarrierRadius` (periapsis \(204\)–\(676\) fm), bo ta
granica zależy od momentu pędu \(L\) orbity, nie od promienia: orbity o
niskim \(L\) łamią przybliżenie retardacyjne przy większym \(r\) niż orbity o
wysokim \(L\). Żaden stały promień nie mógł więc poprawnie uchwycić tej
granicy dla całej populacji trajektorii.

Naprawione zastąpieniem stałego `finalApproachMultiple` warunkiem podwójnym:
pętla sekularna zatrzymuje się, gdy KTÓREKOLWIEK z dwóch niezależnych
ograniczeń zostanie osiągnięte pierwsze — periapsis \(\le\)
`comptonBarrierRadius` (granica ważności modelu punktowej cząstki) LUB
okres/czas-przelotu-światła \(\le150\) (granica ważności przybliżenia
Darwina, zależna od \(L\), margines ~2x nad zmierzonym zakresem awarii)
(`ff48f08`). Zweryfikowane: N=1000 dla obu kanałów CREM, zero błędów
numerycznych, zero cenzury w obu — pierwszy raz jednoczesna pełna unifikacja
granicy Collision między CREM (eksperymenty 1/2) i wszystkimi trzema
wariantami eksperymentów rozproszeniowych (3/4/5), z każdą trajektorią
schodzącą tak głęboko, jak jej własny moment pędu na to bezpiecznie pozwala.

Produkcyjny bieg N=1000, wszystkie pięć eksperymentów, e⁺e⁻, bez pola
zewnętrznego, na finalnym stanie tej serii (`f477866`…`ff48f08`): eksperyment
1 (p-Ps CREM, 14m47s) 1000/1000 czystych, mediana Kaplana-Meiera \(30{,}76\)
ps, średnia \(33{,}62\pm0{,}44\) ps; eksperyment 2 (o-Ps CREM, 14m58s)
1000/1000 czystych, mediana \(31{,}73\) ps, średnia \(34{,}12\pm0{,}44\) ps —
statystycznie ta sama liczba w obu kanałach, jak poprzednio. Eksperyment 3
(wiązka, ostry fokus, 5m56s): `escaped`=725, `collision`=137, `captured`=137,
`unresolved`=1, zero awarii — `sigma(reach cutoff)`=\(446288\pm35421\) barn
wobec referencji czysto kulombowskiej \(438405\) barn, w granicach błędu.
Eksperyment 4 (wiązka, szeroki kąt, 2m53s): `escaped`=1000/0/0/0/0, bez zmian
względem każdego wcześniejszego przebiegu — potwierdza, że cała ta seria
poprawek nie dotyka geometrii próbkowania eksperymentu 4. Eksperyment 5
(oddziaływania, 4m26s): `Collision`=230 (23,0%), `Scattering`=654 (65,4%),
`Para-Positronium`=33 (3,3%), `Ortho-Positronium`=81 (8,1%), `Unresolved`=2 —
para:orto \(33{:}81\approx1{:}2{,}45\), w granicach szumu statystycznego
oczekiwanego \(1{:}3\) przy 114 zdarzeniach związanych. Walidacja pól 33/33
na całej serii.

Pytanie, czy `collisionRadius` jest tym samym promieniem dla dosłownego
`Collision` i dla par związanych, doprowadziło do znalezienia osobnego,
głębszego błędu w gałęzi ekstrapolacji czasu kolapsu schwytanych par w
eksperymencie 5 (blok liczący "Bound-state CREM collapse time"). Wzór
domknięty \(t=(-E/3P)(1-r_*^3/a^3)\) ekstrapolował, aż PÓŁOŚ WIELKA \(a\)
osiągnie \(r_*=\)`comptonBarrierRadius`, podczas gdy dosłowne `Collision`
(i CREM w exp1/2) testują PERYAPSIS. Dla mocno ekscentrycznych schwytanych
par (zmierzone \(e=0{,}92\)–\(0{,}99995\), mediana \(0{,}9955\) — sam kod
już to wie, komentarz przy bramce mocy mówi wprost "e ~ 0.999") peryapsis
\(=a(1-e)\) jest \(180\)–\(20000\times\) mniejsze od \(a\), więc formuła
kazała czekać, aż \(a\) skurczy się do \(193{,}3\) fm — co odpowiada
peryapsis rzędu \(0{,}1\)–\(1\) fm, głęboko pod skalą jądrową. Zmierzone
wprost (N=1000): przy stałym mimośrodzie prawdziwe przecięcie bariery przez
peryapsis następuje, gdy \(a\) jest jeszcze w zakresie \(2400\)–\(3{,}9\cdot
10^6\) fm (mediana \(\approx43000\) fm) — para, która realnie by tam
dotarła, zostałaby wcześniej sklasyfikowana jako `Collision` przez
niezależny test na bieżącej separacji, więc niepoprawiony wzór liczył czas
do celu, którego para nigdy dynamicznie by nie osiągnęła.

Naprawione podstawieniem celu w \(a\): \(a_*=r_*/(1-e)\), z \(e\) liczonym
z bieżących \((E,L)\) pary i zamrożonym na czas ekstrapolacji — to samo
uproszczenie "zamrożonej orbity", jakie niepoprawiony wzór już stosował dla
energii (moment pędu też nigdy nie jest ewoluowany, tylko energia). Przy
okazji zunifikowane wszystkie cztery niezależne miejsca kodu, które osobno
wypisywały ten sam warunkowy wzór (`isPositronium(pair) ?
comptonBarrierRadius : collisionBoundaryRadius`) w jeden helper
`pointParticleBoundaryOf(pair)` (`particle_species.hpp`) — dokładnie ten sam
wzorzec rozjazdu, który już raz naprawiono w `340b67b`.

Zbadany też osobno warunek ratio okres/czas-przelotu-światła (`ff48f08`)
chroniący mechaniczny CREM — na próbie N=1000 (89 zdarzeń związanych) żadne
nie zbliżało się do strefy załamania Darwina przy swoim bieżącym peryapsis
(zmierzony stosunek minimalny \(1936\), próg \(150\)), więc ten konkretny
warunek nie jest dziś żywym ryzykiem w tej gałęzi — ryzykiem był cel w
\(a\), nie brak strażnika ratio. Po naprawie liczba zdarzeń z policzonym
czasem kolapsu wzrosła z \(3\)–\(4\) do \(8\) na próbie porównywalnej
wielkości, choć w samym finalnym biegu produkcyjnym (patrz niżej) akurat
żadne z \(93\) zdarzeń związanych nie przekroczyło pełnego okresu Keplera w
oknie obserwacji — statystyka tej bramki jest osobna od naprawionego celu i
zależy od losowego ziarna, nie od tej poprawki.

Ponownie przeliczone na produkcyjnym przebiegu N=1000 dla wszystkich pięciu
eksperymentów po tej naprawie: eksperyment 1 (p-Ps CREM, 14m26s) 1000/1000
czystych, mediana Kaplana-Meiera \(30{,}26\) ps, średnia \(33{,}55\pm0{,}44\)
ps; eksperyment 2 (o-Ps CREM, 15m04s) 1000/1000 czystych, mediana
\(30{,}78\) ps, średnia \(33{,}51\pm0{,}44\) ps. Eksperyment 3 (5m56s):
`escaped`=739, `collision`=131, `captured`=129, `unresolved`=1, zero
awarii — `sigma(reach cutoff)`=\(426742\pm34757\) barn wobec referencji
czysto kulombowskiej \(438405\) barn. Eksperyment 4 (2m55s):
`escaped`=1000/0/0/0/0, bez zmian. Eksperyment 5 (4m36s): `Collision`=216
(21,6%), `Scattering`=690 (69,0%), `Para-Positronium`=17 (1,7%),
`Ortho-Positronium`=76 (7,6%), `Unresolved`=1 — para:orto
\(17{:}76\approx1{:}4{,}5\), w granicach szumu statystycznego przy zaledwie
93 zdarzeniach związanych (mniej niż w poprzednich biegach — losowa
fluktuacja rozmiaru próby, nie efekt tej poprawki, która nie dotyka
klasyfikacji Collision/Scattering/para/ortho, tylko ekstrapolację czasu
kolapsu już sklasyfikowanych par). Zero awarii na całej serii, walidacja
pól 33/33.

Sprawdzone też diagnostycznie (bez zmiany produkcyjnej), gdzie dziś leży
NASTĘPNA ściana numeryczna mechanicznej pętli CREM, gdyby oba warunki
zatrzymania (promień i ratio) tymczasowo poluzować do nieosiągalnie
głębokich wartości, zostawiając samą podłogę zmiękczenia
(`separationFloor()`) bez zmian, na produkcyjnym `comptonBarrierRadius`. Na
10 ziarnach pętla dochodzi teraz aż do peryapsis \(23\)–\(383\) fm (mediana
\(56\) fm, średnia \(92\) fm) i stosunku okres/czas-przelotu-światła
\(13\)–\(52\) (mediana \(20\), średnia \(23\)) — dopiero tam napotyka
prawdziwą awarię numeryczną. To \(6\)–\(12\times\) głębiej w skali ratio niż
produkcyjny próg zatrzymania \(150\): margines, jaki `ff48f08` przyjął
(~2× nad zmierzonym wtedy zakresem awarii \(37{,}8\)–\(71{,}6\)), okazuje
się w praktyce jeszcze bardziej konserwatywny niż zakładano — prawdziwa
ściana integratora leży wyraźnie głębiej niż strefa, w której produkcyjny
warunek już się zatrzymuje.

To NIE jest propozycja obniżenia progu produkcyjnego. `comptonBarrierRadius`
nie jest marginesem numerycznym do przycinania jak dawne
`finalApproachMultiple` — to wyprowadzona fizycznie granica, przy której
klasyczna elektrodynamika punktowej cząstki przestaje obowiązywać (patrz
"Podłoga zasięgu CREM" wyżej). Próg \(150\) już dziś nie kosztuje nic
produkcyjnie (N=1000 dla obu kanałów CREM: zero awarii, zero cenzury), więc
nie ma powodu go zawężać tylko dlatego, że silnik ma jeszcze zapas. Wynik
jest tu wyłącznie diagnostyczny: potwierdza, że zmiękczenie Plummera
(`f477866`) i warunek ratio (`ff48f08`) razem dają silnikowi solidny bufor
poniżej punktu, w którym model i tak przestaje deklarować ważność.

**Ta sama sonda powtórzona z pytaniem "czy da się zejść w stronę 10 fm
(`nuclearCutoff`)", zamiast tylko notować, gdzie akurat pada.** Cel obu
warunków zatrzymania retargetowany na `nuclearCutoff` zamiast
`comptonBarrierRadius` (`separationFloor()` znów bez zmian), z tym samym
zamiarem diagnostycznym co wyżej — i za pierwszym podejściem wynik był
identyczny co poprzednio: prawdziwa awaria numeryczna, nie licząc miejsca
zatrzymania. Ale `maximumDepth` mechanicznego integratora (budżet
podpodziału kroku w pojedynczej mierzonej orbicie, `crem_trajectory.hpp`)
od zawsze siedzi na `12`, ustalonym w audycie SPRZED naprawy sprzężenia
historii retardowanej z krokiem (`d164f69`) — tamten audyt zmierzył, że
podniesienie tego budżetu kosztuje `445×` czasu i **nadal** kończy się
awarią, więc nikt nie sprawdził go ponownie po tamtej naprawie.

Sprawdzone teraz: `maximumDepth=20` (z `12`), 10 ziaren, budżet zegarowy do
250-350 s na próbę. Wynik zaskakująco pozytywny — mediana głębokości, na
jakiej pojawia się prawdziwa awaria (lub próba wyczerpuje budżet, wciąż w
trakcie zbiegania, bez żadnej awarii), spadła z \(56\)-\(92\) fm
(mediana/średnia poprzedniej sondy) do przedziału \(13\)-\(20\) fm na
wszystkich dziesięciu próbach jednocześnie — \(5\)-\(6\times\) głębiej, i to
konsekwentnie, nie tylko w najlepszym przypadku. Dwa ziarna dotarły do
prawdziwej awarii numerycznej przy \(19{,}17\) i \(13{,}47\) fm; dwa
wyczerpały wewnętrzny budżet zegarowy (nie awaria) przy \(20{,}1\)-\(20{,}3\)
fm; pozostałych sześć ubiła zewnętrzna granica czasu tego pomiaru w trakcie
liczenia KOLEJNEGO punktu kontrolnego, z ostatnim zarejestrowanym peryapsis
\(13{,}4\)-\(17{,}1\) fm i bez żadnego śladu awarii — a więc prawdopodobnie
zaszłyby jeszcze głębiej z cierpliwszym budżetem. Sprawdzone też: podniesienie
`maximumDepth` dalej, z `20` na `28`, na najgorszym ziarnie z tej serii nie
zmieniło NIC — identyczna awaria, identyczna głębokość, co do ostatniej
cyfry — więc `20` już nasyca to, co ten konkretny budżet może dać; głębsza
awaria wymagałaby czegoś innego niż więcej podpodziałów. Sprawdzone też
`--integrator-order 4` na dwóch najgorszych ziarnach: pogorszenie na obu
(awaria płytsza, czas \(3\)-\(4\times\) dłuższy) — potwierdza wcześniejsze
ustalenie, że rząd kompozycji nie pomaga tutaj.

Koszt: \(5\)-\(20\times\) więcej czasu na próbę (z \(8\)-\(15\) s do
\(40\)-\(260\)+ s), bo każdy kolejny punkt kontrolny bliżej dna kosztuje
więcej niż poprzedni — obserwowane bezpośrednio: jeden punkt kontrolny dla
ziarna 11 sam zajął \(87\) s. To NADAL nie jest propozycja zmiany czegokolwiek
produkcyjnego — `comptonBarrierRadius` pozostaje fizyczną, nie numeryczną,
granicą, a próg \(150\) już ma zerowy koszt produkcyjny przy obecnym
`maximumDepth=12`. Wynik jest wart odnotowania z innego powodu: silnik,
po naprawie sprzężenia historii z krokiem, ma jeszcze WIĘCEJ zapasu
numerycznego niż `4c48ffe` sądziło — jedna nietknięta od dawna stała
(`maximumDepth`) sama, bez żadnej innej zmiany, przesunęła zmierzoną ścianę
o rząd wielkości bliżej `nuclearCutoff`, a wcześniejsza ocena tej samej
gałki ("445× kosztu i nadal pada") była prawdziwa tylko względem
NIEISTNIEJĄCEGO już bugu, nie samego integratora. Sonda w
`crem_collapse.hpp`/`crem_trajectory.hpp` cofnięta, `positronium_validation`
33/33 bez zmian.

**Scharakteryzowane dokładnie, co konkretnie zawodzi na dnie tego przedziału
(13–20 fm).** Dodana chwilowa instrumentacja rozróżniająca trzy niezależne
gałęzie awarii w `crem_collapse.hpp` (błędna osculacja na starcie punktu
kontrolnego, sama mechaniczna orbita pomiarowa, strażnik "nie więcej niż
50% energii na jedną orbitę") oraz dokładne miejsce poddania się rekursji w
`ClassicalTrajectoryEngine::advanceAdaptive` (`crem_engine.hpp`). Wynik na
dwóch niezależnych ziarnach (1, 23): za każdym razem ta sama gałąź —
mechaniczna orbita pomiarowa zwraca `NumericalFailure`, ale jej stan
końcowy jest **w pełni skończony** (`runFinalFinite=1`), a poddanie się
zdarza się niemal natychmiast (`run.elapsedTime~10⁻²³` s, dosłownie na
pierwszym kroku tej orbity). To NIE jest przepełnienie zmiennoprzecinkowe
ani rozbieżność fizyki. Strażnik straty energii na orbitę (próg 50%) i
próg periapsis na starcie punktu kontrolnego — obie inne gałęzie awarii —
ani razu się nie uruchomiły.

Bezpośredni odczyt z punktu poddania rekursji (ziarno 23, periapsis
13,3 fm): `depth=20/20` — cały przyznany budżet podpodziału **wyczerpany**
— przy `error=3,046·10⁻⁵` wobec progu `1·10⁻⁵`: **zaledwie 3× za wysoko**,
nie o rzędy wielkości. Precyzja zmiennoprzecinkowa nie jest tu winna: krok
w chwili poddania (`dt=3,28·10⁻³⁰` s) jest wciąż `9·10⁸` razy większy od
rozdzielczości `double` przy tej wartości czasu symulacji. Przyczyna jest
prostsza: tuż nad tą głębokością częstość orbitalna \(\omega=\sqrt{k/(\mu
r^3)}\) rośnie tak szybko, że sam KROK ZEWNĘTRZNY (`2\pi/(128\omega)`)
kurczy się do rzędu `10⁻²⁴` s, zanim jeszcze zacznie się podpodział — a
osiągnięcie tolerancji `10⁻⁵` na skomponowanym kroku Yoshidy przy tej
sztywności wymaga więcej niż 20 połówkowań, nie z powodu jakiejś
osobliwości, tylko dlatego, że błąd lokalny maleje z krokiem *wolniej* niż
budżet głębokości rośnie.

Potwierdzone bezpośrednio, że to naprawdę tylko kwestia budżetu, nie
twardej ściany: `maximumDepth=24` na tym samym ziarnie (23) **przetrwało
znacznie dłużej** niż `maximumDepth=20` na identycznym punkcie startowym —
ale kosztem, który eksploduje szybciej niż liniowo: sam JEDEN kolejny
punkt kontrolny pochłonął ponad 100 s obliczeń (więcej niż cała reszta
przebiegu do tego miejsca razem wzięta) i nie zdążył się rozstrzygnąć
nawet przy budżecie zegarowym 500 s, więc próba została przerwana bez
czystego wyniku. Wcześniejsze zaobserwowane plateau (`28` nie poprawiało
nic ponad `20` na innym ziarnie, 7) nie jest więc uniwersalną granicą
integratora — jest specyficzne dla orbity: niektóre podejścia do periapsis
są "sztywniejsze" od innych i wymagają więcej połówkowań, by osiągnąć tę
samą tolerancję.

Wniosek: granica 13–20 fm nie jest pojedynczą ścianą fizyczną ani
numeryczną w konkretnym promieniu — to punkt, w którym koszt utrzymania
stałej tolerancji `10⁻⁵` metodą czystego połówkowania kroku zaczyna rosnąć
kombinatorycznie, nie punkt, w którym staje się to dosłownie niemożliwe.
Głębsze zejście wymagałoby czegoś jakościowo innego niż większy
`maximumDepth` przy tej samej tolerancji — np. zmiennych regularyzujących
bliski przelot (podobnie jak `osculatingPeriapsis`/`regularizedPeriod` już
robią dla samej sekularnej ekstrapolacji) albo tolerancji łagodniejszej w
tym jednym reżimie — nie jest to jednak zmiana warta wprowadzania: cel
pozostaje wyłącznie diagnostyczny, `comptonBarrierRadius` (193,3 fm) leży
z dużym zapasem powyżej całego zbadanego tu przedziału. Cała
instrumentacja (`crem_collapse.hpp`, `crem_engine.hpp`,
`crem_trajectory.hpp`) cofnięta, `positronium_validation` 33/33 bez zmian.

**Dodany, na stałe: `--radiation-reaction stochastic` — kwantowana emisja
zamiast ciągłego hamowania.** Odpowiedź na pytanie, czy dyskretyzacja
promieniowania (fotony zamiast gładkiej mocy Larmora) da się w ogóle wstawić
do tego silnika. Mechanizm: ta sama moc dipola E1, którą
`leadingElectricDipolePower` już liczy dla każdego modelu, jest teraz
bankowana jako HARACZ procesu Poissona (całka `moc/(ħω_orb)dt`, ω_orb —
ta sama chwilowa częstość oskulacyjna, której już używa dobór kroku) zamiast
być odejmowana jako ciągła siła. Próg następnego fotonu losowany jest
metodą odwróconego dystrybuanty (Exp(1) z 64-bitowego strumienia
`splitMix64`-podobnego, zasianego bitami stanu początkowego — deterministyczne,
odtwarzalne z `--seed`, bez nowego parametru w sygnaturach funkcji). Między
fotonami `chargeReaction` tego modelu jest **dokładnie zerowe** — para
porusza się po gołej trajektorii wzajemnej siły Lorentza, dokładnie
zachowując energię mechaniczną, tak jak prawdziwy emiter nie traci energii
w sposób ciągły między skokami kwantowymi. Sam foton to kopnięcie
zachowujące całkowity pęd pary, usuwające dokładnie `ħω` z energii kinetycznej
ruchu względnego (ten sam, zredukowanej-masy, dwuciałowy obraz, którego
`crem_collapse.hpp` już używa dla "energii orbity"); pęd odrzutu fotonu
(`ħω/c`) nie jest osobno modelowany — przy prędkościach, do których ten tryb
jest pomyślany, jest o rzędy wielkości mniejszy niż już nakładana zmiana
prędkości.

**Sprawdzone bezpośrednio, zanim cokolwiek uznano za działające: czy
podwójnie nie liczy energii.** Pierwsza wersja kredytowała usuniętą energię
też do `s.radiatedEnergy`/`orbitalRadiatedEnergy` — okazało się to błędem:
te pola są już, bezwarunkowo i niezależnie od modelu reakcji, wypełniane
przez kwadraturę strumienia Poyntinga w `integrateElectrodynamicStep`
(prawdziwy, zmierzony strumień pola tej konkretnej trajektorii — nie
przestaje płynąć tylko dlatego, że `chargeReaction` jest zerowe: para nadal
przyspiesza pod samą siłą Lorentza, a równania Maxwella nie wiedzą nic o
księgowości tego trybu). Naprawione: kopnięcie usuwa WYŁĄCZNIE energię
mechaniczną; niezależny pomiar strumienia zostaje nietknięty, dokładnie jak
dla pozostałych modeli, gdzie oba są already porównywane, nie scalane
(`reaction/flux` diagnostyka).

**Znaleziona luka ekspozycji — i domknięta.** Ręcznie policzony haracz dla
pojedynczej orbity przy a=3,1 pm (checkpoint 16 ziarna 42) wynosi `5,5·10⁻⁵`
— więc pojedyncza zmierzona orbita praktycznie nigdy nie przekracza progu
Exp(1). `crem_collapse.hpp` woła `runMechanicalTrajectory` osobno dla
KAŻDEGO punktu kontrolnego, licząc tylko JEDNĄ orbitę na pomiar, po czym
analitycznie PRZESKAKUJE do 200 000 kolejnych — a pierwszy akumulator
haraczu (ten w `crem_trajectory.hpp`) żyje tylko wewnątrz jednego wywołania
i znika, gdy ono się kończy, więc nigdy nie widzi tych 200 000 pominiętych
orbit, gdzie żyje niemal cały prawdziwy haracz.

Domknięte dodaniem DRUGIEGO, niezależnego akumulatora haraczu — tym razem w
`crem_collapse.hpp` samym, całkującego analitycznie tę samą stawkę
`moc/(ħω)` po całym pominiętym odcinku, w zamkniętej formie względem
klasycznej obwiedni `u(n)=u₀(1−Jx)^(-2/3)`, której i tak już używa sam skok
sekularny (J = `jumpParameter`, x = n/orbitsToSkip ∈ [0,1]):

\[
\text{haracz}(J)=\frac{P_0 T_0}{\hbar\omega_0}\cdot n_{\rm skip}\cdot
\frac{3}{J}\Bigl(1-(1-J)^{1/3}\Bigr)
\]

— zweryfikowane względem siłowej kwadratury numerycznej do `10⁻¹²`
względnie w całym zakresie `J∈[0,0.3]`, jaki produkcja kiedykolwiek
generuje. Moment i wielkość każdego fotonu w obrębie skoku wyznacza
odwrócenie tej samej zamkniętej formy (też zweryfikowane numerycznie do
maszynowej precyzji); foton trafia bezpośrednio w elementy oskulacyjne
(kopnięcie energii plus ten sam zamknięty stosunek `k=-(1-e²)/(2+e²)`,
którego już używa ciągła aktualizacja momentu pędu) — bo podczas
pominiętego odcinka nie istnieje żaden mechaniczny `State`, tylko te dwie
liczby. Dla pozostałych modeli zachowanie jest bit w bit identyczne (ta sama
gałąź `else` co zawsze).

**To domknięcie ujawniło coś więcej niż samą lukę — prawdziwy efekt
fizyczny.** `ħω_orb` w skali samego pozytonium jest porównywalne albo
większe od CAŁEJ jego energii wiązania przez większość kolapsu (~9-13 eV na
foton blisko startu, wobec ~6,8 eV wiązania) — pojedyncze fotony NIE są
małym zaburzeniem klasycznej krzywej, tak jak zakłada własne wyprowadzenie
`u₀(1-Jx)^(-2/3)` w reżimie małego haraczu: są duże, rzadkie, skokowe, a
TOR między nimi jest dokładnie zachowywany (patrz akapit wyżej). Zmierzone
bezpośrednio (ziarno 42): modele ciągłe/klasyczne dają ~36-40 ps; ten model
daje **665 ps** dla tego samego ziarna, a na 10 ziarnach medianę **149 ps**
i średnią **276±91 ps** (sigma/średnia≈1,0 — czyli naprawdę szeroki rozkład,
nie ciasny wokół odpowiedzi klasycznej). To rzeczywista, rzędu wielkości
różnica wynikająca z tej konkretnej dyskretyzacji emisji, nie artefakt
zaokrąglenia.

W przeciwieństwie do sondy SED (`--zpf`), która zamknęła się wynikiem
negatywnym bez jakościowej różnicy, ten eksperyment DAJE jakościowo inną,
fizycznie interpretowalną odpowiedź — i to jest sedno sprawy: naiwna
dyskretyzacja emisji rozmiarem fotonu "ħ razy chwilowa częstość klasyczna"
nie redukuje się gładko do odpowiedzi klasycznej, gdy skala fotonu
przestaje być mała względem własnej skali energetycznej układu — dokładnie
ten sam mechanizm, przez który stara (sprzed 1925 r.) półklasyczna teoria
kwantowa okazała się niewykonalna, i dokładnie to, o co pytał użytkownik na
początku tego wątku.

Zweryfikowane: kompilacja czysta, `positronium_validation` 33/33 bez zmian.
Test na 10 ziarnach: 10/10 trajektorii dochodzi do granicy (100% ukończenia,
zero cenzury/awarii) — mechanizm jest stabilny, nie tylko poprawny
punktowo.

**Na wyraźną prośbę: `stochastic` zastąpił `individual` jako domyślny
model produkcyjny (stały inicjalizator `gRadiationReactionModel` w
`positronium.cpp`), nie tylko opcjonalny tryb obok niego.** To odwraca
dotychczasową konwencję tego projektu, w której każdy nowy/eksperymentalny
mechanizm (`--zpf`, `--beam-energy-sigma-ev`) domyślnie zostawał wyłączony
— świadoma decyzja, potwierdzona wprost, nie domyślne zachowanie tej
sesji. Konsekwencja jest realna i warta podkreślenia osobno: **każda
liczba czasu kolapsu cytowana gdziekolwiek wcześniej w tym dokumencie
(rzędu 36-40 ps dla e⁺e⁻, i analogiczne dla mionium/protonium) została
zmierzona pod starym domyślnym `individual`, nie pod obecnym
`stochastic`.** `./positronium` bez żadnych flag daje dziś medianę rzędu
149 ps i znacznie szerszy rozkład (sigma/średnia≈1,0) zamiast ciasnych
~36-40 ps. Stary model nie zniknął — `--radiation-reaction individual`
odtwarza go bit w bit (zweryfikowane: 35,99 ps, identyczne jak przed tą
zmianą) — zmienił się wyłącznie *domyślny wybór* dla kogoś, kto nie poda
żadnej flagi. `positronium_validation` 33/33 niezależnie od tej zmiany
(sama walidacja pól nie zależy od globalnego domyślnego modelu reakcji);
wszystkie 5 eksperymentów sprawdzone dymnie pod nowym domyślnym modelem —
brak awarii, sensowne wyjścia.

**Pełny przebieg produkcyjny N=1000, wszystkie 5 eksperymentów, e⁺e⁻, bez
pola zewnętrznego, ziarno 42, ustawienia domyślne (czyli `stochastic`) —
łącznie ~1 godz. 32 min na 4 wątkach.**

*Eksperymenty 1/2 (kolaps CREM, p-Ps/o-Ps).* 905/1000 (p-Ps) i 899/1000
(o-Ps) trajektorii dotarło do granicy przy domyślnym budżecie zegarowym
20 s/zdarzenie — reszta ocenzurowana tym budżetem, zero
`NumericalFailure` w obu kanałach. Mediana Kaplana-Meiera: p-Ps 118,4 ps,
o-Ps 119,4 ps — statystycznie ta sama liczba w obu kanałach, jak zawsze.
Średnia ukończonych/RMST: p-Ps 161,7±5,7 / 183,3±6,8 ps, o-Ps
162,0±5,7 / 184,5±6,9 ps, rozrzut sigma/średnia ≈1,06 w obu — around
3-5× szerszy rozkład i 3-5× dłuższy medianowy czas niż pod starym
domyślnym `individual` (30,8–31,9 ps), dokładnie zgodnie z tym, co
akapity o `stochastic` wyżej już przewidziały i zmierzyły na mniejszej
próbie. 90,5%/89,9% kompletacji jest poniżej progu 100%, więc — jak samo
narzędzie ostrzega — mediana i RMST/średnia obu kanałów to zakres
wiarygodności, nie liczba bez obciążenia; podniesienie
`--crem-wallclock-budget-s` powyżej domyślnych 20 s zmniejszyłoby cenzurę
kosztem dłuższego przebiegu.

*Eksperyment 3 (kanał krótkiego zasięgu).* Escaped 775, collision 138,
captured 86, unresolved 1, zero `failed`. Dopasowanie kształtu
Rutherforda \(C_R=0{,}691\pm0{,}015\) (95% Wilson \([0{,}662;\,0{,}719]\))
— zauważalnie DALEJ od czystego Coulomba niż pod starym domyślnym modelem
(tam \(C_R\) był bliżej jedynki dla tego kanału). Niezależna reszta
`|dE_reaction-vs-flux|/E_rad` ma medianę \(0{,}999988\) — praktycznie
100% rozjazdu. To nie zaskoczenie: ta reszta porównuje pracę modelu
Landaua-Lifshitza z niezależnie zmierzonym strumieniem, a `stochastic` nie
używa ciągłej siły LL wcale (`chargeReaction` zerowe między fotonami), więc
to porównanie mierzy teraz różnicę między dwoma różnymi modelami fizyki,
nie błąd księgowania — to samo zastrzeżenie o zakresie stosowalności
eksperymentu 3, który już wcześniej był poza reżimem modelu, dotyczy tu
tym bardziej.

*Eksperyment 4 (rozpraszanie szerokokątowe).* Czysty: 1000/1000 escaped,
zero collision/captured/unresolved/failed. \(C_R=0{,}966\pm0{,}024\) (95%
Wilson \([0{,}918;\,1{,}012]\)) — **statystycznie nierozróżnialne od
pomiaru pod starym domyślnym `individual`** (tam też \(C_R=0{,}966\pm
0{,}024\)). Zgodne z oczekiwaniem: przy rozpraszaniu szerokokątowym
trajektorie nigdy nie schodzą wystarczająco głęboko, by wybór modelu
reakcji promieniowania cokolwiek zmienił.

*Eksperyment 5 (interakcja/klasyfikacja).* Collision 216 (21,6%),
Scattering 691 (69,1%), Para-Ps 28 (2,8%), Ortho-Ps 65 (6,5%) — razem 93
zdarzeń związanych, para:orto = 28:65, blisko izotropowego oczekiwania
1:3. Rozkład klasyfikacji nie zależy silnie od modelu reakcji (klasyfikacja
zapada wcześnie, nie podczas długiego kolapsu), więc te udziały procentowe
pozostają zgodne z wcześniejszymi pomiarami pod `individual` w granicach
szumu statystycznego jednej próby N=1000.

Wykresy wszystkich pięciu eksperymentów (`distributions/*.pdf`)
przeliczone i zacommitowane razem z tym akapitem.

**Dopisany kierunek fotonu, na zadane pytanie: "czy fotony są wysyłane w
przypadkowych kierunkach?"** Odpowiedź brzmiała: nie były wysyłane w ŻADNYM
kierunku — `OsculatingElements` przechowywał moment pędu jako gołą liczbę
(`double specificAngularMomentum`), bez orientacji, więc kopnięcie zmieniało
tylko *wielkości* E i L, nigdy kierunek. Dodane: śledzenie orientacji orbity
jako wektora jednostkowego (`angularMomentumDirection`, inicjalizowany z
`noetherAngularMomentum` rzeczywistego wylosowanego stanu początkowego),
oraz fizycznie poprawny rozkład kątowy emisji — nie sugerowany przez
użytkownika spin (S=0/1 steruje w tym modelu czym innym: wyrównaniem
dipoli/regułami wyboru anihilacji, to inna oś niż orbitalny moment pędu),
tylko prawdziwy wzór promieniowania **wirującego** dipola E1:
\(dP/d\Omega\propto(1+\cos^2\theta)\) względem osi momentu pędu — maksimum
WZDŁUŻ osi, nie w płaszczyźnie orbity (to byłby wzór \(\sin^2\theta\)
pojedynczego, liniowo oscylującego dipola, inny przypadek fizyczny).

Odwrócenie dystrybuanty tego rozkładu sprowadza się do sześciennego
równania \(\mu^3+3\mu+(4-8u)=0\) (\(\mu=\cos\theta\)) — rozwiązanego w
zamkniętej formie wzorem Cardana, zweryfikowanym numerycznie względem
dystrybuanty do \(10^{-15}\) bezwzględnie na całym zakresie \(u\in[0,1)\).
Składowa odrzutu wzdłuż osi momentu pędu przechyla płaszczyznę orbity — to
dokładnie \(\delta L=\mathbf r\times\delta\mathbf v\), którego ta
reprezentacja (same elementy oskulacyjne, bez anomalii prawdziwej) nie
umie policzyć dokładnie, bo nie wie, W KTÓRYM miejscu orbity foton akurat
wystrzelił. Rozwiązanie w tym samym duchu co współczynnik `k` wyżej:
zastąpienie nieznanego chwilowego promienia półosią wielką (skalą
uśrednioną po orbicie), a nieznanego azymutu przechyłu — losowaniem
jednostajnym (zamiast go zgadywać) — co daje **izotropowy błądzenie
losowe** orientacji płaszczyzny, nie systematyczną precesję, dokładnie
właściwy obraz uśredniony po nieznanej fazie.

Zweryfikowane bezpośrednio: na pełnej trajektorii (ziarno 42, p-Ps) 3
fotony wystrzeliły, kąty przechyłu 0,006–0,113 rad (kilka stopni, sensowna
skala), `angularMomentumDirection` pozostał wektorem jednostkowym po
każdym obrocie Rodriguesa. Statystyka zbiorcza nietknięta w granicach szumu
innego ziarna (10 trajektorii: mediana 151,6 ps wobec 148,8 ps sprzed tej
zmiany) — jak oczekiwano, bo kierunek jest ortogonalny do już istniejącej
fizyki wielkości E/L, tylko dodaje spójną orientację, której wcześniej nie
było wcale. `positronium_validation` 33/33 bez zmian.

**Sprawdzone na zadane pytanie "czy wszystkie prawa zachowania są tu
spełnione": nie, i pęd liniowy naprawiony.** Energia — tak, dokładnie
(zamknięty rachunek, każdy foton usuwa dokładnie \(\hbar\omega\) z energii
mechanicznej i dokładnie tyle trafia do `radiatedEnergyTotal`). Pęd
liniowy — **nie był w ogóle egzekwowany**, i przy okazji złapany błąd we
własnym komentarzu: kod twierdził, że pęd fotonu \(\hbar\omega/c\) jest "o
rzędy wielkości mniejszy" niż już nakładana zmiana prędkości. Zmierzone
wprost: stosunek \(p_{foton}/p_{orbitalny}\) wychodzi dokładnie jako
(zredukowana długość fali Comptona pary)/\(a\), niezależnie od prędkości —
\(0{,}007\) na starcie orbity, ale już \(0{,}25\) przy \(a=3{,}1\) pm i
\(>1\) poniżej \(\approx0{,}77\) pm — czyli realny, znaczący ułamek na
większości głębokości, jaką ten model osiąga, w dodatku dobrze wewnątrz
własnej zadeklarowanej granicy ważności (`comptonBarrierRadius`). Nieprawdziwy
komentarz poprawiony na miejscu.

Naprawione w `crem_collapse.hpp` (nie w wersji mechanicznej
`crem_trajectory.hpp` — ta nie ma pozycji/orientacji potrzebnej do tej
poprawki i rzadko jest ćwiczona produkcyjnie). Warunki początkowe CREM są
przygotowywane przy dokładnie zerowym pędzie całkowitym (rozkład prędkości
względnej wg stosunku mas daje \(m_1v_1+m_2v_2=0\) tożsamościowo), a każdy
model ciągły utrzymuje to zero z konstrukcji (siła własna klasycznego
dipola jest zgodna z trzecią zasadą Newtona względem pola, które sama
promieniuje) — więc zachowanie pędu sprowadza się tu do prostszego pytania:
czy cały układ dostaje jednorodne kopnięcie równe i przeciwne do pędu
fotonu. Teraz dostaje: pełny trójwymiarowy kierunek fotonu (już
próbkowany dla przechyłu płaszczyzny, teraz też użyty do samego odrzutu)
przesuwa nowo dodaną, trwałą prędkość środka masy (`centreOfMassVelocity`,
start dokładnie w zerze), a wynikająca stąd zmiana energii kinetycznej
środka masy (\(v_{cm}\!\cdot\!p_{foton}+p_{foton}^2/(2M)\), policzona
dokładnie, nie odrzucona jako wyraz drugiego rzędu) jest doliczana do
budżetu orbitalnego ponad samą energię fotonu — żeby energia całkowita
(orbitalna+CM) nadal spadała dokładnie o tyle, ile trafia do
`radiatedEnergyTotal`.

Zweryfikowane bezpośrednio: na tej samej trajektorii (ziarno 42, p-Ps) 3
fotony, `cmEnergyKick/photonEnergy` rośnie z głębokością (\(6{,}3\cdot
10^{-6}\to4{,}2\cdot10^{-5}\to5{,}3\cdot10^{-4}\)) — mały, ale rosnący
efekt, dokładnie zgodny z policzonym wyżej stosunkiem pędów.
\(|v_{cm}|\) rośnie od zera do \(\approx298\) km/s (wciąż nierelatywistyczne,
\(\approx0{,}001c\)) po trzecim fotonie. Statystyka zbiorcza (10
trajektorii) nietknięta w granicach szumu: mediana \(147{,}8\) ps, średnia
\(272{,}9\pm91{,}8\) ps — statystycznie ta sama liczba co przed tą
poprawką, jak należało oczekiwać, bo poprawka jest efektem drugiego rzędu
poza najgłębszą częścią kolapsu. `positronium_validation` 33/33 bez zmian.

**Zbadano moment pędu — i przechył płaszczyzny orbity okazał się być tym
samym pędem fotonu liczonym drugi raz.** Wcześniejszy mechanizm przechylał
`angularMomentumDirection` traktując odrzut jako kopnięcie \(r\times\Delta
v\) ruchu WZGLĘDNEGO (\(r\sim\) półoś wielka, \(\Delta v\sim
p_{foton}/\mu\)) — niezależnie od poprawki pędu liniowej powyżej, która
już równoważy TEN SAM pęd fotonu jednorodnym kopnięciem obu cząstek.
Sprawdzone wprost (Python, trzy stosunki mas: \(1{:}1\), \(1836{:}1\),
\(3{:}1\)): jednorodne kopnięcie środka masy daje moment siły
\(\sum_i r_i\times(m_i\Delta v_{cm})=\Delta v_{cm}\times\sum_i m_ir_i=0\)
**dokładnie**, bo \(\sum_i m_ir_i=0\) tożsamościowo względem środka masy —
resztki rzędu \(10^{-17}\), czysty szum numeryczny. Jednorodne pchnięcie
przez własny środek masy nie może układu skręcić. Przechył wydawał więc
ten sam pęd fotonu dwa razy: raz jako poprawny odrzut środka masy (wyżej),
raz jako niezależne, nieskoordynowane z nim kopnięcie ruchu względnego —
usunięty z kodu jako błąd, nie przeważony na nowo.

Co zostaje jako **prawdziwe** źródło skręcenia płaszczyzny: nie orbitalny
moment pędu fotonu (ten wymagałby nieznanej anomalii prawdziwej, jak już
było), tylko jego **spin** — \(\pm\hbar\) wzdłuż własnego kierunku
propagacji, uniwersalny, dokładny fakt dla dowolnego bezmasowego bozonu
spinu 1, nie przybliżenie orbit-averaged. Nie zaimplementowany, bo to nie
mała poprawka: własny punkt startowy CREM to wartość Bohra/SED
\(L=\hbar\) (`physical_constants.hpp`), więc \(|L_{foton}|/|L_{orbitalny}|
=\hbar/\hbar=1\) **dokładnie** już przy pierwszym fotonie, z samej
konstrukcji warunku początkowego — nie zmierzona małość, jak przy
\(p_{foton}/p_{orbitalny}\) wyżej, tylko przeciwieństwo: zaburzenie rzędu
jedności za każdym razem, którego żadna klasyczna adiabatyczna księgowość
(łącznie ze współczynnikiem `k`) nie jest w stanie wchłonąć. Model siedzi
dokładnie na granicy kwantowej, gdzie moment pędu realnego związanego
układu jest dyskretną drabiną rozstawioną co \(\hbar\), nie ciągłym
wektorem klasycznym — żadna sztuczka księgowa tego nie zamyka.

Zweryfikowane po usunięciu przechyłu: czysta kompilacja (zero ostrzeżeń),
`positronium_validation` 33/33, test dymny (ziarno 7, o-Ps,
`--radiation-reaction stochastic`) — 3/3 trajektorii kolapsuje bez awarii
numerycznych, fotony nadal strzelają poprawnie (przykład: foton #1 zmienia
\(L\) z \(5{,}55\cdot10^{-5}\) na \(1{,}56\cdot10^{-5}\) przez współczynnik
`k`, zgodnie z wcześniej zmierzonym efektem "energia fotonu porównywalna z
energią wiązania"), bez żadnych linii `TILT` w logu. Kierunek płaszczyzny
orbity jest teraz niezmienny wobec odrzutu fotonu — spójne z poprawką pędu
liniowego, a nie uproszczenie mechanizmu, który i tak nigdy nie działał
poprawnie.

**Wyznaczony na zadane pytanie "określ współczynnik k dla nowego modelu
promieniowania" — i okazał się być cicho błędny przy niemal każdym
fotonie produkcyjnym.** Sam \(k=-(1-e^2)/(2+e^2)\) jest poprawny i
niezmieniony: to własność zależności od \(r\) klasycznej siły reakcji
dipolowej, nie zależy od tego, czy strata jest księgowana w sposób ciągły
czy w porcjach fotonowych (ten sam argument, co niezależność od modelu
reakcji, już wyprowadzony wyżej). Błąd był w tym, JAK
`crem_collapse.hpp` stosował ten współczynnik do pojedynczego fotonu:
\(L{\cdot}{=}(E_{po}/E_{przed})^k\), formuła "zamrożonego k", jest tylko
przybliżeniem pierwszego rzędu (dla małego skoku) relacji różniczkowej
\(d(\ln L)/d(\ln|E|)=k(e)\), z której pochodzi. Gałąź zbiorcza/deterministyczna
używa tego samego przybliżenia, ale jej skok na checkpoint jest
ograniczony do \(jumpParameter\le0{,}30\) (stosunek energii \(<{\sim}1{,}33\)),
gdzie błąd jest łagodny; pojedynczy foton stochastyczny potrafi unieść dużo
większą wielokrotność bieżącej energii orbitalnej za jednym razem
(zmierzone gdzie indziej w tym pliku: aż do \(\approx18{,}5\times\)).
Zmierzone wprost: przy takiej skali formuła zamrożonego k nie jest lekko
nietrafiona, jest jakościowo zła — sprawdzone na prawdziwym przebiegu
produkcyjnym (ziarno 42, o-Ps, 5 fotonów): **każdy pojedynczy** dawał
ujemne \(e^2\) po skoku, po cichu przycinane do zera przez istniejący
w tym pliku strażnik `std::max(0.0, ...)`, kasujące prawdziwą informację
o mimośrodzie przy praktycznie każdym fotonie, po czym błąd propagował się
do każdej kolejnej oceny `kHere`.

Naprawione przez dokładne rozwiązanie relacji różniczkowej zamiast jej
ekstrapolacji: rozdzielenie zmiennych w \(s=1-e^2\) i \(x=\ln|E|\) daje
całkę pierwszą w postaci zamkniętej, \((1-e^2)^3/(e^4|E|^3)=\text{const}\),
dokładnie zachowaną wzdłuż dowolnej trajektorii spełniającej
\(dL/L=k(e)\,dE/E\) dla tego \(k(e)\). Zweryfikowane dwoma sposobami:
algebraiczna spójność z formułą mimośrodu, z której została wyprowadzona,
oraz niezależnie — przez bezpośrednie całkowanie tego równania różniczkowego
metodą RK4 (zgodność \(10^{-12}\), czyli błąd dyskretyzacji samego RK4, nie
formuły zamkniętej). Rozwiązywane dla każdego fotonu bisekcją na wynikowym
monotonicznym sześcianie w `crem_collapse.hpp` (wzór Cardana też by to
domknął, ale bisekcja omija casus irreducibilis — trzy rzeczywiste
pierwiastki wymagające wyboru gałęzi — a zdarzenia fotonowe są na tyle
rzadkie na trajektorię, że 80 kroków bisekcji nic nie kosztuje).

Statystyka zbiorcza pozostaje w przybliżeniu niezmieniona (10 ziaren,
o-Ps: mediana \(155{,}7\) ps wobec \(147{,}8\)–\(151{,}6\) ps zapisanych
wcześniej w tym dokumencie, 0 awarii numerycznych w obu przypadkach), mimo
że sama wartość mimośrodu przy niemal każdym fotonie zmieniała się o
wielkość rzędu jedności: czas kolapsu wyznacza głównie całka
energii/hazardu, której ta poprawka nie dotyka, a mimośród zasila
odległość peryapsis oraz warunek wyjścia po stosunku okres/czas
przelotu światła — i zarówno wartość dokładna, jak i stara (przycięta do
zera) były już małe (bliskie okręgu) przy każdym zmierzonym zdarzeniu w
tym przebiegu, więc różnią się o mniej niż \(0{,}01\) w bezwzględnym
\(e^2\), nawet tam gdzie stara wartość przed przycięciem była formalnie
niefizyczna (ujemna). `elements.specificAngularMomentum` jest przenoszone
między checkpointami, więc to nie jest poprawka jednorazowa — każda
kolejna ocena mimośrodu, `k` i warunku peryapsis po fotonie dziedziczy
teraz poprawną wartość zamiast błędnej, do końca tej trajektorii.
Zweryfikowane: czysta kompilacja (zero ostrzeżeń), `positronium_validation`
33/33, 0 awarii numerycznych na 10 ziaren.

**Podjęte na zadane polecenie: wstawiony spin fotonu, wcześniej opisany
jako zaburzenie rzędu jedności, którego "żadna klasyczna adiabatyczna
księgowość nie wchłonie" — i zmierzony jako przeżywalny.** Obie sekcje
wyżej (kierunek "niezmienny wobec odrzutu" i wielkość aktualizowana
współczynnikiem `k`) są teraz zastąpione jednym, prawdziwym kopnięciem
wektorowym: \(L_{para}\mathrel{-}=h\hbar\,\hat n_{foton}\) (\(h=\pm1\) to
skrętność fotonu), nakładanym na rzeczywisty (nie tylko specyficzny)
wektor momentu pędu pary i rozkładanym z powrotem na wielkość i kierunek.
Dlaczego to ZASTĘPUJE, a nie DODAJE do wyniku ze współczynnika `k`: `k`
pochodzi z klasycznego momentu siły reakcji, a z zachowania momentu pędu
zastosowanego do klasycznego pola "moment siły scałkowany po orbicie" i
"to, co unosi ciągłe pole" to ta sama wielkość, nie dwa oddzielne wkłady —
więc po skwantowaniu emisji w prawdziwy foton o znanym spinie, doliczanie
`k` na dodatek podwajałoby dokładnie to, co już policzone, tym samym
błędem, którym był usunięty wcześniej przechył dla pędu liniowego.

Skrętność \(h\) losowana jest z warunkowego (względem już wylosowanego
kąta emisji \(\theta\)) rozkładu standardowego dla przejścia dipolowego
\(\Delta m=\pm1\): \(P(h{=}{+}1|\theta)=(1+\cos\theta)^2/[2(1+\cos^2\theta)]\),
\(P(h{=}{-}1|\theta)=(1-\cos\theta)^2/[2(1+\cos^2\theta)]\)
— nie nowe założenie, te dwa prawdopodobieństwa sumują się dokładnie do
wzorca \((1+\cos^2\theta)\), z którego \(\theta\) było już losowane
(sprawdzone: przy \(\theta=0\) i \(\theta=\pi\) skrętność w połączeniu z
kierunkiem fotonu ZAWSZE daje wektor momentu pędu wzdłuż osi orbitalnej,
zgodnie z \(\Delta m=1\); przy \(\theta=\pi/2\) to równy rozkład 50/50 w
płaszczyźnie orbity). Czego ten mechanizm NIE łapie: orbitalnego momentu
pędu fotonu względem pary (potrzebna nieznana anomalia prawdziwa — ta sama
luka, która wykluczyła przechył). Sprawdzone wprost: wartość oczekiwana
składowej osiowej \(\langle h\cos\theta\rangle\) po całym rozkładzie
kątowym wynosi \(1/2\), nie \(1\) — sam spin odzyskuje więc średnio tylko
połowę reguły wyboru \(\Delta m=1\), reszta to ta sama niedostępna część
orbitalna.

Zmierzone, nie tylko wyargumentowane, że jest to przeżywalne: partia
produkcyjna (o-Ps, 80 trajektorii na dwóch ziarnach) dała 1 awarię
numeryczną (orbita czasem robi się, przez pojedyncze duże kopnięcie,
na tyle mimośrodowa, że opuszcza reżim związany — prawdziwa konsekwencja
zaburzenia rzędu jedności, nie błąd, i rzadsza niż mogłoby się wydawać) —
wobec 0 awarii, ale 3 ucięć budżetem czasowym dla tej samej wielkości
partii pod starym kodem (z zastosowanym `k`). Sam czas kolapsu ledwo
drgnął (mediana orto \(130\)–\(155\) ps w kilku ziarnach, wobec
\(147{,}8\)–\(151{,}6\) ps już zapisanych wcześniej), bo czas kolapsu
wyznacza całka energii/hazardu, której ta zmiana nie dotyka, a mimośród —
teraz naprawdę osiągający wartości jak \(e^2=0{,}9\), których stara ścieżka
z `k` nigdy nie produkowała — wpływa głównie na to, KTÓRY warunek wyjścia
trajektoria trafi i jak szybko, nie na to, czy w ogóle trafi.
Zweryfikowane: czysta kompilacja (zero ostrzeżeń), `positronium_validation`
33/33.

### Wynik audytu kompletności fizycznej

Model **nie jest dokładnym odwzorowaniem fizycznego układu elektron–pozyton**
ani zamkniętym, w pełni kowariantnym rozwiązaniem klasycznej elektrodynamiki
dwóch ładunków z samooddziaływaniem. Jest hybrydą trzech warstw:

1. zweryfikowanych składników klasycznych: relatywistycznej relacji
   pęd–prędkość, siły Lorentza, pól Liénarda–Wiecherta, tensora dipolowego,
   równań Maxwella na siatce oraz zachowania ładunku;
2. kontrolowanych przybliżeń: skończonego profilu Gaussa, regulatora
   krótkiego zasięgu, zredukowanej samosiły Landaua–Lifshitza, skończonej
   historii retardowanej i dyskretnych powierzchni kontrolnych;
3. modeli fenomenologicznych lub kwantowych dołączonych poza klasyczną
   dynamiką: klasyfikacji para/orto, czasów życia, zaniku 2γ/3γ i widma
   Ore'a–Powella.

Elektron i pozyton są w części orbitalnej traktowane symetrycznie: mają równe
masy, przeciwne ładunki oraz te same równania ruchu, przy czym znak ładunku
zmienia siłę Lorentza i precesję. Przypisanie każdej cząstce klasycznego,
trwałego momentu \(|\boldsymbol\mu|=(g/2)\mu_B\) przy mierzonym
\(g=2{,}00231930436256\) jest założeniem modelu, nie klasycznym wyprowadzeniem
własności elektronu. Moment jest wyprowadzany z masy, ładunku i \(g\) danego
gatunku, więc nie może rozjechać się z nimi; \(g\) nie jest przy tym bliskie
dwóm dla każdej cząstki — proton niesie \(g=5{,}5857\), bo jest złożony. Rzeczywisty
moment elektronu, spin 1/2, splątanie, wymiana, energia wiązania, anihilacja i
poprawki radiacyjne należą do QED i nie mogą zostać dokładnie odtworzone przez
ten integrator klasyczny.

Bieżący `positronium_validation` potwierdza własności algebraiczne i numeryczne,
ale jego końcowe `PASS` oznacza wyłącznie przejście ustawionych progów regresji.
Nie jest certyfikatem kompletności fizycznej. Audyt przeprowadzony ponownie
**19 sierpnia 2026 r.** dał następujący obraz.

Kowariancja i operatory pojedyncze (z `positronium_validation`):

- reszta kowariancji pojedynczego pola Liénarda–Wiecherta \(9{,}89\cdot10^{-16}\)
  i reszta siły po boostcie \(1{,}10\cdot10^{-6}\) — obie na poziomie
  numerycznym;
- reszta skończonego strumienia promieniowania po boostcie \(7{,}0\cdot10^{-4}\)
  i skumulowanego czteropędu promieniowania \(2{,}3\cdot10^{-6}\). Były to
  \(0{,}2084\) i \(0{,}2084\), zdecydowanie najgorsze liczby audytu; przyczyna
  została **znaleziona i naprawiona**, a droga do niej jest zapisana niżej, bo
  cztery z sześciu hipotez okazały się błędne.

  Przyczyną był **brakujący czynnik Dopplera** \(\kappa=1-\hat n\!\cdot\!
  \boldsymbol\beta\). Całka po sferze liczyła \(|E_{rad}|^2R^2\), czyli moc
  **odbieraną** na jednostkę czasu obserwatora, która niesie \(\kappa^{-6}\),
  podczas gdy moc **emitowana** na jednostkę czasu emitera niesie
  \(\kappa^{-5}\). Brakowało dokładnie jednego \(\kappa\). Dla źródła w ruchu
  dawało to błąd pierwszego rzędu w \(\beta\): przy \(\beta=0{,}35\) moc
  promieniowania wychodziła \(1{,}1955\) razy większa od spoczynkowej, czego
  zabrania niezmienniczość wzoru Larmora, a \(|p|/(E/c)\) wychodziło
  \(0{,}467\) tam, gdzie struktura czterowektorowa wymaga \(\beta=0{,}350\).
  Zgodność jest przy tym ścisła: \(1/(1-\beta\langle\cos\theta\rangle)
  =1/(1-0{,}35\cdot0{,}46712)=1{,}19544\) wobec zmierzonych \(1{,}1955\).

  Po przywróceniu \(\kappa\) reszta strumienia spada 298-krotnie,
  skumulowanego czteropędu 100 000-krotnie, a \(|p|/(E/c)\) wynosi
  \(0{,}35002\) wobec przewidywanych \(0{,}35003\). Źródła spoczywające są
  nietknięte, bo dla nich \(\kappa\equiv1\): stosunek do wzoru Larmora
  pozostaje \(1{,}00007\), a produkcyjna energia wypromieniowana jest bit
  w bit ta sama.

  Hipotezy obalone po drodze, zapisane, żeby ich nie powtarzać: **zbyt krótka
  historia opóźniona** (promień \(10^6a_0\) wymaga \(1{,}8\cdot10^{-13}\) s,
  historia daje \(7{,}2\cdot10^{-19}\) s, ale zejście do \(1\,a_0\) zmienia
  resztę tylko z \(0{,}2084\) na \(0{,}1695\)); **rozdzielczość kwadratury**
  (50 do 770 kierunków rusza piątą cyfrę); **transformacja przyspieszeń**
  (przyspieszenie własne, niezmiennik, zgadza się do \(2{,}6\cdot10^{-8}\));
  **niecentryczne próbkowanie powłoki** (wycentrowanie sfery na opóźnionej
  pozycji środka masy zmienia \(|p|/(E/c)\) z \(0{,}46712\) na \(0{,}46702\)).
  Rozstrzygnął dopiero skan promienia: stosunek \(|p|/(E/c)\) był stały co do
  czwartej cyfry przez pięć rzędów \(R\), więc nic się między powierzchniami
  nie gromadziło, strumień równał się tempu emisji i musiał spełniać
  \(|p|/(E/c)=\beta\);
- reszta ewolucji tensora dipolowego \(3{,}86\cdot10^{-3}\) — **nie jest to
  błąd całkowania**. Jest odporna na krok (ośmiokrotne zagęszczenie przy stałym
  czasie całkowitym nie zmienia jej ani na piątej cyfrze) i na tolerancję
  (identyczna od \(10^{-6}\) do \(10^{-12}\)), stała w czasie, a przy tym
  dokładnie pierwszego rzędu w prędkości orbitalnej: połowienie \(\beta_{orb}\)
  połowi ją, a wyzerowanie zbija ją 288-krotnie. Pokrywa się co do pięciu cyfr
  z raportowaną obok **luką reprezentacji**, bo jest tą samą wielkością
  widzianą z drugiej strony. Model przechowuje dipol **własny** i odtwarza
  tensor laboratoryjny jako \(\mathrm{boost}(\{0,\mathbf m\},\mathbf v)\), a ta
  parametryzacja nie wyraża każdego tensora: złożenie boostu układu z własną
  prędkością orbitalną cząstki to boost **razy obrót Wignera**, a obrócony
  tensor nie ma już tej postaci dla żadnego \(\mathbf m\). Ograniczenie jest
  więc strukturalne, rzędu \(\beta_{boost}\beta_{orb}\), i żadne zagęszczanie
  kroku go nie ruszy.

**Zbadane ponownie i potwierdzone rygorystycznie** (przy okazji audytu
fizyki, po naprawie precesji BMT — patrz niżej): rozłożone na dwie osobne
ścieżki i zweryfikowane numerycznie z osobna. Złożenie **dwóch kolejnych**
boostów samego tensora (`lorentzBoostDipole` wywołane dwa razy — raz
prędkością orbitalną, raz boostem układu) sprawdza się **dokładnie**
przeciw pojedynczemu boostowi złożoną (relatywistycznie) prędkością —
różnica to czyste zaokrąglenie maszynowe (`~10⁻¹⁸`), żadnego śladu obrotu.
To ma sens: tensor pola/dipola transformuje się jako
\(F'^{\mu\nu}=\Lambda^\mu_\alpha\Lambda^\nu_\beta F^{\alpha\beta}\) pod
**dowolnym** złożonym \(\Lambda\), więc dwa kolejne wywołania tej samej
formuły dają dokładnie poprawną odpowiedź niezależnie od tego, czy złożenie
dwóch nie-współliniowych boostów samo w sobie "zawiera" rotację. Rozbieżność
bierze się **wyłącznie** z pośredniego kroku przez czterowektor spinu
(`dipoleFourVector`→boost→`properDipoleFromFourVector`) — odtworzone
niezależnie z prawdziwymi stałymi testu, luka wyszła \(0{,}394\%\), tego
samego rzędu co zmierzone \(0{,}386\%\). Potwierdza to: reprezentacja
"pojedynczy wektor `properDipole` + jeden czysty boost" jest z definicji
niezdolna wyrazić każdy tensor osiągalny przez złożenie dwóch
nie-współliniowych boostów — to nie luka do domknięcia mniejszą poprawką,
tylko dowiedziona granica wybranej parametryzacji stanu.

**Sprawdzone też wprost, czy to ma jakikolwiek wpływ produkcyjny — nie ma.**
Operacja "przenieś cały stan symulacji do innej ramki obserwatora"
(`boostFourVector`/`boostEvent`) istnieje **wyłącznie** wewnątrz tego
jednego testu w `maxwell_validation.hpp`; żadna ścieżka produkcyjna nigdy
jej nie wykonuje (symulacja zawsze integruje w jednej, ustalonej ramce
laboratoryjnej). Prawdziwa naprawa wymagałaby rozszerzenia stanu o jawny
stopień swobody rotacji (kwaternion lub macierz 3×3 na cząstkę) plus osobnej
derywacji, jak ta rotacja narasta pod ciągłym, nie-współliniowym
przyspieszeniem — porównywalnej trudności do derywacji precesji BMT niżej —
i dotknęłaby rzędu 40 miejsc w kodzie (`properDipole`: 33, `DipoleTensor`/
`lorentzBoostDipole`: 10, poza samym plikiem walidacyjnym). Przy zerowym
wpływie na jakikolwiek zgłaszany wynik fizyczny — koszt tej wielkości nie
jest uzasadniony. Pozostawione jako świadomie zaakceptowane, teraz w pełni
udowodnione ograniczenie, nie zadanie do zrobienia.

**Sprawdzone osobno: czy w torze produkcyjnym jest gdzieś oddziaływanie
natychmiastowe (action-at-a-distance) zamiast retardowanego.** Prześledzony
bezpośrednio cały łańcuch wywołań składający się na siłę faktycznie pchającą
cząstki w `integrateElectrodynamicStep` (jedyny produkcyjny krok
integracji): pole ładunek-ładunek (`lienardWiechertField`), pole dipola
magnetycznego i elektrycznego (`retardedMagneticDipoleField`/
`retardedElectricDipoleField`), siła gradientowa na własnym dipole
(`covariantDipoleGradientForce`→`fieldFromOtherParticleAt`, te same trzy
pola retardowane) i reakcja promieniowania (różnica skończona z
`retardedExternalForces`, nie `allExternalForces`, w stanach `before`/
`after`) — wszystkie retardowane, przez `history`. Jedyne nieretardowane
składniki w torze produkcyjnym to jednorodne pole zewnętrzne
(`gExternalMagneticField`) i pole punktu zerowego (`gZeroPointField`) — oba
jawnie udokumentowane jako pola **tła**, nie źródłowane przez żadną z
cząstek, więc retardacja nie ma tu zastosowania fizycznego (nie ma od kogo
liczyć opóźnienia).

Jedno miejsce wymagało bliższego sprawdzenia: `chargeDipoleForces` liczy
siłę wprost z bieżącej, nie retardowanej separacji. Ale ta funkcja ma
dokładnie jedno wywołanie w całym kodzie — wewnątrz `allExternalForces`
(natychmiastowej sumy sił) — a `allExternalForces` ma z kolei dokładnie
trzy wywołania, wszystkie poza dynamiką: `causalInitialHistory`
(`crem_engine.hpp:50`, jednorazowy zalążek historii przed pierwszym
prawdziwym krokiem, natychmiast poprawiany iteracją Picarda względem wzoru
retardowanego), `integrateConservativeMidpoint`
(`electrodynamics.hpp:2147`, test odwracalności wyłącznie pod
`POSITRONIUM_ENABLE_FIELD_VALIDATION`) i `makeFrame`
(`crem_trajectory.hpp:28`, czysto diagnostyczna energia Schotta do raportu,
nie napędza ruchu).

Wniosek: żadna siła, która faktycznie zmienia pęd cząstki w symulacji
produkcyjnej, nie jest natychmiastowa. To weryfikacja, nie znalezisko —
kod już wcześniej twierdził to samo w kilku miejscach, ale teraz
prześledzone bezpośrednio przez cały łańcuch wywołań, nie tylko
zacytowane.

**Sprawdzone też: czy wszędzie są używane rzeczywiście zmierzone momenty
magnetyczne cząstek.** Architektura już wcześniej broni się przed
niespójnością tego typu — moment magnetyczny jest wyłącznie **wyprowadzany**
z g-faktora (`magneticMoment(species) = 0,5·g·magneton(species)`), nigdy
tabelaryzowany osobno; komentarz przy tym w `particle_species.hpp` wprost
opisuje wcześniejszy incydent ("g-factor of 2.0023 sitting beside a moment
of mu_B"). Prześledzone: `firstMagneticMoment`/`secondMagneticMoment` i
`firstGFactor`/`secondGFactor` mają dokładnie po dwa miejsca przypisania w
całym kodzie (inicjalizacja statyczna i `applyPairFromOption` przy zmianie
`--pair`), oba przez tę samą funkcję, żadnej równoległej, potencjalnie
rozjeżdżającej się wartości nigdzie indziej — zero trafień na literalny
`g=2,0` czy twardo wpisaną wartość magnetonu poza `physical_constants.hpp`.

Same liczby zweryfikowane bezpośrednio względem physics.nist.gov (CODATA
2022): masy elektronu, mionu i protonu — **zgodne co do wszystkich
cytowanych cyfr**. G-faktory mionu (\(2{,}00233184123\)) i protonu
(\(5{,}5856946893\)) — **zgodne dokładnie**. Bohr magneton — zgodny.
Jedno realne znalezisko: g-faktor elektronu w kodzie (`2,00231930436256`)
to wartość CODATA **2018**, nie 2022, jak twierdził komentarz w
`particle_species.hpp` ("g-factors from the same source"). Aktualna CODATA
2022 to `2,00231930436092` — przesunięcie o \(8{,}2\cdot10^{-13}\)
względnie (11.-12. cyfra znacząca), z aktualizacji uwzględniającej pomiar
Fan/Gabrielse 2023. Naprawione (`physical_constants.hpp`,
`particle_species.hpp`): stała zaktualizowana, cytat w komentarzu
poprawiony. Efekt na jakikolwiek zgłaszany wynik: **żaden** — przesunięcie
jest o rzędy wielkości poniżej najciaśniejszej tolerancji, jaką ten projekt
gdziekolwiek sprawdza (\(10^{-8}\) i wyżej). Zweryfikowane:
`positronium_validation` 33/33 bez zmian po poprawce.

**Sprawdzone też: czy warunki startowe pięciu eksperymentów statystycznych
odpowiadają praktyce eksperymentalnej.** Losowanie kąta uderzenia w
przekroju czynnym (eksperymenty 3, 4) i próbki parametru zderzenia po
`√(losowa liczba)` (jednorodnie po **polu**, nie po promieniu) to poprawna
technika Monte Carlo dla ekstrakcji przekroju czynnego niezależnie od
rzeczywistego profilu wiązki — nie jest to uproszczenie do naprawienia.
Eksperyment 5 już wcześniej próbkuje energię kinetyczną obu cząstek z
rozkładu Gaussa (`sampleKinetic`) i parametr zderzenia z rozkładu
półnormalnego \(|\mathcal N(0,\sigma)|\), odzwierciedlając skończoną
rozdzielczość energetyczną i ogniskowanie realnej wiązki. Eksperymenty 1/2
(kolaps CREM) startują z ustalonego \(L=\hbar\) na promieniu Bohra pary —
to nie próbka eksperymentalna, tylko warunek początkowy modelu klasycznego
odtwarzającego zadaną energię wiązania (patrz komentarz przy
`positroniumBohrRadius` w `physical_constants.hpp`), więc rozkład tu nie
ma odpowiednika fizycznego do którego by dążyć.

Jedno realne odstępstwo: eksperymenty 3 i 4 uderzają wiązką o **ustalonej**
energii środka masy \(K_{CM}\) — zamierzone, bo ich celem jest porównanie
zmierzonego przekroju czynnego z formułą Rutherforda przy jednej, znanej
energii (sama \(\sigma(\theta)\) jest zdefiniowana przy ustalonym
\(K_{CM}\)), ale to jednocześnie odbiega od realnej wiązki, która zawsze ma
skończoną rozdzielczość energetyczną. Dodana opcjonalna poprawka: nowa
flaga `--beam-energy-sigma-ev` (domyślnie `0`, zachowując dotychczasowe
zachowanie bit w bit), która pozwala próbkować \(K_{CM}\) danego zdarzenia
z rozkładu Gaussa o zadanym \(\sigma\) wokół `--beam-energy-ev`, tym samym
próbkowaniem odrzucającym co `sampleKinetic` w eksperymencie 5 (odrzuca
próbki \(\le 0\), do 1000 prób, po czym pada z powrotem na wartość
średnią). Nowe pole `BeamEvent::incidentEnergy` niesie faktycznie
wylosowaną energię danego zdarzenia; widmo strat energii
(`fiducialEnergyLossesEv`, jedyne miejsce liczące \(\Delta E\) wprost)
poprawione, by liczyć różnicę względem tej rzeczywistej energii zdarzenia,
nie względem ustalonej średniej konfiguracji — inaczej przy \(\sigma>0\)
każde zdarzenie dostawałoby błędny, systematyczny offset. Baner startowy
wypisuje \(\sigma\) obok \(K_{CM}\), gdy jest niezerowa.

Zweryfikowane: kompilacja czysta (zero ostrzeżeń), `positronium_validation`
33/33 bez zmian. Test dymny `--phenomenon 3 --beam-energy-ev 20
--beam-energy-sigma-ev 5 --runs 8` przechodzi bez awarii, baner poprawnie
pokazuje `K_CM = 20 eV (Gaussian sigma = 5 eV)`, widmo strat energii i
przekrój czynny liczą się dalej. Przy \(\sigma=0\) (domyślnie) próbkowanie
Gaussa jest pomijane przez krótkie spięcie przed jakimkolwiek losowaniem
(`if (!(sigma > 0.0)) return configuration.centreOfMassKineticEnergy;`),
więc strumień losowy pozostaje nietknięty — zachowanie domyślne jest
identyczne co przed tą zmianą.

**Dodana referencja: amplitudy Bhabhy/Motta obok Rutherforda.** Wykres
różniczkowego przekroju czynnego (eksperymenty 3, 4) od zawsze zestawiał
model trajektorii z klasyczną formułą Rutherforda i osobno zastrzegał w
tekście, że dane są "classical low-velocity model, not Bhabha scattering" —
prawdziwe pytanie, jak daleko klasyczny model faktycznie odbiega od pełnej
relatywistycznej QED, pozostawało otwarte. Dodana druga, niezależna krzywa
referencyjna: pełna, drzewowa (jednopętlowa nie, tylko diagramy na poziomie
drzewa) amplituda QED, obliczona osobno od trajektorii klasycznej i
narysowana obok Rutherforda, nie zamiast niego.

Dwa różne procesy, zależnie od pary. Dla dwóch RÓŻNYCH gatunków (np.
e⁻+p̄) istnieje tylko jeden diagram: wymiana fotonu w kanale t między
dwiema nietożsamymi liniami fermionowymi (uogólnione rozpraszanie Motta) —
wierzchołek foton-fermion łączy linię fermionową wyłącznie samą ze sobą,
więc kanał anihilacyjny nie istnieje dla par mieszanych. Dla pary
cząstka-własna antycząstka (e⁺e⁻, μ⁺μ⁻, p+p̄ — jedyne takie pary, jakie
program obsługuje, rozpoznawane po `mass1==mass2`, dokładnie warunku, przy
którym anihilacja do wirtualnego fotonu i z powrotem do tych samych dwóch
cząstek w stanie końcowym jest w ogóle możliwa) dochodzi drugi diagram —
anihilacja w kanale s — i oba diagramy interferują, bo kończą się tym samym
stanem końcowym. e⁺e⁻→e⁺e⁻ konkretnie to właśnie nazwa "rozpraszanie
Bhabhy"; ta sama struktura dwóch diagramów przenosi się bez zmian na
μ⁺μ⁻ i p+p̄.

Obie formuły (kwadrat amplitudy, uśredniony po spinach) wyprowadzone ręcznie
ze standardowej techniki śladów Diraca (Σuū=p̸+m, Σvv̄=p̸−m) i zweryfikowane
niezależnie: jawny numeryczny rachunek spinorowy (konkretne macierze Diraca,
konkretne spinory u/v w bazie helicity, siłowa suma po wszystkich 16
kombinacjach spinów) dla szerokiego skanu mas, energii i kątów — zgodność
na poziomie \(10^{-9}\) względnie lub lepiej wszędzie. Człon interferencyjny
to dokładnie ten fragment, który podręczniki zwykle podają bez wyprowadzenia
— pierwsza ręczna próba była **błędna** (brakujący czynnik podniesiony przy
kontrakcji przez slash pędu), a to właśnie weryfikacja numeryczna złapała
błąd, zanim trafił do kodu — dokładnie ta sama dyscyplina audytu, co reszta
wyprowadzeń fizycznych tego projektu. Obie formuły w granicy bezmasowej
redukują się do podręcznikowych wyników (rozpraszanie e-μ:
\(|M|^2=2e^4(s^2+u^2)/t^2\); e⁺e⁻→μ⁺μ⁻: \(|M|^2=2e^4(t^2+u^2)/s^2\);
interferencja Bhabhy: \(-2\,\mathrm{Re}(M_tM_s^*)=4e^4u^2/(st)\); źródła:
Peskin i Schroeder, *An Introduction to Quantum Field Theory*, wzór 5.61 i
wynik e⁺e⁻→μ⁺μ⁻ z rozdziału 5.1; samo rozpraszanie Bhabhy: H.J. Bhabha,
*Proc. R. Soc. A* 154, 195 (1936)), a przy \(K_{CM}=20\) eV (energia
produkcyjna wiązki) obie zgadzają się z już istniejącą, klasyczną formułą
Rutherforda tego kodu na poziomie \(10^{-4}\) — dokładnie rzędu poprawki
\(\beta^2\sim4\cdot10^{-5}\) spodziewanej przy tej energii, sprawdzone
zarówno dla pary o równych masach (e⁺e⁻), jak i różnych (e⁻+proton).

Implementacja: nowy moduł `modules/qed_reference.hpp`
(`qedElasticDifferentialCrossSection`, z dyspozytorem wybierającym Bhabhę
albo Motta po `mass1==mass2`), zero zależności od reszty silnika poza
`physical_constants.hpp` — czysto referencyjny, nigdzie nie wpływa na
całkowanie trajektorii. Krzywa liczona metodą Simpsona (9 punktów na
koszyk kątowy w `cos(theta)`, bo w przeciwieństwie do Rutherforda QED nie
ma prostej analitycznej całki koszykowej) i rysowana na panelu 1
(`3_1_1`/`4_1_1_differential_cross_section.pdf`) trzecim kolorem palety
Okabe-Ito (`#CC79A7`, czerwonawy fiolet — jedyny z siedmiu jeszcze
nieużyty), z legendą podpisującą się "tree-level Bhabha" albo "tree-level
Mott" zależnie od pary. Zastrzeżenie w tekście konsoli poprawione, by nie
przeczyć nowej krzywej: dane trajektorii nadal są klasyczne, ale QED jest
teraz obecna na wykresie jako referencja, nie tylko przywoływana w zdaniu
zaprzeczającym.

Zweryfikowane: kompilacja czysta, `positronium_validation` 33/33 bez zmian.
Testy dymne `--phenomenon 4 --pair electron,positron` (gałąź Bhabhy) i
`--phenomenon 4 --pair electron,proton` (gałąź Motta, różne masy) oraz
`--phenomenon 3` przechodzą bez awarii; wykresy sprawdzone wizualnie —
krzywa QED pokrywa się z Rutherfordem na tej skali, zgodnie z oczekiwaniem
przy \(\beta\sim0{,}006\), a legenda poprawnie przełącza etykietę między
parami.

Relacja moment–spin niosła błąd czynnika \(g\). Model definiuje
\(\boldsymbol\mu=\gamma\mathbf S\) z \(\gamma=gq/2m\) i przechowuje
\(|\boldsymbol\mu|=(g/2)\,\text{magneton}\), ale cztery miejsca liczyły
\(\gamma\) jako \(q/2m\): wewnętrzny moment pędu w `particleFieldTotals`
i w `noetherAngularMomentum` (obie postaci \(\boldsymbol\mu/\gamma\)) oraz
odpowiedź na moment siły \(d\boldsymbol\mu=\gamma\boldsymbol\tau\,dt\)
w `applyDipoleRadiationTorque` i w `pushStateWithGridField`. Raportowany
wewnętrzny moment pędu wychodził przez to \(1{,}001160\,\hbar\) zamiast
\(0{,}500000\,\hbar\), czyli był zawyżony dokładnie o \(g=2{,}0023\),
a odpowiedź dipola na moment siły była o tyle samo za słaba. Wszystkie cztery
miejsca korzystają teraz ze wspólnego `gyromagneticRatio(charge, mass, gFactor)`,
żeby \(g\) nie mogło wypaść przy kolejnym wywołaniu. Niezależne residuum
`raw dJ` poprawiło się o 3,7%.

Podwojenie `particle-field dJ` po tej poprawce w relacji
\(\boldsymbol\mu=\gamma\mathbf S\) zostało **zbadane i jest artefaktem
normalizacji, nie regresją**. Licznik domknięcia nie zmienił się ani o bit
i wynosi \(9{,}57308412\cdot10^{-38}\) przed poprawką i po niej. Zmienił się
mianownik: skala \(|J_{\rm pocz}|\) zawiera wkład spinowy, czyli dokładnie ten
człon, który był zawyżony o \(g\), więc spadła z \(2{,}11399975\cdot10^{-34}\)
na \(1{,}05698206\cdot10^{-34}\) — z \(2{,}0046\,\hbar\) na
\(1{,}0023\,\hbar\). Iloraz podwoił się z \(4{,}53\cdot10^{-4}\) na
\(9{,}06\cdot10^{-4}\) bez żadnej zmiany fizycznej, a dawna wartość była po
prostu schlebiana zawyżonym mianownikiem.

Próg podniesiono z \(10^{-3}\) na \(2\cdot10^{-3}\), żeby kontrola zachowała
dotychczasową surowość **bezwzględną**; pozostawienie go uczyniłoby ją
dwukrotnie ostrzejszą jako niezamierzony skutek uboczny poprawki jednostek.

Bilans energii jest **tożsamością algebraiczną**, nie testem zachowania.
Dotyczy to obu miejsc, w których się pojawia: pozycji `identity resid` /
`identity |dP|` / `identity |dJ|` w trybie diagnostycznym oraz linii
`pair-field IDENTITY` w walidacji. Ta druga była do niedawna opisana jako
`pair-field dE/dP/dJ` i dawała `0 / 0 / 1,1e-26`, co czytało się jak ścisłe
zachowanie energii, pędu i momentu pędu; jest to jednak suma teleskopująca,
bo `boundField*` jest akumulowane co krok właśnie jako reszta domykająca.
Sprawdzone pomiarem: całkowite złamanie antysymetrii sondy nie rusza tych
liczb wcale.

Zestaw kontroli rozszerzono o pięć pozycji towarzyszących modułowi
`two_body_kinematics.hpp`: `two-body-role-invariance` sprawdza, że zamiana ról
w parze nie zmienia kinematyki wejściowej ani swobodnej,
`two-body-lorentz-boost` — że boost i powrót odtwarzają oba stany,
`two-body-causality` — przyczynowość ustawienia dwuciałowego,
`coherent-magnetic-dipole` — residuum pochodnej w koherentnym sektorze
magnetycznym, a `adaptive-depth-rejection` — zachowanie integratora przy
odrzucaniu kroku na granicy głębokości podziału.

Wzorzec „wielkość wypisywana, ale niesprawdzana" został przejrzany
**systematycznie**: dla każdej ze 135 wielkości w wyjściu walidacji sprawdzono,
czy jej nazwa występuje w jakimkolwiek predykacie. Trzynaście nie występowało.
Dziesięć z nich jest niebramkowanych słusznie — pięć to pomiary wydajności
(bramkowanie ich kazałoby suite wywalać się na obciążonej maszynie), cztery to
obiekty pomocnicze i pętla raportująca, a `gaussBefore` i `quadrupoleMagnitude`
są kontekstem, którego znormalizowane odpowiedniki mają własne bramki.
`cflStep` jest egzekwowany strukturalnie, bo każdy krok w suite jest ułamkiem
`courantTimeStep()`. Pozostałe trzy — trzy wejścia bramki mieszającej modele
reakcji — dostały progi, a razem z nimi `llValidity`, sprawdzany dotąd wyłącznie
przez `isfinite`.

Jedna z tych bramek jest nietypowa i warto wiedzieć dlaczego. `smooth`, czyli
zgodność oszacowań trzeciej pochodnej momentu dipolowego przy dwóch krokach,
sięga \(2{,}611\) dla mionium — tam ta pochodna jest po prostu szumem. Progu
na samą wartość postawić się więc nie da; sprawdzane jest zamiast tego, czy
**bramka to zauważyła i wyzerowała wagę**. Przez cztery pary:
\(3{,}97\cdot10^{-4}\to1\), \(1{,}66\cdot10^{-3}\to1\),
\(2{,}611\to0\), \(2{,}79\cdot10^{-2}\to0{,}552\).

Ten sam wzorzec dotyczył całego sektora reakcji promieniowania. `finiteReactionBenchmark` testował wyłącznie
`isfinite`, więc dziewięć liczb opisujących trzy modele reakcji mogło przyjąć
dowolną wartość, a suite i tak raportowałby komplet zdanych kontroli. Mają teraz progi.

Przy okazji ujawniła się pułapka normalizacji. `reaction flux off/LL/C` dzieli
rozjazd pracy reakcji przez **energię wypromieniowaną**, a ta w sondzie trwającej
\(8\cdot10^{-22}\) s wynosi dla pary protonów \(1{,}2\cdot10^{-29}\) J. Iloraz
sięga wtedy \(1893\) przy \(0{,}026\) dla e⁺e⁻, co czyta się jak katastrofa,
a jest wyłącznie znikającym mianownikiem. Ta sama wielkość odniesiona do energii
mechanicznej orbity wynosi \(9{,}2\cdot10^{-7}\) dla e⁺e⁻,
\(5{,}6\cdot10^{-11}\) dla mionium i \(7{,}1\cdot10^{-13}\) dla protonium —
znikoma i **malejąca** z masą, czyli wniosek dokładnie przeciwny. Bramkowana
jest ta druga postać, dopisana jako `reaction mism/E_mech`.

Niezależnym pomiarem jest linia `raw dE/dP/dJ`, która wyklucza `boundField*`.
Do niedawna **była wypisywana i nigdy nie sprawdzana**, więc mogła dryfować
dowolnie przy komplecie zdanych kontroli; ma teraz własne progi. Jej składowa pędowa
wymaga przy tym ostrożności tego samego rodzaju co reszta bilansu pola: wysyłana
sonda jest antysymetryczna, przez co pęd wychodzi 1,3e-08, czyli pięć rzędów
lepiej niż energia — po złamaniu symetrii wychodzi 8,1e-04, czyli tyle samo co
energia (8,4e-04) i moment pędu (1,8e-04).


`boundField*` liczone jest jako brakująca reszta bilansu po każdym kroku, więc
wykresy `diagnostic_*_balance` i pozycje `identity resid` / `identity |dP|` /
`identity |dJ|` mierzą wyłącznie zaokrąglenia i błąd interpolacji punktu
końcowego. Niezależne są trzy miary poniżej, zmierzone dla trzech ziaren na
kanał w trybie diagnostycznym. `E_rel` normalizuje rezerwuar pola związanego do
skali energii samej orbity, co jest uczciwsze niż dzielenie go przez energię
wypromieniowaną, która bywa dowolnie mała:

| kanał | \(\lvert E_{bound}\rvert/\lvert E_{rel}\rvert\) | \(\lvert E_{Schott}\rvert/E_{rad}\) | \(\lvert dE_{LL\text{-}flux}\rvert/E_{rad}\) |
| --- | --- | --- | --- |
| orbity związane (1, 2) | \(3{,}6\cdot10^{-5}\)–\(1{,}4\cdot10^{-4}\) | \(7\cdot10^{-5}\)–\(2{,}9\cdot10^{-3}\) | \(5{,}1\cdot10^{-4}\)–\(8{,}2\cdot10^{-3}\) |
| rozpraszanie szerokie (4) | \(2{,}1\cdot10^{-4}\)–\(2{,}5\cdot10^{-4}\) | \(3{,}3\cdot10^{-6}\)–\(4{,}9\cdot10^{-5}\) | \(0{,}026\)–\(0{,}28\) |
| **kanał krótkiego zasięgu (3)** | **1,4–9,4** | **3,4–3,7** | **6,2–6,4** |

Wniosek jest ostrzejszy niż w poprzednim audycie, który podawał jeden zakres
\(0{,}03\)–\(0{,}09\) i nie ujawniał rozwarstwienia. Kanały związany i szeroko
rozpraszający są czyste we wszystkich trzech miarach: rezerwuar pola bliskiego
to \(10^{-4}\) energii orbity, człon Schotta ułamek promila energii
wypromieniowanej, a praca reakcji LL zgadza się ze strumieniem dalekiego pola
do promila. **Kanał krótkiego zasięgu jest natomiast poza zakresem
stosowalności modelu wedle jego własnej diagnostyki**: rezerwuar pola
związanego przekracza energię orbity nawet dziewięciokrotnie, człon Schotta jest
trzy i pół raza większy od energii wypromieniowanej, a praca LL rozmija się ze
strumieniem o czynnik sześć. To nie jest efekt małego mianownika — bezwzględnie
przy \(E_{rel}=-7{,}8\) eV wychodzi \(E_{bound}=+11{,}0\) eV i
\(E_{Schott}=-1{,}5\) eV wobec \(0{,}45\) eV wypromieniowanych. Liczb
z eksperymentu 3 nie wolno traktować jako pomiaru energetycznego.

Względem audytu z 14 sierpnia zmieniło się przy tym pięć rzeczy: dołożono człon
interferencyjny w indywidualnej reakcji LL (residuum \(0{,}501\to0{,}003\),
czas kolapsu skrócony dokładnie dwukrotnie), przeniesiono początek układu
kwadrupola elektrycznego do środka masy, wyprowadzono promień regularyzacji
dipola z wybranej pary, dodano wybór pary w czasie działania oraz niezależne
losowanie energii obu cząstek w eksperymencie 5.

W konsekwencji ilościowo wiarygodne są testy jednostkowe danego operatora w
podanym zakresie rozdzielczości. Długoczasowe trajektorie, promieniowanie,
reakcja promieniowania, klasyfikacja para/orto i zachowanie przy cutoffie są
wynikami przyjętego modelu i wymagają osobnej analizy zbieżności oraz
wrażliwości. Bez niezależnej całki tensora energii–pędu na wspólnej
hiperpowierzchni, testu kowariancji dążącego do zera z krokiem i
rozdzielczością oraz porównania z QED nie wolno interpretować ich jako
precyzyjnych przewidywań pozytonium.

### Ulepszenia dokładności o małym koszcie

Całka dalekiego pola używa obecnie domyślnie 50-punktowej kwadratury Lebiediewa
stopnia 11. Reguła zachowuje symetrię oktagonalną sfery i zastępuje dla
produkcyjnych 50 kierunków siatkę Fibonacciego; ta ostatnia pozostaje tylko dla
niestandardowych liczebności w teście zbieżności. Reszta sumy wag wynosi
\(8{,}48\cdot10^{-16}\), a reszta pierwszego i drugiego momentu sferycznego
\(2{,}12\cdot10^{-16}\). Błąd względem gęstej referencji spadł z
\(5{,}51\cdot10^{-3}\) do \(8{,}59\cdot10^{-4}\) bez zwiększenia liczby
próbek.

Walidator raportuje teraz dwa bilanse. `pair-field dE/dP/dJ` nadal obejmuje
rekonstruowany rezerwuar `boundField*`, natomiast `raw dE/dP/dJ` całkowicie go
pomija i dlatego jest niezależnym testem niedomknięcia. Dla krótkiego testu
wartości raw wynoszą odpowiednio około \(3{,}74\cdot10^{-4}\),
\(1{,}34\cdot10^{-8}\) i \(3{,}00\cdot10^{-4}\).

Historia początkowa przechodzi dwie iteracje Picarda. Pierwsza parabola jest
budowana z siły chwilowej, po czym przyspieszenie jest dwukrotnie przeliczane
z pełnego oddziaływania retardowanego i historia jest odbudowywana. Ogranicza
to niedopasowanie przy \(t=0\) bez wydłużania właściwej trajektorii. Kroki
różnicowe pochodnych dipolowych są związane z lokalnym odstępem historii przez
współczynnik 2 zamiast wcześniejszego 8, co zmniejsza błąd obcięcia bez
rekurencyjnego obliczania pola BMT. Sam odstęp historii nie idzie już jednak
za krokiem całkowania — patrz „Siatka historii retardowanej".

Dostępny jest alternatywny model `coherentElectricDipole`: wspólne pole reakcji
Abrahama–Lorentza jest wyznaczane z tej samej trzeciej pochodnej elektrycznego
momentu dipolowego pary. Kanał M1 jest zawsze liczony koherentnie z
\(\mathbf m=\mathbf m_1+\mathbf m_2\), zatem zachowuje zarówno konstruktywną,
jak i destruktywną interferencję obu dipoli. Siły na elektron i pozyton są
dokładnie przeciwne;
test wielomianu sześciennego daje względny błąd pochodnej
\(8{,}72\cdot10^{-15}\). Model pozostaje opcjonalny, ponieważ opisuje wprost
tylko dominujący kanał E1; domyślna trajektoria nadal używa indywidualnego LL.
Walidator podaje także wskaźnik \(|F_{LL}|/|F_{ext}|\), obecnie około
\(4{,}00\cdot10^{-5}\), aby kontrolować założenie małej redukcji rzędu.

Walidator uruchamia identyczny krótki przebieg bez reakcji, z LL oraz z
koherentnym E1+M1. Dla kolejnych wariantów niezależna reszta energii `raw`
wynosi \(3{,}68\cdot10^{-6}\), \(4{,}11\cdot10^{-6}\) i
\(4{,}30\cdot10^{-6}\), a niezgodność pracy reakcji ze strumieniem odpowiednio
1, \(0{,}487\) i \(0{,}299\). Koszt koherentnego wariantu jest zbliżony do LL,
lecz jego różnica między krokiem \(\Delta t\) i \(\Delta t/2\) wynosi
\(8{,}40\cdot10^{-6}\), wobec \(2{,}29\cdot10^{-11}\) dla LL. Nie spełnia więc
przyjętego kryterium 1% ani kryterium zbieżności i nie został ustawiony jako
domyślny. Eksperymentalny selektor `automatic` sprawdza zwartość źródła,
dominację E1, gładkość \(\dddot{\mathbf p}\) i małość LL. Zamiast progu
zero–jeden wyznacza ciągłą wagę czterech funkcji przejścia i miesza obie siły.
W krótkim benchmarku zachowuje wynik LL: `raw=4,11e-6`, niezgodność
reakcja–strumień `0,487` i resztę kroku `2,29e-11`. Nie jest domyślny, ponieważ
nie wykazał jeszcze poprawy względem LL, a czysty wariant koherentny nie
spełnił kryterium 1%.

Do testów części zachowawczej dodano samouzgodnioną mapę niejawnego punktu
środkowego, działającą również dla ujemnego kroku. Pełny test przód–tył bez
promieniowania i historii retardowanej daje resztę \(2{,}75\cdot10^{-19}\).
Nie zmienia to faktu, że produkcyjny model przyczynowy z promieniowaniem jest
fizycznie nieodwracalny.

Rozwijany backend pola Maxwella używa gładkiego profilu Gaussa, sferycznego w
chwilowym układzie spoczynkowym cząstki, o promieniu regularizacji
\(r_c=0{,}01a_0\). W układzie laboratoryjnym profil jest skrócony Lorentzowsko,
a jego gęstość pomnożona przez \(\gamma\), dzięki czemu całkowity ładunek
pozostaje równy \(q\); gęstość prądu wynosi \(\mathbf J=\rho\mathbf v\).
Parametr \(r_c\) jest odcięciem numerycznym planowanej siatki adaptacyjnej, a
nie fizycznym promieniem elektronu.

Profil jest teraz opisany jednym modelem kowariantnej materii rozszerzonej.
W chwilowym układzie spoczynkowym definiowane są gęstość właściwa, tensor
polaryzacji–magnetyzacji \(M^{\mu\nu}\) oraz symetryczny tensor materii
\(T_{\rm matter}^{\mu\nu}\); wszystkie trzy są następnie transformowane tą
samą macierzą Lorentza. Swobodny czteroprąd ma postać

\[
J^\mu_{\rm free}=qf(\mathbf x)u^\mu,
\]

a przestrzenne składowe \(M^{ij}\) zawierają magnetyzację. Boost automatycznie
generuje składowe \(M^{0i}\), czyli polaryzację elektryczną poruszającego się
dipola. Tensor materii zawiera również izotropowe naprężenie kohezyjne. Jest
ono wyznaczane z równowagi statycznej profilu Gaussa przez całkę radialnej
gęstości siły \(|\rho E|\), nie dobierane jako swobodny parametr.

Test algebraiczny daje względny dryf normy czteroprądu
\(1{,}96\cdot10^{-16}\), dokładną antysymetrię \(M^{\mu\nu}\) oraz dokładną
symetrię \(T_{\rm matter}^{\mu\nu}\). Energia własna pola nie została jeszcze
pominięta w dynamice: poniżej opisano jej jawne rozdzielenie między masę gołą
i elektromagnetyczną oraz numeryczną kalibrację samosiły.

Dla użytego profilu Gaussa analityczna energia własna wynosi

\[
U_{\rm self}=\frac{q^2}{8\pi^{3/2}\varepsilon_0r_c}.
\]

Masa fizyczna jest rozdzielana na

\[
m_{\rm phys}=m_{\rm bare}+U_{\rm self}/c^2.
\]

Dla \(r_c=0{,}01a_0\) część elektromagnetyczna stanowi
\(m_{\rm em}/m_{\rm phys}=0{,}00150219\), więc
\(m_{\rm bare}/m_{\rm phys}=0{,}99849781>0\). Program sprawdza dodatniość
masy gołej; zmniejszenie \(r_c\) może ostatecznie naruszyć ten warunek i wtedy
model klasycznej chmury przestaje być dopuszczalny bez dodatkowej
renormalizacji tensora materii.

Pole własne dyskretnej chmury jest kalibrowane na izolowanym bloku o tej samej
geometrii, kroku siatki, prędkości i fazie podkomórkowej, a następnie odejmowane
od pola używanego przez pusher. Pole drugiej cząstki nie jest odejmowane.
Test ośmiu przesunięć podkomórkowych wykazał, że niepoprawiona samosiła może
osiągać około \(1{,}09\%\) referencyjnej siły wzajemnej przy odległości
\(2{,}5r_c\), więc korekta nie jest zaniedbywalna. Krótki test sprzężony używa
teraz dynamicznej tablicy samopola. Jednorazowa kalibracja zawiera \(4^3=64\)
próbki fazy komórki w chwilowym układzie spoczynkowym; periodyczna interpolacja
trójliniowa usuwa skok przy przejściu przez ścianę komórki. Pole jest następnie
transformowane tensorem elektromagnetycznym do aktualnej prędkości i skalowane
znakiem ładunku. Pusher pobiera nową korektę w każdym półkroku, więc nie używa
już pola własnego zamrożonego w stanie początkowym. Względna nieciągłość na
szwie fazy wynosi \(1{,}38\cdot10^{-7}\), a test zmiany znaku ładunku daje
resztę równą zero.

Backend zawiera obecnie trójwymiarową hierarchię trzech zagnieżdżonych bloków
Maxwella o stosunku zagęszczenia 2. Najdrobniejszy krok przestrzenny wynosi
\(\Delta x=r_c/2\). Każdy blok przechowuje wszystkie składowe
\(\mathbf E\), \(\mathbf B\), \(\rho\) i \(\mathbf J\), a krok czasu spełnia
trójwymiarowy warunek CFL. Równania rotacji są całkowane metodą punktu
środkowego. Centralne operatory różnicowe spełniają dyskretnie
\(\nabla_h\!\cdot(\nabla_h\!\times)=0\), natomiast ograniczenie elektryczne

\[
\nabla_h\cdot\mathbf E=\rho/\varepsilon_0
\]

jest przywracane zgodną projekcją Poissona. Nanoszenie każdego obciętego
profilu Gaussa jest renormalizowane, aby jego dyskretna całka była dokładnie
równa ładunkowi cząstki.

Rozpoczęto migrację właściwego solvera na siatkę Yee. Nowy rdzeń przechowuje
normalne składowe \(\mathbf E\) i \(\mathbf J\) na ścianach, składowe
\(\mathbf B\) na komplementarnych ścianach oraz \(\rho\) w środkach komórek.
Równania Faradaya i Ampère’a używają par zgodnych różnic naprzód/wstecz, dzięki
czemu dyskretna dywergencja rotacji znika algebraicznie.

Prąd swobodny CIC jest nanoszony lokalnym operatorem dekompozycji trajektorii
typu Esirkepova. Droga w jednym kroku jest dzielona symetrycznie na odcinki
jednoosiowe i dalej na fragmenty mieszczące się w pojedynczej komórce. Każdy
fragment przekazuje dokładnie ten ładunek przez przeciętą ścianę Yee, który
znika z jednej pary wag CIC i pojawia się w drugiej. W rezultacie

\[
\frac{\rho^{n+1}-\rho^n}{\Delta t}+\nabla_h\cdot\mathbf J^{n+1/2}=0
\]

jest spełnione lokalnie, bez globalnego rozwiązania Poissona dla prądu.
Wieloosiowy test daje względną resztę ciągłości
\(6{,}25\cdot10^{-16}\). Po ośmiu krokach fali próżniowej względna reszta
\(\nabla_h\cdot\mathbf B\) wynosi \(2{,}41\cdot10^{-16}\), również bez
projekcji.

Relatywistyczny profil Gaussa jest na siatce Yee reprezentowany przez 125
znormalizowanych elementów materialnych. Ich przesunięcia między krokami
obejmują zmianę kontrakcji Lorentza, a każdy element używa tego samego lokalnego
operatora prądowego. Polaryzacja \(\mathbf P\) jest nanoszona na ściany
elektryczne, magnetyzacja \(\mathbf M\) na komplementarne krawędzie, po czym
źródła związane są tworzone jako
\(\rho_b=-\nabla_h\cdot\mathbf P\) oraz
\(\mathbf J_b=\Delta_t\mathbf P+\nabla_h\times\mathbf M\). Wspólna orientacja
operatorów incydencji zapewnia dokładne
\(\nabla_h\cdot(\nabla_h\times\mathbf M)=0\). Łączny test ruchomej chmury,
zmiennej kontrakcji oraz momentu magnetycznego daje względną resztę ciągłości
około \(1{,}8\cdot10^{-14}\) i zachowuje całkowity ładunek do dokładności
maszynowej.

Dwupoziomowa hierarchia Yee wykonuje dwa kroki poziomu dokładnego na jeden krok
poziomu grubego. Rejestr strumieni jest zamykany przez zachowawczą restrykcję:
gęstość ładunku jest średnią ośmiu objętości, a normalne składowe pól i prądu
są średnimi czterech współpokrywających się ścian. Jest to pełnopokryciowy test
referencyjny refluxingu, w którym każda ściana rodzica jest kontrolowaną granicą
AMR. Po sześciu krokach względna niezgodność restrykcji wynosi zero,
\(\nabla_h\cdot\mathbf B=0\), a ruchomy profil źródła zachowuje równanie
ciągłości z resztą około \(2{,}9\cdot10^{-14}\).

Zewnętrzny blok Yee ma także warstwę CPML. Każda z sześciu kierunkowych
pochodnych rotacji Faradaya i sześciu pochodnych rotacji Ampère’a ma osobną
zmienną pamięci splotowej. Profile \(\sigma\), \(\kappa\) i \(\alpha\) są
stopniowane sześciennie wzdłuż osi właściwej dla danej pochodnej; CPML nie jest
nakładane na wewnętrzne granice AMR. W teście impulsu po 180 krokach we wnętrzu
pozostaje \(0{,}0066\) energii początkowej, wobec \(0{,}4085\) w identycznym
bloku okresowym, a wewnętrzna reszta \(\nabla_h\cdot\mathbf B\) pozostaje
równa zeru.

Pętla cząstka–pole Yee zbiera każdą składową pola z jej własnych ścian lub
krawędzi metodą CIC. Relatywistyczny pusher Boris–Vay wyznacza nowe pędy i
położenia obu cząstek, po czym pełne trajektorie są przekazywane do operatora
Esirkepova. Dopiero po utworzeniu \(\rho^{n+1}\), \(\mathbf J^{n+1/2}\) oraz
prądu polaryzacji–magnetyzacji wykonywany jest krok Maxwella. Ta ścieżka nie
wywołuje projekcji Poissona: ograniczenie Gaussa jest propagowane przez
dyskretne równanie ciągłości. Test czterech pełnych kroków daje maksymalną
względną resztę ciągłości \(5{,}2\cdot10^{-13}\), ładunek pary zgodny z zerem
do dokładności maszynowej i maksymalną prędkość \(0{,}025c\). Osobny test mapy
pędu przód–tył daje zerową resztę w użytej precyzji.

W ten sposób wszystkie elementy nowego backendu — źródła, Yee, AMR, CPML i
pusher — mają wspólną ścieżkę referencyjną bez projekcji. Starszy solver
centrowany pozostaje dostępny dla dotychczasowych wykresów i testów
porównawczych; nie jest używany wewnątrz nowej pętli Yee.

Krótki tryb walidacyjny ma już dwukierunkowe sprzężenie cząstka–pole. Pole jest
interpolowane trójliniowo (CIC) do położeń środka profilu, a pęd aktualizuje
relatywistyczny, odwracalny w czasie pusher Boris–Vay. Nowe położenia wyznaczają
\(\rho^{n+1}\). Źródło prądowe jest składane jednorazowo z prądu konwekcyjnego,
pełnego prądu związanego oraz podłużnego dopełnienia Helmholtza wyznaczonego z
pary \(\rho^n,\rho^{n+1}\). Dopełnienie jest dyskretnym gradientem, więc nie
zmienia części poprzecznej ani rotacji naniesionego prądu, a z konstrukcji
spełnia na siatce

\[
\frac{\rho^{n+1}-\rho^n}{\Delta t}+\nabla_h\cdot\mathbf J^{n+1/2}=0.
\]

Nie jest to projekcja pola po kroku ani korekta bilansu po fakcie: rozwiązanie
zgodnego równania skalarnego jest częścią operatora nanoszenia źródła przed
aktualizacją Maxwella. Test sprzężony daje względną resztę ciągłości
\(2{,}04\cdot10^{-10}\), a rot dodanej części podłużnej względem skali
gradientu prądu \(8{,}24\cdot10^{-17}\). Osobna projekcja elektryczna nadal
kontroluje dryf dyskretnego prawa Gaussa samego solvera pola.

Pole siatkowe nie musi już startować od zera ani od chwilowego pola Coulomba.
Inicjalizator rozwiązuje dla każdej komórki równanie czasu retardowanego z tej
samej historii stanów co dynamika cząstek, a następnie sumuje pełne pola
Liénarda–Wiecherta elektronu i pozytonu: składnik prędkościowy \(1/R^2\) oraz
przyspieszeniowy \(1/R\). W jądrze profil punktowy jest mnożony przez dokładny
gaussowski udział ładunku zamkniętego

\[
F(R)=\operatorname{erf}\!\left(\frac{R}{\sqrt2r_c}\right)
-\sqrt{\frac{2}{\pi}}\frac{R}{r_c}
 \exp\!\left(-\frac{R^2}{2r_c^2}\right),
\]

więc pole pozostaje skończone i przechodzi w rozwiązanie punktowe daleko od
chmury. Elektryczna projekcja zachowuje poprzeczną część retardowaną, uzgadniając
część podłużną z dyskretnym rozkładem \(\rho\). Początkowe pole momentów
magnetycznych powstaje jako \(\mathbf B_\mu=\nabla_h\times\mathbf A_\mu\), co
gwarantuje brak monopoli magnetycznych konstrukcyjnie. Jest ono obecnie
regularnym polem quasi-statycznym; retardowane promieniowanie zmiennego momentu
wymaga historii \(\dot{\boldsymbol\mu}\) i \(\ddot{\boldsymbol\mu}\).

Test inicjalizacji przy niezerowych prędkościach i przyspieszeniach daje
względną resztę prawa Gaussa \(8{,}40\cdot10^{-6}\), względną resztę
\(\nabla_h\cdot\mathbf B\) równą \(3{,}08\cdot10^{-9}\) oraz niezerowy sygnał
składnika przyspieszeniowego \(1{,}46\cdot10^{-9}\) względem pola całkowitego.

Moment magnetyczny jest równocześnie obracany równaniem Thomasa–BMT w polu
siatkowym i zachowuje swoją normę. Ten backend pozostaje krótkim testem
walidacyjnym; standardowa wizualizacja nadal używa pól Liénarda–Wiecherta.
Na potrzeby backendu siatkowego dodano synchronizację poziomów i otwarte,
absorbujące warunki brzegowe opisane poniżej.

Poziomy są synchronizowane od najdrobniejszego do najgrubszego przez
objętościową restrykcję ośmiu komórek potomnych do jednej komórki nadrzędnej.
Podczas podcykli powłoka poziomu drobnego otrzymuje pole coarse→fine
interpolowane pomiędzy początkiem i końcem kroku grubego; po restrykcji powłoka
jest uzgadniana ponownie w tej samej chwili fizycznej. Na zewnętrznych komórkach
każdego bloku działa niesplittingowa warstwa convolutional PML (CPML). Dla
każdego kierunku i obu pól przechowywana jest osobna zmienna pamięci splotu.
Sześcienne profile \(\sigma\) i \(\kappa\), uzupełnione przesunięciem
częstotliwościowym \(\alpha\), modyfikują kierunkowe pochodne w rotacjach
Maxwella zamiast bezpośrednio tłumić amplitudy \(\mathbf E\) i \(\mathbf B\).
Pamięć CPML jest zerowana po przesunięciu łatki, aby nie przenosić historii
absorbera do innego obszaru fizycznego.

Dalekie pole i wypływ są definiowane wyłącznie na najgrubszym, wspólnym
poziomie, dzięki czemu fala nie jest księgowana osobno przez każdą łatkę.
Całki objętościowe mają postać kompozytową: komórki poziomu grubego przykryte
przez poziom drobny są pomijane. Na wewnętrznej powierzchni warstwy program
całkuje wektor Poyntinga oraz tensor
naprężeń Maxwella. Dostępne są więc szybkości wypływu energii, pędu i momentu
pędu, a energia faktycznie rozproszona w warstwie jest księgowana niezależnie.
Test czystego impulsu płaskiego \(E_y=cB_z\) osobno całkuje charakterystyki
biegnące w kierunkach \(+x\) i \(-x\). Dla ośmiokomórkowej CPML zmierzony
stosunek energii fali powracającej do wychodzącej wynosi
\(2{,}25\cdot10^{-32}\), a po teście w wewnętrznej części domeny pozostaje
\(3{,}08\cdot10^{-3}E_0\). Jest to test padania normalnego; kąty ukośne i
szerokie widmo wymagają osobnej serii walidacyjnej.

Bilans backendu siatkowego nie używa energii Darwina ani energii Schotta.
W każdej komórce całkowane są bezpośrednio

\[
U_{EM}=\int\left(\frac{\varepsilon_0E^2}{2}
+\frac{B^2}{2\mu_0}\right)dV,\quad
\mathbf P_{EM}=\varepsilon_0\int\mathbf E\times\mathbf B\,dV,
\]

\[
\mathbf J_{EM}=\varepsilon_0\int
\mathbf r\times(\mathbf E\times\mathbf B)\,dV.
\]

Do nich dodawane są relatywistyczna energia kinetyczna, pęd, orbitalny i
wewnętrzny moment pędu profili cząstek oraz skumulowane strumienie przez
powierzchnię absorbującą. Po aktywacji magnetyzacji test regresyjny wymaga
względnych reszt mniejszych niż \(10^{-3}\) dla energii, \(10^{-10}\) dla
pędu i \(10^{-3}\) dla momentu pędu. Obecne wartości wynoszą odpowiednio około
\(3{,}64\cdot10^{-4}\), \(2{,}16\cdot10^{-16}\) i
\(4{,}54\cdot10^{-4}\).

Dostępny jest operator nanoszący pełny, poruszający się dipol. Transformacja
\(M^{\mu\nu}\) daje laboratoryjną magnetyzację \(\mathbf M\) oraz indukowaną
polaryzację

\[
\mathbf P=\frac{\mathbf v\times\mathbf M}{c^2}.
\]

Na siatkę trafia więc pełny czteroprąd związany

\[
\rho_b=-\nabla_h\cdot\mathbf P,
\qquad
\mathbf J_b=\nabla_h\times\mathbf M+\partial_t\mathbf P.
\]

Pochodna polaryzacji jest liczona z sąsiednich stanów czasowych zgodnych z
predyktorem punktu środkowego. Niezależny test translacji dipola, wykonany
przed projekcją prądu swobodnego, daje względną resztę równania ciągłości
\(1{,}80\cdot10^{-16}\) oraz całkowity ładunek związany
\(|Q_b/e|=3{,}03\cdot10^{-21}\). Dla nieruchomego dipola pozostaje
\(\mathbf J_b=\nabla_h\times\mathbf M\), a test daje
\(|\nabla_h\cdot\mathbf J_b|\) na względnym poziomie około
\(6{,}75\cdot10^{-17}\). W sprzężonym teście prąd związany jest aktywny razem
z wzajemną siłą i momentem siły, liczonymi z tej samej gęstości Lorentza:

\[
\mathbf F_\mu=\int\mathbf J_M\times\mathbf B\,dV,
\qquad
\boldsymbol\tau_\mu=\int(\mathbf r-\mathbf r_c)\times
(\mathbf J_M\times\mathbf B)\,dV.
\]

Moment siły aktualizuje wewnętrzny moment pędu
\(\mathbf L_\mu=\boldsymbol\mu/\gamma_\mu\), który jest teraz częścią bilansu
\(\mathbf J\). Moc \(\int\mathbf J_b\cdot\mathbf E\,dV\), pomniejszona o
pracę translacyjną \(\mathbf F_\mu\cdot\mathbf v\), jest księgowana jako praca
więzów utrzymujących stałą normę klasycznego momentu. Jednostronny test bez
siły reakcji dawał około 20% reszty energii; po włączeniu pary wzajemnych
członów i zmniejszeniu kroku czasowego reszta wynosiła około 0,66%. Obecny
predyktor wyznacza stan cząstek oraz magnetyzację w połowie kroku. Pole
połówkowe służy do pełnego kroku pushera, \(\rho^{n+1}\) pochodzi z położenia
końcowego, a prąd konwekcyjny, \(\mathbf M\) i \(\mathbf P\) — ze stanu
środkowego. Po tej
zmianie reszta energii spadła do około 0,0364% bez korekty bilansu po kroku.

Test rozdzielczości zmniejsza względną resztę prawa Gaussa z około
\(2{,}32\cdot10^{-4}\) dla \(\Delta x=r_c\) do
\(2{,}31\cdot10^{-6}\) dla \(\Delta x=r_c/2\). Odbicie CPML nie jest już
szacowane z całej późnej energii: powyższy test rozdziela fale wychodzącą i
powracającą za pomocą zmiennych charakterystycznych na płaszczyźnie kontrolnej.

Solver Maxwell–Yee/AMR/CPML jest opcjonalnym backendem walidacyjnym i nie jest
kompilowany do zwykłego programu `positronium`. Osobny walidator buduje się i
uruchamia poleceniami:

```bash
make validation
./positronium_validation
```

Tryb raportuje ograniczenia Gaussa, ciągłości i \(\nabla\cdot\mathbf B\), pełne
bilanse cząstka–pole, test magnetyzacji, zbieżność przestrzenną, pochłanianie
granicy oraz koszt obliczeń. Dla bieżącej konfiguracji na testowej maszynie
otrzymano \(\Delta t_{CFL}=4{,}84\cdot10^{-22}\,\mathrm s\), około
\(2{,}07\cdot10^9\) kroków na pikosekundę i konserwatywnie około
\(8{,}56\cdot10^7\) sekund obliczeń na pikosekundę. Ten pomiar obejmuje także
koszt projekcji i testów pomocniczych, ale jednoznacznie wyklucza prostą zamianę
backendu wielopikosekundowej animacji.

Walidacja zbieżności ewoluuje tę samą periodyczną falę płaską do ustalonego
czasu na siatkach 12³ i 24³, zmniejszając równocześnie \(\Delta x\) i
\(\Delta t\) zgodnie z CFL. Względne błędy wynoszą odpowiednio
\(1{,}3437\cdot10^{-1}\) i \(3{,}3687\cdot10^{-2}\), co daje obserwowany rząd
zbieżności \(p=1{,}996\), zgodny z metodą punktu środkowego i centralnymi
pochodnymi drugiego rzędu.

Osobny test kowariancji porównuje pole Liénarda–Wiecherta ładunku poruszającego
się z prędkością \(0{,}35c\) z bezpośrednim boostem Lorentza spoczynkowego pola
Coulomba w tym samym zdarzeniu. Względna różnica pól wynosi
\(9{,}89\cdot10^{-16}\), a reszta niezmienników
\(\mathbf E\cdot\mathbf B\) i \(E^2-c^2B^2\) wynosi
\(1{,}83\cdot10^{-15}\). Test ten doprowadził również do zastąpienia
stałopunktowego wyznaczania czasu retardowanego metodą Newtona z pochodną
\(1-\mathbf n\cdot\boldsymbol\beta\).

Geometria produkcyjna zawiera teraz dwie niezależnie ruchome łatki o
\(\Delta x=r_c/2\), wycentrowane na elektronie i pozytonie, oraz jeden wspólny
blok dalekiego pola o \(\Delta x=4r_c\). Dla początkowej separacji \(a_0=100r_c\)
oba źródła znajdują się we własnych łatkach i zarazem we wspólnej domenie;
interpolator wybiera najdrobniejsze dostępne pole. Test przesunięcia jednej
cząstki potwierdza, że porusza się tylko właściwa łatka i pokrycie nie zostaje
utracone.

Rozgałęziona hierarchia ma wspólny zegar poziomu dalekiego pola. Jeden jego
krok odpowiada ośmiu podkrokom każdej łatki cząstkowej. W każdym podkroku
warunek brzegowy far→fine jest interpolowany w czasie, a po zakończeniu kroku
osiem komórek na oś jest uśrednianych objętościowo fine→far. Obie gałęzie są
restrykowane do tego samego bloku, po czym ich powłoki zostają ponownie
uzgodnione. Nakładanie łatek jest wykrywane i zatrzymuje krok, ponieważ w tej
sytuacji potrzebna jest jedna scalona łatka bliskiego pola, a nie podwójne
nanoszenie źródła.

Test wspólnej fali obejmuje 20 restrykowanych komórek poziomu dalekiego pola.
Po pełnym kroku rozgałęzionym względny defekt kompozytowej energii wynosi
\(6{,}06\cdot10^{-6}\), a pędu \(6{,}07\cdot10^{-6}\). Oba defekty są jawnie
zwracane przez operator sprzężenia, dzięki czemu nie znikają z diagnostyki
bilansu.

Podkroki dwóch gałęzi są wykonywane współbieżnie jako niezależne zadania CPU;
synchronizacja następuje dopiero przed wspólną restrykcją. Test uruchamia ten
sam krok również szeregowo i porównuje wszystkie komórki trzech bloków. Reszta
równoległość–szereg wynosi zero, a na dostępnych czterech logicznych CPU
zmierzono przyspieszenie około \(1{,}60\times\). Wynik czasowy jest informacją
diagnostyczną, nie sztywnym kryterium testu, ponieważ zależy od obciążenia
systemu.

Nie jest to jeszcze zgoda na domyślne przełączenie `Visual Simulation`.
Pozostaje produkcyjne przyspieszenie obliczeń. Najdrobniejszy CFL wymaga około
\(2{,}07\cdot10^9\) kroków na pikosekundę, więc wielopikosekundowa animacja bez
równoległości domenowej/GPU trwałaby lata. Do czasu usunięcia tych ograniczeń
`Visual Simulation` używa backendu Liénarda–Wiecherta, a tryb `maxwell` służy
do krótkich, samospójnych eksperymentów walidacyjnych.

Bieżąca maszyna nie udostępnia platformy obliczeniowej CUDA, HIP, SYCL ani
OpenCL (`clinfo --list` jest puste), mimo obecności urządzeń graficznych w PCI.
Z tego powodu kernel GPU i końcowe przełączenie nie są oznaczone jako
zweryfikowane. Aktywacja pełnego backendu w `Visual Simulation` wymaga najpierw
implementacji i regresji kerneli rotacji, projekcji oraz nanoszenia źródeł na
wspieranym urządzeniu; samo skompilowanie nieprzetestowanej gałęzi warunkowej
nie stanowiłoby walidacji fizycznej.

Pierwsza część tej optymalizacji jest już dostępna. Każdy poziom ma strefę
ochronną; po jej opuszczeniu środek łatki przesuwa się o całkowitą liczbę
komórek, pola w części wspólnej są przenoszone bez interpolacyjnego rozmycia,
a nowo odsłonięte komórki są zerowane. Poziomy poruszają się niezależnie, więc
mała zmiana położenia przesuwa tylko najdrobniejszą łatkę. Test zachowuje
\(0{,}99999248\) energii lokalnego pakietu przy regriddingu; strata pochodzi z
obciętego ogona Gaussa.

Ewolucja poziomów ma również rekurencyjne podcykle `2:1`: jeden krok poziomu
grubego odpowiada dwóm krokom następnego i czterem krokom poziomu
najdrobniejszego. Po każdym kroku grubym następuje restrykcja fine→coarse.
Po dodaniu czasowej interpolacji coarse→fine test wspólnej fali próżniowej daje
kompozytowy względny dryf energii \(2{,}87\cdot10^{-3}\), a maksymalna
znormalizowana różnica pola na uzgodnionej powłoce wynosi zero w precyzji
testu. Jawny refluxing strumieni ścianowych pozostaje niemożliwy w obecnym
układzie pól centrowanych w komórkach i będzie wymagał przejścia na
konserwatywny układ ścianowo-krawędziowy. Nowe dwie ruchome łatki pozostają
częściami jednego obiektu hierarchii i nie mogą zostać uruchomione jako dwie
niezależne domeny, które utraciłyby wzajemne pole i wspólny bilans.

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

### 2. Retardowane pola Liénarda–Wiecherta

Oddziaływanie orbitalne ładunek–ładunek jest liczone z pełnych pól
Liénarda–Wiecherta drugiej cząstki w czasie retardowanym. Na każdą cząstkę
działa siła Lorentza

\[
\mathbf F_i(t)=q_i\left[\mathbf E_j(t_{\rm ret})
+\mathbf v_i(t)\times\mathbf B_j(t_{\rm ret})\right],
\]

gdzie \(t_{\rm ret}=t-R(t_{\rm ret})/c\). Historia stanów jest przechowywana
w krótkim buforze. Położenie źródła jest interpolowane kubicznym wielomianem
Hermite’a wyznaczonym przez położenia i prędkości na końcach przedziału;
prędkość i przyspieszenie są analitycznymi pochodnymi tego samego wielomianu.
Zapewnia to ciągłość kinematyczną podczas iteracyjnego wyznaczania czasu
retardowanego. Przed początkiem symulacji przyjmowana jest gładka kontynuacja
stanu początkowego ze stałym początkowym przyspieszeniem. Pole zawiera zarówno
składnik prędkościowy \(1/R^2\), jak i przyspieszeniowy \(1/R\). W granicy
małych prędkości i zaniedbywalnego opóźnienia odzyskuje przyciąganie Coulomba.

### 3. Przybliżenie Darwina w diagnostyce

Lagrangian Darwina nie steruje już ruchem orbitalnym, ponieważ dublowałby
część magnetyczną retardowanych pól Liénarda–Wiecherta. Nadal służy jako
niskoprędkościowy punkt odniesienia w diagnostycznej energii i pędzie
kanonicznym. Oprócz członu Coulomba zawiera on

\[
L_D=\frac{kq_eq_p}{2c^2r}\left[
\mathbf v_e\cdot\mathbf v_p+
(\mathbf v_e\cdot\hat{\mathbf r})
(\mathbf v_p\cdot\hat{\mathbf r})\right].
\]

Jest to zachowawcze, chwilowe przybliżenie do rzędu \(v^2/c^2\). Nie stanowi
ścisłej energii całego układu retardowanego, ponieważ ta wymaga jawnego
uwzględnienia energii pola elektromagnetycznego w przestrzeni.

Odpowiadający człon energii, zapisany za pomocą prędkości, ma postać

\[
E_D=\frac{kq_eq_p}{2c^2r}\left[
\mathbf v_e\cdot\mathbf v_p+
(\mathbf v_e\cdot\hat{\mathbf r})
(\mathbf v_p\cdot\hat{\mathbf r})\right].
\]

Jest on uwzględniany wyłącznie w raportowanej energii diagnostycznej.

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
są losowane na początku, a następnie dynamicznie obracają się pod wpływem
lokalnego pola magnetycznego.

Pole momentu magnetycznego jest obecnie oceniane w czasie retardowanym. Z
potencjału

\[
\mathbf A=\frac{\mu_0}{4\pi}\left[
\frac{\boldsymbol\mu(t_r)\times\mathbf n}{r^2}
+\frac{\dot{\boldsymbol\mu}(t_r)\times\mathbf n}{cr}\right]
\]

wyznaczane są człony bliskiego pola \(1/r^3\), indukcyjny \(1/r^2\) i
promieniujący \(1/r\), a także indukowane pole elektryczne oraz ruchowy człon
\(-\mathbf v_\mu\times\mathbf B_\mu\). Pochodne momentu są obliczane z
retardowanej historii symetrycznie w jej wnętrzu i jednostronnie na końcach.
Pole lokalne jest sumą tego pola dipola oraz retardowanego pola
Liénarda–Wiecherta drugiego ładunku. Precesja jest opisana relatywistycznym
równaniem Thomasa–BMT. W czasie laboratoryjnym ma ono postać

Promieniowanie zmiennego momentu jest księgowane z mocą
\(P_\mu=\mu_0|\ddot{\boldsymbol\mu}|^2/(6\pi c^3)\). Odpowiadający mu
moment reakcyjny jest wyznaczany z \(\boldsymbol\mu\times
\dddot{\boldsymbol\mu}\), działa na orientację dipola, a energia jest
jednocześnie przenoszona z sektora wewnętrznego do energii promieniowania.
Bilans obejmuje także konwekcyjny pęd promieniowania dipolowego i wynoszony
moment pędu.

\[
\frac{d\boldsymbol\mu_i}{dt}
=\frac{q_i}{m_i}\boldsymbol\mu_i\times\mathbf B_{\mathrm{BMT},i},
\]

gdzie, dla \(a=(g-2)/2\),

\[
\mathbf B_{\mathrm{BMT}}=
\left(a+\frac1\gamma\right)\mathbf B
-a\frac{\gamma}{\gamma+1}
(\boldsymbol\beta\!\cdot\!\mathbf B)\boldsymbol\beta
-\left(a+\frac1{\gamma+1}\right)
\frac{\boldsymbol\beta\times\mathbf E}{c}.
\]

Model używa mierzonej wartości \(g\) danego gatunku, czyli
\(a=(g-2)/2=+0{,}00115965\) dla elektronu, a nie klasycznego \(g=1\)
(\(a=-1/2\)). Uwzględnia precesję Thomasa, wpływ pola elektrycznego i podłużnej
składowej pola magnetycznego. Znaki precesji elektronu i pozytonu są przeciwne.
Obrót jest wykonywany wzorem Rodriguesa w dwóch symetrycznych półkrokach
i zachowuje \(|\boldsymbol\mu_i|=(g/2)\mu_B\).

Niskoprędkościowym punktem odniesienia dla sprzężenia poruszającego się
ładunku z poruszającym się dipolem drugiej cząstki jest człon lagrangianu

\[
L_{q\mu}=q\,(\mathbf v_q-\mathbf v_\mu)\cdot\mathbf A_\mu,
\qquad
\mathbf A_\mu=\frac{\mu_0}{4\pi}\,
w(r)\frac{\boldsymbol\mu\times\mathbf r}{r^3}.
\]

W aktualnej ścieżce ruchu nie jest to jednak samodzielny, dokładny lagrangian
całego układu. Siła Lorentza korzysta z pól retardowanych, a siła nośnika
dipola z gradientu energii w jego chwilowym układzie spoczynkowym. Człon
Darwina występuje wyłącznie w diagnostyce energii i pędu Noethera. Z tego
powodu równość sił akcji i reakcji oraz energia kanoniczna tego przybliżonego
lagrangianu nie stanowią dowodu zachowania pełnego układu retardowanego.

Punktowe wzory magnetyczne są osobliwe dla \(r\to0\), dlatego pole orbitalne,
energia i siła dipolowa są płynnie tłumione funkcją

\[
w(r)=\frac{1}{1+(a/r)^6},\qquad
a^3=\frac{\mu_0}{4\pi}\frac{\mu_1\mu_2}{E/2}K_6.
\]

Przejście regulatora leży daleko poniżej skali atomowej. Jego wpływ na
najkrótszy, modelowy fragment trajektorii zależy od wybranej pary i orientacji
dipoli; nie jest ono utożsamiane z dodatkową granicą kończącą integrację.
Regularizacja jest nakładana na potencjał wektorowy, a wszystkie pola dipolowe
są wyznaczane z jednego związku

\[
\mathbf B_{\rm reg}=\boldsymbol\nabla\times\mathbf A_{\rm reg}
=w\mathbf B_0+\boldsymbol\nabla w\times\mathbf A_0.
\]

Energia dipol–dipol ma postać
\(U_{dd}^{\rm reg}=-\boldsymbol\mu_1\cdot\mathbf B_{2,\rm reg}\), a jej pełny
gradient daje siłę. Regularizacja jest założeniem numeryczno-modelowym, a nie
nowym prawem fizycznym.

### 5. Pełne pole dalekie ładunków i strumień Maxwella

Promieniowanie nie jest już redukowane do wzoru dipola elektrycznego. Dla obu
ładunków wyznaczane są pełne pola Liénarda–Wiecherta na sferycznym froncie
świetlnym: człon prędkościowy proporcjonalny do \(1/R^2\) oraz człon
przyspieszeniowy proporcjonalny do \(1/R\). Pola elektronu i pozytonu są
najpierw sumowane, dlatego interferencja wszystkich multipoli sektora
ładunkowego jest zachowana również wtedy, gdy rozmiar układu nie jest mały
względem długości fali.

Na sferze kontrolnej o promieniu \(10^6a_0\) wykonywana jest kwadratura po 50
równomiernie rozłożonych kierunkach. Z wektorów

\[
\mathbf S=\frac{1}{\mu_0}\mathbf E\times\mathbf B,
\qquad
T_{ij}=\varepsilon_0\left(E_iE_j-\tfrac12\delta_{ij}E^2\right)
+\frac{1}{\mu_0}\left(B_iB_j-\tfrac12\delta_{ij}B^2\right)
\]

program całkuje strumienie energii, pędu i momentu pędu:

\[
\dot E_{rad}=\oint\mathbf S\cdot\mathbf n\,dA,\qquad
\dot{\mathbf P}_{rad}=-\oint\mathbf T\cdot\mathbf n\,dA,\qquad
\dot{\mathbf J}_{rad}=-\oint\mathbf r\times
  (\mathbf T\cdot\mathbf n)\,dA.
\]

Są one całkowane w czasie metodą trapezów. Wspólny czas źródłowy oznacza, że
sfera jest frontem świetlnym kolejnych chwil emisji; dzięki temu nie jest
potrzebne przechowywanie całej historii aż do czasu przelotu przez dużą sferę.
„Pełne” dotyczy tu pola punktowych ładunków. Dynamiczne dipole magnetyczne
mają osobny, retardowany sektor indukcyjny i radiacyjny.

### 6. Indywidualna reakcja promieniowania Landaua–Lifshitza

Każdy ładunek otrzymuje osobną samosiłę w zredukowanym rzędzie. Zewnętrzna
czterosiła \(f_{\rm ext}^\mu\) zawiera wyłącznie pola drugiej cząstki i
oddziaływania dipolowe. Samosiła jest obliczana jako

\[
f_{\rm self}^\mu=\tau_i
\left(\delta^\mu_{\ \nu}-\frac{u_i^\mu u_{i\nu}}{c^2}\right)
\frac{d f_{\rm ext}^\nu}{d\tau},
\qquad
\tau_i=\frac{q_i^2}{6\pi\varepsilon_0m_ic^3}.
\]

Pochodna po czasie własnym jest wyznaczana centralnie wzdłuż przewidywanej
linii świata. Projekcja zapewnia \(u_\mu f_{\rm self}^\mu=0\), więc reakcja nie
zmienia masy spoczynkowej. Redukcja rzędu usuwa niezależne przyspieszenie z
równania i nie dopuszcza rozwiązań samoprzyspieszających typowych dla
niezredukowanego równania Lorentza–Abrahama–Diraca.

Poprzednie kolektywne pole proporcjonalne do \(\dddot{\mathbf p}\) zostało
usunięte. Wzajemne, interferencyjne oddziaływanie radiacyjne pozostaje w
retardowanych polach Liénarda–Wiecherta, natomiast samosiła zawiera wyłącznie
indywidualną reakcję danego ładunku.

### 7. Bilans pola związanego i kontrola reakcji radiacyjnej

Każdy zaakceptowany półkrok zapisuje osobno zmianę wielkości Noethera
cząstek, strumień przez daleką sferę oraz rezerwuar pola związanego i
interferencyjnego wewnątrz powierzchni kontrolnej:

\[
E_N+E_{rad}+E_{bound}=\mathrm{const},\quad
\mathbf P_N+\mathbf P_{rad}+\mathbf P_{bound}=\mathrm{const},\quad
\mathbf J_N+\mathbf J_{rad}+\mathbf J_{bound}=\mathrm{const}.
\]

Rezerwuar nie jest dodatkową siłą ani korektą trajektorii. Niezależnie
całkowana jest różnica między pracą indywidualnych sił Landaua–Lifshitza a
koherentnym strumieniem ładunkowym, więc zamknięcie bilansu nie ukrywa
niezgodności obu opisów.

Pomocniczo zachowana jest także energia bliskiego pola Schotta.

Do diagnostycznego bilansu energii dołączana jest niskoprędkościowa suma
indywidualnych energii Schotta

\[
E_S=-\sum_i\frac{q_i^2}{6\pi\varepsilon_0c^3}
\mathbf a_i\cdot\mathbf v_i,
\]

bez kolektywnych członów krzyżowych, które należą do wzajemnego pola
retardowanego. Reprezentuje ona odwracalną wymianę energii z bliskim polem
każdego ładunku. Raportowana wielkość kontrolna to

\[
E_{diag}=K_e+K_p+U_C+E_D+U_{dd}^{\rm reg}+E_{rad}+E_S.
\]

Wartość Schotta nie jest dodawana drugi raz do pełnego bilansu, ponieważ jest
częścią rezerwuaru pola związanego. Rezerwuar zawiera także zmianę wzajemnego
pola retardowanego, interferencję i residuum dyskretnego kroku.

### 8. Wielkości Noethera

Dla zachowawczej części przybliżonego lagrangianu program wyznacza pędy
kanoniczne. Obejmują one pęd relatywistyczny cząstki, pochodne członu Darwina
po prędkościach oraz przeciwne wkłady \(q\mathbf A_\mu\) sprzężenia
ładunek–dipol. Raportowane są

\[
\mathbf P_N=\mathbf p^{\rm can}_e+\mathbf p^{\rm can}_p,
\]

oraz

\[
\mathbf J_N=\mathbf r_e\times\mathbf p^{\rm can}_e
+\mathbf r_p\times\mathbf p^{\rm can}_p
+\frac{\boldsymbol\mu_e}{\gamma_e}
+\frac{\boldsymbol\mu_p}{\gamma_p}.
\]

Są to wielkości Noethera modelu przybliżonego, opisujące sektor cząstek i
efektywnego pola bliskiego. Diagnostyka porównuje teraz sumy
\(\mathbf P_N+\mathbf P_{rad}\) oraz
\(\mathbf J_N+\mathbf J_{rad}\), gdzie drugie składniki są skumulowanymi
strumieniami tensora Maxwella przez front kontrolny. Wartości te służą
diagnostyce oraz wykresom PDF; program nie zapisuje już archiwów ROOT.

### 9. Stochastyczne (skwantowane) promieniowanie dipolowe elektryczne

Sekcje 1–8 opisują **ciągłe** modele reakcji radiacyjnej: siła działa na
każdym kroku integratora, energia i pęd znikają płynnie. `--radiation-
reaction stochastic` (`ChargeRadiationReactionModel::stochasticElectricDipole`,
domyślny model produkcyjny) to inny wybór modelowania tego samego zjawiska:
promieniowanie jest emitowane w dyskretnych kwantach \(\hbar\omega\), zgodnie
z tym, jak rzeczywiście przenosi je pole — nie ciągły strumień, tylko
strumień fotonów. Continuous force jest tu wyłączona całkowicie
(`particleMultipoleRadiation` nie nakłada reakcji dla tego modelu); cała
fizyka poniżej żyje w sekularnym estymatorze `estimateCremCollapse`
(`modules/crem_collapse.hpp`), tym samym, który dla modeli ciągłych na
zmianę mierzy jeden obieg mechanicznie i analitycznie pomija do
\(200\,000\) kolejnych po zamkniętej obwiedni

\[
u(n)=u_0(1-Jx)^{-2/3},\qquad x=n/n_{skip}\in[0,1],\qquad J\le0{,}30,
\]

gdzie \(u=|E|\) (patrz "Związane pozytonium" wyżej dla tej samej techniki
elementów oskulacyjnych w modelach ciągłych). Poniższe podpunkty opisują,
co dokładnie ten model dokłada do tej wspólnej infrastruktury.

**A. Kiedy strzela foton: proces Poissona i luka ekspozycji.** Klasyczna moc
promieniowania \(P_{E1}(a,e)\) (ta sama `electricDipoleRadiatedPower`, na
której opiera się dipol koherentny) definiuje hazard \(\lambda=P/\hbar
\omega_{orb}\), gdzie \(\hbar\omega_{orb}=\hbar\cdot2\pi/T\) jest energią
"referencyjnego" fotonu przy bieżącej częstości orbitalnej — odpowiednik
korespondencyjny Bohra, nie założenie o konkretnym przejściu kwantowym.
Fotony są rozmieszczane techniką *thinning*: akumulowana jest całka
hazardu \(\int\lambda\,dt\), a zdarzenie następuje, gdy przekroczy ona
świeżo wylosowany próg \(\mathrm{Exp}(1)\) (metoda odwrotnej dystrybuanty).
Naiwna wersja — akumulator zerowany przy każdym pomiarze mechanicznym,
widzący tylko JEDEN zmierzony obieg na checkpoint — nie strzelała
**wcale**: hazard pojedynczego obiegu (\(\sim5{,}5\cdot10^{-5}\) przy
\(a=3{,}1\) pm) prawie nigdy nie przekracza progu, podczas gdy
\(200\,000\) pomijanych obiegów (których akumulator nigdy nie widział)
niosłoby \(\sim11\) fotonów oczekiwanego hazardu. Domknięte drugim,
niezależnym akumulatorem, całkującym hazard analitycznie po całym
pomijanym odcinku wzdłuż tej samej obwiedni \(u(n)\):

\[
\text{hazard}(J)=\frac{P_0T_0}{\hbar\omega_0}\,n_{skip}\,
\frac{3}{J}\Bigl(1-(1-J)^{1/3}\Bigr),
\]

zweryfikowane wobec siłowej kwadratury numerycznej (Python) do
\(10^{-12}\) względnie w całym zakresie \(J\), jaki kod kiedykolwiek
produkuje. Pozycja fotonu WEWNĄTRZ pomijanego odcinka (ułamek \(x\), a
stąd chwilowa energia orbity w chwili emisji) jest odzyskiwana odwróceniem
tej samej całki w zamkniętej formie.

**B. Energia fotonu.** \(E_\gamma=\hbar\omega_{orb}\cdot
(u(n_\gamma)/u_0)^{3/2}\) — referencyjny kwant przeskalowany do chwilowej
częstości orbitalnej w miejscu emisji. Odejmowana jest dokładnie z
budżetu orbitalnego i dokładnie tyle samo trafia do `radiatedEnergyTotal` —
zamknięty rachunek, sprawdzony na starcie tej pracy pod kątem podwójnego
liczenia: `radiatedEnergy`/`orbitalRadiatedEnergy` są już wypełniane
bezwarunkowo przez kwadraturę strumienia (`integrateElectrodynamicStep`),
bo para wciąż przyspiesza pod gołą siłą Coulomba nawet przy wyłączonej
reakcji — ręczny, dodatkowy zapis tej samej energii został znaleziony i
usunięty. Zmierzone bezpośrednio: pojedynczy foton potrafi nieść energię
porównywalną lub większą od całej bieżącej energii orbity (do
\(\approx18{,}5\times\) w zmierzonych przebiegach) — studnia kulombowska
jest bezdenna, więc to czyni orbitę głębiej związaną, nie łamie zachowania
energii, ale unieważnia każde przybliżenie zakładające mały skok (patrz
punkt E).

**C. Kierunek emisji.** Kąt \(\theta\) względem osi \(\hat{\mathbf L}\)
(orbitalny moment pędu, **nie** oś spinu kwantowego — to była wyjściowa
wątpliwość, która uruchomiła tę pracę) losowany jest z fizycznego wzorca
**wirującego** dipola elektrycznego,

\[
\frac{dP}{d\Omega}\propto1+\cos^2\theta,
\]

nie ze wzorca \(\sin^2\theta\) pojedynczego dipola liniowego (maksimum w
płaszczyźnie) — para krążąca po orbicie to dipol wirujący, nie oscylujący
liniowo. Odwrócenie dystrybuanty tego wzorca sprowadza się do sześcianu
zredukowanego \(\mu^3+3\mu+(4-8u)=0\) (\(\mu=\cos\theta\)), rozwiązywanego
w postaci zamkniętej wzorem Cardana; zweryfikowane numerycznie wobec
dystrybuanty do \(10^{-15}\) bezwzględnie. Azymut jest niezależnym losowaniem
jednostajnym; wspólna baza ortonormalna prostopadła do \(\hat{\mathbf L}\)
jest budowana raz i użyta zarówno dla kierunku fotonu, jak i (patrz E) dla
azymutu skrętności.

**D. Pęd liniowy i odrzut.** Foton niesie realny pęd \(\hbar\omega/c\).
Ponieważ warunki początkowe CREM są przygotowywane przy dokładnie zerowym
pędzie całkowitym (rozkład prędkości względnej wg stosunku mas), a każdy
model ciągły utrzymuje to zero z konstrukcji, zachowanie pędu sprowadza
się do jednorodnego kopnięcia całego układu: pełny, próbkowany w C
kierunek fotonu przesuwa trwałą prędkość środka masy
(`centreOfMassVelocity`), a wynikająca stąd zmiana energii kinetycznej CM
(\(v_{cm}\!\cdot\!p_\gamma+p_\gamma^2/(2M)\), policzona dokładnie, nie
odrzucona jako drugi rząd) jest doliczana do budżetu orbitalnego ponad
energię samego fotonu, żeby energia całkowita nadal spadała dokładnie o
\(E_\gamma\). Zmierzone: stosunek \(p_\gamma/p_{orbitalny}\) wychodzi
dokładnie jako (zredukowana długość fali Comptona pary)/\(a\), niezależnie
od prędkości — \(0{,}007\) na starcie orbity, ale już \(0{,}25\) przy
\(a=3{,}1\) pm i \(>1\) poniżej \(\approx0{,}77\) pm, więc realny, nie
pomijalny efekt w większości głębokości, jaką model osiąga.

**E. Moment pędu.** Dwie osobne wielkości, kierunek i wartość bezwzględna
\(|\mathbf L|\), są dziś aktualizowane JEDNYM mechanizmem, nie dwoma:

- *Odrzucony pierwszy pomysł (przechył orbitalny).* Wcześniejsza wersja
  przechylała \(\hat{\mathbf L}\) traktując odrzut jako kopnięcie
  \(r\times\Delta v\) ruchu względnego, niezależnie od odrzutu środka masy
  z punktu D. Sprawdzone numerycznie (trzy stosunki mas): jednorodne
  kopnięcie środka masy daje **dokładnie zerowy** moment siły
  (\(\sum_im_i\mathbf r_i=0\) względem środka masy, tożsamościowo) —
  jednorodne pchnięcie przez własny środek masy nie może układu skręcić.
  Przechył liczył więc ten sam pęd fotonu dwa razy; usunięty.
- *Odrzucone drugie podejście (współczynnik `k` z modeli ciągłych).*
  Klasyczny, orbitalnie uśredniony związek między tempem strat energii i
  momentu pędu dla reakcji dipolowej E1,
  \(k(e)=-(1-e^2)/(2+e^2)\) (metodą Petersa 1964, zaadaptowaną dla dipola
  EM), dawał \(L\mathrel{{*}{=}}(E_{po}/E_{przed})^{k}\) — przybliżenie
  pierwszego rzędu relacji różniczkowej \(d(\ln L)/d(\ln|E|)=k(e)\), z
  której pochodzi. Poprawione najpierw *w ramach tego samego podejścia*:
  rozdzielenie zmiennych w \(s=1-e^2\), \(x=\ln|E|\) daje całkę pierwszą w
  postaci zamkniętej, \((1-e^2)^3/(e^4|E|^3)=\mathrm{const}\), dokładnie
  zachowaną wzdłuż trajektorii — zweryfikowaną algebraicznie i niezależnie
  całkowaniem RK4 tego samego równania różniczkowego (zgodność
  \(10^{-12}\)). To ujawniło, że stara, zamrożona formuła dawała **ujemne**
  \(e^2\) po niemal każdym fotonie produkcyjnym (foton potrafi unieść do
  \(18{,}5\times\) bieżącej energii, gdzie przybliżenie małego skoku jest
  jakościowo złe), po cichu przycinane do zera przez istniejący strażnik.
  Mimo poprawienia, formuła `k` została ostatecznie odrzucona jako
  mechanizm STOSOWANY do stanu (zostaje tylko jako diagnostyka w logu
  `CREM_DEBUG`): moment siły reakcji scałkowany po orbicie i to, co unosi
  ciągłe pole, to ta sama wielkość z zachowania momentu pędu — łączenie jej
  z realnym momentem pędu skwantowanego fotonu (niżej) byłoby ponownie
  podwójnym liczeniem.
- *Zastosowany mechanizm: spin fotonu.* Foton to bezmasowy bozon spinu 1,
  więc niesie dokładnie \(\pm\hbar\) momentu pędu wzdłuż WŁASNEGO kierunku
  propagacji \(\hat{\mathbf n}\) (ten sam `photonDirection` z punktu C) —
  fakt uniwersalny, nie oszacowanie orbitalnie uśrednione. Skrętność
  \(h=\pm1\) losowana jest z rozkładu warunkowego (względem już
  wylosowanego \(\theta\)) dla standardowego przejścia dipolowego
  \(\Delta m=\pm1\):
  \[
  P(h{=}{+}1\,|\,\theta)=\frac{(1+\cos\theta)^2}{2(1+\cos^2\theta)},\qquad
  P(h{=}{-}1\,|\,\theta)=\frac{(1-\cos\theta)^2}{2(1+\cos^2\theta)},
  \]
  które sumują się dokładnie do wzorca \(1+\cos^2\theta\) z punktu C — nie
  nowe założenie, tylko jego rozdzielczość polaryzacyjna (sprawdzone: przy
  \(\theta=0\) i \(\theta=\pi\) skrętność w połączeniu z kierunkiem fotonu
  ZAWSZE daje wektor wzdłuż osi orbitalnej, zgodnie z \(\Delta m=1\); przy
  \(\theta=\pi/2\) to równy rozkład 50/50 w płaszczyźnie orbity).
  Zastosowane jako rzeczywiste kopnięcie wektorowe na pełnym (nie
  specyficznym) wektorze momentu pędu pary,
  \[
  \mathbf L_{para}\mathrel{-{=}}h\hbar\,\hat{\mathbf n},
  \]
  po czym wynik jest rozkładany z powrotem na wielkość
  \(|\mathbf L|\) i kierunek \(\hat{\mathbf L}\) — jeden mechanizm
  aktualizujący oba naraz, zamiast dwóch osobnych.

  Czego to NIE łapie: orbitalnego momentu pędu fotonu względem pary
  (\(\mathbf r\times\mathbf p_\gamma\), potrzebna nieznana anomalia
  prawdziwa — ta sama luka informacyjna reprezentacji samych elementów
  oskulacyjnych, która wykluczyła przechył). Sprawdzone wprost: wartość
  oczekiwana składowej osiowej \(\langle h\cos\theta\rangle\) po całym
  rozkładzie kątowym wynosi \(1/2\), nie \(1\) — sam spin odzyskuje więc
  średnio tylko połowę reguły wyboru \(\Delta m=1\); reszta jest
  niedostępna w tej architekturze.

  Zmierzone empirycznie: partia produkcyjna (o-Ps, \(80\) trajektorii, dwa
  ziarna) dała \(1\) awarię numeryczną wobec \(0\) awarii, ale \(3\) ucięć
  budżetem czasowym dla tej samej wielkości partii pod poprzednim kodem.
  Czas kolapsu ledwo drgnął (mediana orto \(130\)–\(155\) ps w kilku
  ziarnach, wobec \(147{,}8\)–\(151{,}6\) ps zmierzonych wcześniej z
  formułą `k`) — czas kolapsu wyznacza głównie całka energii/hazardu z
  punktów A–B, której ten mechanizm nie dotyka; mimośród (teraz naprawdę
  osiągający wartości jak \(e^2\approx0{,}9\)) wpływa na to, KTÓRY warunek
  wyjścia trajektoria trafi i jak szybko. Ta jedyna awaria doczekała się
  własnego, pełnego śledztwa — patrz punkt H.

**F. Co model świadomie zostawia otwarte.** (1) Orbitalny moment pędu
fotonu względem pary — wymaga anomalii prawdziwej, niedostępnej w
reprezentacji samych elementów oskulacyjnych; ilościowo to połowa reguły
\(\Delta m=1\) (punkt E). (2) Ścieżka mechaniczna
(`crem_trajectory.hpp`, używana rzadko produkcyjnie — głównie tryb
wizualny) ma własną, odrębną wersję hazardu i odrzutu fotonu, która NIE
otrzymała poprawek pędu liniowego ani spinu z tej sekcji: brakuje jej
pozycji/orientacji, jakiej te poprawki wymagają, i jest rzadko ćwiczona
produkcyjnie. (3) Zamrożone `k` pozostaje mechanizmem STOSOWANYM (nie
tylko diagnostyką) dla modeli ciągłych i dla gałęzi zbiorczej tego samego
estymatora (jej skok na checkpoint jest ograniczony do \(30\%\) energii,
gdzie błąd przybliżenia jest łagodny) — poprawka z punktu E dotyczy
wyłącznie pojedynczych zdarzeń fotonowych.

**G. Niezależna weryfikacja zastosowana na każdym etapie:** siłowa
kwadratura numeryczna w Pythonie dla całki hazardu (\(10^{-12}\)
względnie); porównanie z dystrybuantą dla odwrócenia Cardana kąta emisji
(\(10^{-15}\) bezwzględnie); całkowanie RK4 niezależne od wyprowadzenia
algebraicznego dla całki pierwszej \(k(e)\) (\(10^{-12}\)); sprawdzenie
trzech stosunków mas dla tożsamości \(\sum m_i\mathbf r_i=0\)
(\(10^{-17}\), szum numeryczny); odtwarzanie pojedynczej trajektorii z
partii przez jej własne ziarno (punkt H) do namierzenia rzadkich awarii
zamiast zgadywania ich przyczyny; partie produkcyjne po \(30\)–\(230\)
trajektorii mierzące rzeczywisty wskaźnik awarii, nie tylko argument
teoretyczny; `positronium_validation` 33/33 po każdej zmianie.

**H. Jedyna awaria numeryczna z partii testowej punktu E: znaleziona,
błędnie wyjaśniona za pierwszym razem, i naprawiona za drugim.** Historia
w trzech krokach, każdy zmierzony, żaden założony.

*Krok 1 — namierzenie konkretnej trajektorii.* Partia \(80\) trajektorii
(punkt E) dała \(1\) awarię, bez wskazania, która. `runCremCollapseExperiment`
wyprowadza ziarno trajektorii o indeksie \(i\) w partii jako
`splitMix64(masterSeed+i)` (`crem_collapse.hpp`), więc dowolne pojedyncze
zdarzenie z wielowątkowej partii — inaczej nie do odróżnienia w
przeplatanym logu wielu wątków naraz — daje się odtworzyć osobnym,
jednowątkowym przebiegiem `--runs 1 --seed (masterSeed+i)`. Przeszukanie
\(30\) wartości \(i\) dla partii ziarno \(99\) namierzyło winowajcę:
indeks \(8\), czyli ziarno \(107\).

*Krok 2 — pierwsze wyjaśnienie, zgadnięte, a nie zmierzone, i błędne.*
Pierwsza wersja tego dokumentu przypisała awarię temu, że "orbita, przez
pojedyncze duże kopnięcie spinowe, robi się na tyle mimośrodowa, że
opuszcza reżim związany" — czyli `elements.specificEnergy>=0`. Brzmiało
prawdopodobnie (kopnięcie spinowe jest w końcu zaburzeniem rzędu
jedności), ale nie zostało sprawdzone bezpośrednio na miejscu.

*Krok 3 — sprawdzone bezpośrednio, i inne.* Dodana diagnostyka
(`CREM_DEBUG`, dwa punkty `DIAG`, ciche bez tej zmiennej środowiskowej)
pokazała, że `elements.specificEnergy` w chwili awarii wynosiło
\(-1{,}1508\cdot10^{13}\) — głęboko ujemne, nie w pobliżu zera. Awarię
wywoływał inny, wcześniej istniejący strażnik:
`maxRelativeLossPerOrbit=0{,}5` w `crem_collapse.hpp`. Pojedynczy kop
spinowy zepchnął orbitę do \(e^2\approx0{,}945\); następny, pojedynczy
zmierzony mechanicznie obieg (ten sam pomiar, z którego każdy model
reakcji czerpie swoje tempo strat, oparty na prawdziwym strumieniu
retardowanego pola, niezależny od tego, który model reakcji jest
aktywny) wypromieniował w TYM JEDNYM obiegu \(56\%\) energii wiązania —
powyżej progu \(50\%\), przy którym ten sam plik już dokumentuje:
"sekularna inspirala nie może stracić dużej części własnej energii
wiązania w JEDNYM obiegu; gdyby mogła, uśrednianie po orbicie w ogóle by
nie obowiązywało".

*Krok 4 — pierwszy wniosek z kroku 3, też sprawdzony ponownie, na
wyraźne żądanie, i też okazał się za wcześnie wyciągnięty.* "Strażnik
zadziałał poprawnie, nie ma czego naprawiać" było naturalnym, ale
niesprawdzonym wnioskiem z samego faktu, że strażnik zareagował na
prawdziwe, duże, fizyczne zjawisko (a nie na jawny NaN). Sprawdzone,
co strażnik faktycznie CHRONI: dla `isStochastic` wielkość, którą testuje
(`deltaEnergyPerOrbit`), nie jest używana do NICZEGO poza samym
strażnikiem i osobno zabezpieczonym, nigdy niewracającym do stanu
diagnostycznym stosunkiem Larmora tuż niżej — do wyznaczenia
`orbitsToSkip` (czyli do ekstrapolacji na wiele pominiętych obiegów, to,
przed czym strażnik ma chronić) służy zamiast niej `expectedLossPerOrbit`,
analityczne tempo Larmora dla BIEŻĄCYCH elementów oskulacyjnych. Awaria,
przed którą strażnik chroni — skażony pomiar jednego obiegu
ekstrapolowany na wiele pominiętych — nie może więc tą ścieżką w ogóle
zajść dla tego modelu, niezależnie od wielkości zmierzonej straty. Sam
pomiar też nie jest podejrzany w sposób, w jaki był przypadek uzasadniający
strażnik pierwotnie (`--radiation-reaction coherent`, degenerujący się
stencil trzeciej pochodnej momentu dipolowego siły reakcji, dE/E o
\(22\) rzędy wielkości w jednym checkpoincie) — `stochasticElectricDipole`
nigdy tego stencila nie liczy (nie nakłada żadnej ciągłej siły reakcji,
patrz punkt A), a `orbitalRadiatedEnergy` pochodzi ze strukturalnie
niezwiązanej całki strumienia (`electromagneticFieldFluxRates` przez
`particleMultipoleRadiation`).

*Naprawa i weryfikacja.* Część progowa strażnika (nie część `isfinite`,
która nadal obowiązuje bezwarunkowo, dla każdego modelu, jako ochrona
przed prawdziwym uszkodzeniem numerycznym) jest teraz wyłączona dla
`isStochastic`. Sprawdzone, nie założone: dokładnie ta sama trajektoria
co w kroku 1 (ziarno \(107\)) teraz kończy się poprawnie
(\(172{,}78\) ps, zgodnie ze skalą już zanotowaną), bez żadnego
zadziałania strażnika. Powtórzona ta sama partia \(80\) trajektorii na
dwóch ziarnach, plus dwie kolejne partie (\(100\) i \(50\) trajektorii,
inne ziarna) dają \(0\) awarii numerycznych na \(230\) — wobec \(1/80\)
z nienaprawionym strażnikiem. Czas kolapsu ledwo drgnął (mediana orto
\(100\)–\(155\) ps w kilku ziarnach, wobec \(147{,}8\)–\(151{,}6\) ps
sprzed obu zmian z punktu E), bo wyznacza go całka energii/hazardu,
której żadna z tych zmian nie dotyka.

*Co z tego wynika ogólniej, poza tym jednym strażnikiem.* Dwa błędne
pierwsze wnioski z rzędu, każdy sprawdzony i poprawiony dopiero na
wyraźne żądanie ponownej weryfikacji, a nie z własnej inicjatywy — warte
zapisania wprost, nie tylko naprawienia po cichu.

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
  \(A=(\boldsymbol\mu_1\cdot\boldsymbol\mu_2)/(\mu_1\mu_2)\), czyli cosinusa
  kąta między momentami. Normalizacja przez iloczyn momentów obu ról, a nie
  przez \(\mu_B^2\), utrzymuje \(A\in[-1,1]\) także dla pary, w której momenty
  różnią się o rzędy wielkości.

Klasyfikacja przebiega następująco:

1. Ruch do wewnątrz i \(r_{min}<10^{-14}\,\mathrm m\) daje
   `Direct collision`.
2. W przeciwnym przypadku \(E_{rel}\ge 0\) daje `Scattering`.
3. Dla \(E_{rel}<0\), gdy \(A\ge0{,}5\), wynikiem jest
   `Para-positronium`.
4. Pozostały związany przypadek daje `Ortho-positronium`.

Próg orientacji dipoli jest fenomenologiczną regułą programu. Nie jest
klasycznym odpowiednikiem kwantowego singletu i trypletu.

`Cutoff time` w trybie wizualnym oznacza wyłącznie czas dojścia klasycznej
trajektorii do granicy rozdzielczości modelu. Dla zderzenia bezpośredniego
jest nią promień chmury ładunku \(0{,}01a_0\), a dla pozostałych przebiegów
\(10^{-14}\,\mathrm m\). Jeśli granica nie zostanie osiągnięta, program
wyświetla `not reached`. Wielkość ta nie jest
czasem życia ani czasem anihilacji i nie jest używana na wykresach zaniku
pozytonium.

## Tryb wizualny i statystyczny

Pierwsze pytanie programu wybiera jeden z dwóch trybów pracy:

1. `Visual simulation` wykonuje pojedynczy przebieg. Wariant `Line` pokazuje
   animację 3D z liniami torów, cząstkami i momentami magnetycznymi. Wariant
   `Dot` pokazuje wyłącznie próbki położenia: ciemnoniebieskie dla elektronu
   i pomarańczowe dla pozytonu.
2. `Statistical analysis` wykonuje domyślnie **1000 trajektorii w każdym
   eksperymencie statystycznym**.
   Eksperyment 5 (`Interactions`) wystrzeliwuje 100 par naprzeciw siebie
   z półgaussowskim parametrem zderzenia i **energią losowaną osobno dla każdej
   cząstki** z tego samego rozkładu \(N(\mu,\sigma)\), po czym
   klasyfikuje każdą trajektorię jako `Collision`, `Scattering`,
   `Para-Positronium` albo `Orto-Positronium` i podaje zestawienie liczności.
   Próg zderzenia to granica rozdzielczości modelu
   \(0{,}01\,a_0=529\,\mathrm{fm}\), a nie granica punktowa
   \(10\,\mathrm{fm}\): ta ostatnia leży poniżej klasycznej bariery
   dipol-dipol przy \(193\,\mathrm{fm}\), więc była nieosiągalna;
   szczegóły w `README.ap`, sekcja 5.4a.

   Niezależne losowanie obu energii oznacza, że **środek masy pary się
   porusza**. Dla układu dwóch ciał ze środkiem masy w spoczynku energie nie
   mogą być niezależne: zerowy pęd całkowity wymusza \(|p_1|=|p_2|\) i ustala
   drugą energię, gdy zna się pierwszą i obie masy. `--interaction-energy-ev`
   oznacza więc teraz energię kinetyczną **jednej cząstki w układzie
   laboratoryjnym**, a nie energię zderzenia. Energia zderzenia \(K_{CM}\)
   jest niezmiennikiem wyliczanym z obu pędów i raportowanym per zdarzenie.
   Po wyznaczeniu \(K_{CM}\) trajektoria jest transformowana i całkowana w
   układzie COM. Na sferze dopasowania wspólny pęd obu cząstek wynika z
   \(K_{CM}+k|q_1q_2|/r\), a pędy są równe i przeciwne. Dzięki temu geometria
   sfery, zegar orbit związanych i kryterium wychwytu należą do tej samej ramy;
   translacja środka masy nie może zmienić stanu związanego w niezwiązany.
   Parametr zderzenia (poprzeczny do osi wiązki) zachowuje wartość przy tym
   boostcie, a izotropowe kierunki dipoli są losowane w COM.

   Zależność od pary jest przy tym silna i warto ją znać przed doborem
   \(\mu\). Przy równych masach pędy prawie się znoszą, środek masy prawie nie
   rusza i \(\langle K_{CM}\rangle\) wychodzi około dwa razy większe od
   średniej na cząstkę — dla e⁺e⁻ przy \(N(0{,}6;0{,}4)\) eV jest to 1,27 eV.
   Przy p+e⁻ proton ma przy tej samej energii pęd 43 razy większy, więc niemal
   cała jego energia idzie w ruch środka masy i na zderzenie zostaje w praktyce
   sama energia elektronu: \(\langle K_{CM}\rangle = 0{,}66\) eV.
   Panele kinematyki fotonów są dokładnymi krzywymi referencyjnymi, a nie
   wynikiem Monte Carlo: generator anihilacji jest kwantową receptą niezależną
   od modelu klasycznego, a próbkowanie go odtwarzało jedynie rozkład, z
   którego losuje. Jego spójność jest sprawdzana w `positronium_validation`
   jako test `annihilation-generator`. Zestaw
   paneli zależy od wybranego eksperymentu; program nie wymusza tych samych
   czterech histogramów dla zjawisk o innej fizyce.

### Związane pozytonium

Dla przygotowanego p-Ps lub o-Ps eksperymentalną charakterystyką rozpadu jest
\(\Gamma=1/\tau\), a nie przekrój czynny. Obecny CREM nie zawiera kwantowego
operatora anihilacji. Statistical 1/2 oblicza zamiast tego operacyjny klasyczny
czas kolapsu, i to **mierzy** go mechanicznie, a nie ekstrapoluje wzorem.

Trajektoria nie jest orbitą zamkniętą, tylko inspiralą: okres kurczy się jak
\(T\propto a^{3/2}\), więc pojedyncza liczba go nie opisuje. Program raportuje
okres oskulacyjny

\[
T=2\pi\sqrt{\frac{\mu a^{3}}{k e^{2}}}
\]

na obu końcach przebiegu oraz liczbę obiegów między nimi. Bezpośrednie
całkowanie wszystkich \(\sim10^{5}\) obiegów byłoby zbyt kosztowne, więc jeden
reprezentatywny obieg jest rozwiązywany pełnym silnikiem CREM, jego **zmierzona**
strata energii i momentu pędu daje bieżące tempo dysypacji, a ograniczona liczba
kolejnych obiegów jest pomijana analitycznie przy zamrożonym tempie, zanim
nastąpi ponowny pomiar. To standardowa technika elementów oskulacyjnych. Bieg
kończy się, gdy **perycentrum** osiągnie \(0{,}1\,a_0\); pozostały odcinek do
\(0{,}01\,a_0\) jest świadomie obcinany, bo przy \(t\sim a^3\) odpowiada za
0,1% czasu kolapsu.

Każde zdarzenie ma budżet zegarowy (`--crem-wallclock-budget-s`). Trajektoria,
która go wyczerpie, jest **cenzurowana prawostronnie**, a nie odrzucana: wiadomo
o niej, że czas kolapsu przekracza osiągnięty czas symulowany. Krzywa przeżycia
jest dlatego estymatorem iloczynowym Kaplana-Meiera z błędami Greenwooda,
liczonym ze wszystkich użytecznych trajektorii. Gdy cenzurowana jest cała
próba, estymator pozostaje skończony: \(S(t)=1\) do największego czasu cenzury,
RMST jest równe temu horyzontowi, a mediana pozostaje nieosiągnięta. Nie jest to
ekstrapolacja czasu życia poza dane.

To rozróżnienie jest istotne ilościowo. Cenzura **nie** jest niezależna od
mierzonej wielkości: czas kolapsu skaluje się jak \(a^3\), a liczba obiegów do
scałkowania jak \(a^{3/2}\), więc budżet zatrzymuje przede wszystkim orbity
najszersze, czyli najwolniej zapadające się. Uśrednianie samych przebiegów
ukończonych zaniża wynik — w teście kontrolnym na tym samym ziarnie dawało
1,72 ps wobec 3,60 ps przy pełnej kompletacji. Ponieważ cenzura jest
informatywna, sam Kaplan-Meier też nie jest nieobciążony: po jej rozpoczęciu
estymator nie widzi kolejnych kolapsów, więc \(S(t)\) jest podtrzymywane i RMST
zawyża. Oba obciążenia mają **przeciwne znaki**, więc program podaje obie
liczby wraz z frakcją ukończenia i poniżej 90% ostrzega, że wartość jest tylko
zakresem wiarygodności. Przy 100% kompletacji RMST i zwykła średnia są sobie
równe. Wnioski ilościowe wymagają budżetu dającego kompletację bliską 100%.

Po poprawieniu skali orbity na \(a_{Ps}\) oraz dołożeniu brakującego członu
wzajemnej reakcji promieniowania zmierzony czas kolapsu wynosi **33,9 ps** przy
stanie związanym o energii \(-6{,}3\,\mathrm{eV}\) i \(L=1{,}03\,\hbar\).
Dla porównania: przy poprzedniej, zbyt ciasnej orbicie było 3,6 ps, a po samej
poprawce skali — 76 ps, przy czym ta ostatnia wartość zawierała jeszcze
dwukrotne niedopromieniowanie.

Ponownie zmierzone na produkcyjnym przebiegu N=1000, ziarno 42, e⁺e⁻, bez pola
zewnętrznego (18 minut na 4 wątkach), już po przełączeniu celu ekstrapolacji z
`collisionBoundaryRadius` na barierę Comptona (`e0e5b4d`): **1000/1000
trajektorii dotarło do granicy, zero cenzury**, mediana Kaplana-Meiera
\(30{,}76\) ps, średnia (= RMST przy pełnej kompletacji) \(33{,}88\pm0{,}44\)
ps, rozrzut sigma/średnia \(0{,}410\) — nierozróżnialne od wartości sprzed tej
zmiany co do trzeciej cyfry znaczącej, dokładnie zgodnie z argumentem
"obcięcie ostatniego odcinka to ~0,1% czasu kolapsu" z sekcji o podłodze CREM.
Ten sam przebieg dla orto-pozytonium (eksperyment 2, też 18 minut, też
1000/1000) daje **statystycznie tę samą liczbę**, \(33{,}88\) ps — bezpośrednie
potwierdzenie zastrzeżenia, że oba kanały różnią się tu jedynie wymuszonym
wyrównaniem dipoli, którego sprzężenie (\(\approx10^{-5}\) potencjału
Coulomba) jest za słabe, by rozróżnić klasyczną inspiralę.

Ten sam bieg daje trzy dalsze porównania z zamkniętą elektrodynamiką i
pomiarem, każde na całej próbie N=1000: stosunek zmierzonego czasu kolapsu do
zamkniętego wzoru klasycznej inspirali \(\langle t_{\rm CREM}/t_{\rm
classical}\rangle=1{,}018\pm0{,}014\) (`1_1_3`), stosunek zmierzonej mocy
promieniowania do wzoru Larmora dla koherentnego dipola elektrycznego
\(\langle P_{\rm CREM}/P_{\rm Larmor}\rangle=0{,}9927\pm0{,}0116\) (`1_1_4`)
— oba w granicach 2% jedynki, zero parametrów swobodnych — oraz klasyczne
sprzężenie dipol-dipol przygotowanej pary \(\langle|U_{dd}|/h\rangle=7{,}875\)
GHz, pokrywające \(3{,}87\%\) mierzonego rozszczepu nadsubtelnego o-Ps/p-Ps
\(203{,}3941\) GHz (`1_2_2`) — reszta to anihilacja wirtualna i człon
kontaktowy Fermiego, których model klasyczny nie niesie.

Ponownie przeliczone po serii ośmiu poprawek warstwy sekularnej CREM
(`ecf7380`…`0ebc7d2`, N=1000, e⁺e⁻, bez pola zewnętrznego, bez ustalonego
ziarna, 16,4 i 16,9 minuty na 4 wątkach dla p-Ps/o-Ps): **1000/1000
trajektorii dotarło do granicy w obu kanałach, zero cenzury i zero
`NumericalFailure`/`not decaying`** — pierwszy raz od uruchomienia tego
przebiegu bez ani jednej awarii numerycznej na całej próbie w obu kanałach
naraz. Mediana Kaplana-Meiera p-Ps \(31{,}89\) ps, średnia \(34{,}40\pm0{,}44\)
ps (rozrzut sigma/średnia \(0{,}407\)); o-Ps mediana \(31{,}90\) ps, średnia
\(34{,}42\pm0{,}44\) ps (rozrzut \(0{,}408\)) — statystycznie ta sama liczba
w obu kanałach, jak poprzednio, i w granicach szumu innego ziarna
nierozróżnialne od przebiegu sprzed tej serii poprawek. Trzy dalsze
porównania: \(\langle t_{\rm CREM}/t_{\rm classical}\rangle=1{,}018\pm0{,}014\)
w obu kanałach (`1_1_3`), \(\langle P_{\rm CREM}/P_{\rm
Larmor}\rangle=0{,}9987\pm0{,}012\) (p-Ps) i \(0{,}999\pm0{,}014\) (o-Ps)
(`1_1_4`), sprzężenie dipol-dipol \(7{,}976\) GHz / \(3{,}92\%\) (p-Ps) i
\(7{,}655\) GHz / \(3{,}76\%\) (o-Ps) rozszczepu nadsubtelnego (`1_2_2`) —
wszystko zgodne z poprzednim przebiegiem w granicach 1-2 punktów procentowych,
dokładnie tego rzędu co różnica losowego ziarna. Poprawki tej serii dotyczyły
wyłącznie tego, jak głęboko sekularna księgowość CREM potrafi bezpiecznie
zajść pod barierą Comptona (fundament, warstwa sekularna, sygnał energii),
nie samej fizyki mierzonej w normalnym zasięgu produkcyjnym — zero awarii na
1000/1000 w obu kanałach jest bezpośrednim, produkcyjnym potwierdzeniem tej
poprawy odporności, nie tylko wynikiem diagnostyki opisanej wyżej.

Czwarta poprawka tej serii (`f477866`…`ff48f08`, opisana szczegółowo w
"Podłoga zasięgu CREM: bariera Comptona zamiast granicy zderzenia" wyżej)
usunęła kink w pierwszej pochodnej siły na granicy zacisku (zmiękczenie
Plummera) i zastąpiła stały próg zatrzymania `finalApproachMultiple` warunkiem
podwójnym: periapsis \(\le\) bariera Comptona LUB okres/czas-przelotu-światła
\(\le150\), zależnym od momentu pędu orbity, nie od stałego promienia.
Ponownie zmierzone na produkcyjnym przebiegu N=1000 dla obu kanałów: p-Ps
1000/1000 czystych (zero błędów, zero cenzury), mediana Kaplana-Meiera
\(30{,}76\) ps, średnia \(33{,}62\pm0{,}44\) ps; o-Ps 1000/1000 czystych,
mediana \(31{,}73\) ps, średnia \(34{,}12\pm0{,}44\) ps — statystycznie
nierozróżnialne od siebie i od przebiegu sprzed tej poprawki, dokładnie jak
przewiduje argument "ucięty odcinek to ~0,1% czasu kolapsu": zmiana progu
zatrzymania nie mogła zmienić mierzonej fizyki, tylko jej numeryczną
niezawodność.

Zmierzony czas życia p-Ps to 125,1 ps, więc pozostała rozbieżność to czynnik
\(\approx3{,}7\). Świadomie nie próbujemy jej domknąć: silnik odtwarza teraz
zamknięty wzór klasycznej inspirali z dokładnością \(1{,}8\%\)
(`N_1_3_collapse_time_vs_theory.pdf`), więc dalsze zbliżanie się do 125 ps
wymagałoby zepsucia elektrodynamiki, a nie jej poprawiania. Klasyczna inspirala
i anihilacja \(2\gamma\) to różne procesy. Czas życia o-Ps pozostaje poza
zasięgiem: 142 ns odpowiadałby orbicie \(a=33a_0\), czyli \(n\approx4\), a nie
stanowi podstawowemu.

Wynik nie jest kwantowym czasem anihilacji. Zewnętrzne pomiary nie są wejściem
CREM; tworzą wyłącznie krzywe porównawcze

\[
P_{\rm exp}(t)=\tau_{\rm exp}^{-1}e^{-t/\tau_{\rm exp}},
\]

z parametrami pochodzącymi z cytowanych pomiarów. Osobny generator produktów
zachowuje energię i pęd w układzie spoczynkowym pozytonium.

Dla obu kanałów wyświetlany jest ten sam zestaw czterech rozkładów i dwóch
diagnostyk:

- krzywa przeżycia kolapsu CREM (Kaplan-Meier) i krzywa eksperymentalna;
- widmo czasu anihilacji z mierzonej stałej rozpadu, jako odniesienie skali;
- zmierzony czas kolapsu wobec zamkniętego wzoru klasycznej inspirali;
- stosunek zmierzonej mocy promieniowania do wzoru Larmora dla dipola
  koherentnego;
- rozkład mocy promieniowania uśrednionej po trajektorii;
- klasyczne sprzężenie dipol-dipol wobec mierzonego rozszczepu nadsubtelnego.

Dwa środkowe panele i ostatni są testami: pierwsze dwa wobec zamkniętej
elektrodynamiki klasycznej, ostatni wobec pomiaru. Odniesienia nie zawierają
żadnego składnika CREM, więc odchyłka jest wypowiedzią o silniku, a nie o
pozytonium.

Są to wielkości mierzalne, ale na poziomie „truth”: idealna próżnia, stan w
spoczynku, brak rozdzielczości i akceptancji detektora, zderzeń z materiałem,
pick-off, quenchingu i poprawek radiacyjnych.

### Niezwiązana para e⁺e⁻

Statystyka stanów niezwiązanych jest osobnym eksperymentem wiązkowym. Dla
zadanej energii kinetycznej \(K_{CM}\) program losuje parametr zderzenia równomiernie po
powierzchni,

\[
b=b_{max}\sqrt{U},\qquad A=\pi b_{max}^2,
\]

i całkuje trajektorię od sfery dopasowania do ponownego wyjścia z obszaru
oddziaływania. Każda z \(N\) prób ma wagę \(A/N\).

Wspólny moduł kinematyki wyznacza dokładny relatywistyczny pęd COM dla obu
mas. Na sferze dopasowania cząstki dostają równe i przeciwne pędy, a nie
masowo dzieloną prędkość względną. Dzięki temu zadane \(K_{CM}\), całkowity
pęd równy zeru i wynik są niezmiennicze przy zamianie kolejności ról, również
dla pary proton–elektron.

Wyświetlane są:

- różniczkowy przekrój elastyczny \(d\sigma/d\Omega\) z błędami dwumianowymi;
- skumulowany przekrój akceptancji \(\sigma(\theta\ge\theta_{min})\);
- karta ustawień, przekroju fiducjalnego i liczby nieudanych trajektorii.

Ten sam schemat obsługuje dwa eksperymenty: 4 to rozpraszanie szerokie, 3 to
kanał krótkiego zasięgu. Dla eksperymentu 3 obowiązuje ostrzeżenie o zakresie
stosowalności podane przy opisie jego plików wyjściowych — z tego kanału nie
wolno czytać wielkości energetycznych, bo bilans energii modelu się w nim
rozjeżdża.

Dla kanału krótkiego zasięgu dochodzi modelowe widmo strat energii
\(d\sigma/d\Delta E\), gdzie \(\Delta E=K_{CM}-E_{out}\). Zakres osi obejmuje
energię wiązki i jest w razie potrzeby poszerzany do pełnej próbki. Dla
dalekiego rozpraszania elastycznego program
celowo nie rysuje wąskiego histogramu \(E_{out}\): bez modelu odpowiedzi
detektora byłaby to delta przy \(K_{CM}\), a jej mikroskopijna szerokość w
obliczeniu pochodzi przede wszystkim z dryfu numerycznego.

Koszyki kąta są dobrane tak, aby w granicy Rutherforda miały porównywalną
oczekiwaną statystykę. Na wykresach znajduje się krzywa Rutherforda jako
kontrola granicy małych prędkości. Całkowity przekrój Coulombowski jest
rozbieżny w przód, dlatego program zawsze raportuje przekrój z jawnym progiem
kątowym. Zdarzenie dochodzące do granicy Collision (`BeamConfiguration::
collisionRadius`: bariera Comptona \(193{,}3\) fm dla e⁺e⁻, `nuclearCutoff`
\(=10^{-14}\,\mathrm m\) dla każdej innej pary — ta sama warunkowa reguła co w
CREM i w eksperymencie 5, ujednolicona w `340b67b`) trafia do osobnego kanału
`reach cutoff`; jego przekrój nie jest przekrojem anihilacji QED.
Pełna relatywistyczna interpretacja danych e⁺e⁻ wymaga amplitudy rozpraszania
Bhabhy, której klasyczny integrator jeszcze nie zawiera.

Trajektorie wiązkowe są liczone równolegle, natomiast ROOT tworzy wykresy po
ich zakończeniu w głównym wątku. Sto zdarzeń daje szybki podgląd; do gładkiego
\(d\sigma/d\Omega\) należy zwiększyć `--runs` co najmniej do około 1000.

Klasyczne tryby Visual i Statistical używają jednego obiektu
`ClassicalTrajectoryEngine`. Wspólne są równania ruchu, historia retardowana,
reakcja promieniowania Landaua–Lifshitza, precesja momentów oraz księgowanie
strumieni pola. Statistical zmienia wyłącznie sposób próbkowania warunków
początkowych, kryterium zakończenia eksperymentu wiązkowego i agregację
wyników; nie ma własnego integratora fizycznego. Test regresyjny uruchamia dwa
odbiorniki od identycznego stanu i wymaga bitowo identycznego wyniku.

Walidacja integratora zawiera również niezależny przebieg po boostcie
Lorentza (0{,}35c), prostopadłym do początkowej orbity. Zdarzenia na liniach
świata są porównywane po transformacji czasoprzestrzennej, a czterosiła po
transformacji czterowektorowej. Skończona całka promieniowania na sferze ma
osobne residuum: sfera przy stałym czasie nie transformuje się w sferę przy
stałym czasie drugiego układu, dlatego nie jest ono ścisłym testem
czterowektora. Pełny test promieniowania wymagałby wspólnej kowariantnej
hiperpowierzchni świata.

Stan cząstki ma również pełną dekompozycję tensora dipolowego na elektryczny
moment \(\mathbf p\) i magnetyczny moment \(\boldsymbol\mu\). Transformacja
Lorentza miesza obie części zgodnie z

\[
\mathbf p'=\gamma\left(\mathbf p+\frac{\mathbf V\times\boldsymbol\mu}{c^2}\right)
-\frac{\gamma^2}{\gamma+1}\frac{\mathbf V(\mathbf V\cdot\mathbf p)}{c^2},
\]

oraz analogicznym wzorem dla \(\boldsymbol\mu'\). Test kontroluje transformację
tam i z powrotem oraz dwa niezmienniki tensora. Produkcyjne warunki początkowe
nie mają trwałego elektrycznego dipola w układzie własnym; dotychczasowy
ruchowy człon pola pozostaje jedynym sprzężeniem, aby nie liczyć go podwójnie.

Elektryczna składowa tensora ma własne opóźnione pole punktowego dipola:
statyczny człon \(1/r^3\), indukcyjny \(1/r^2\) i radiacyjny \(1/r\).
Pochodne \(\dot{\mathbf p}\) i \(\ddot{\mathbf p}\) są liczone z tej samej
historii co moment magnetyczny. Pole jest dołączane do lokalnego pola precesji
oraz siły Lorentza na drugą cząstkę. Test boostu trajektorii zawiera niezerowy
moment magnetyczny i indukowany moment elektryczny.

Niezależnym stopniem swobody jest obecnie moment magnetyczny w chwilowym
układzie spoczynkowym cząstki. Jest on zamieniany na czterowektor dipolowy
ortogonalny do czteroprędkości i ewoluowany względem czasu własnego
czterowymiarowym równaniem Thomasa–BMT z mierzonym \(g\) danego gatunku,
osobno dla każdej z dwóch ról. Krok RK4
kończy projekcja \(S\cdot u=0\) oraz odtworzenie stałej normy własnej. Pełny
tensor jest następnie boostowany do układu laboratoryjnego i wyznacza jednocześnie
\(\mathbf p\) oraz \(\boldsymbol\mu\). Poprzedni ręczny człon
\(-\mathbf v_\mu\times\mathbf B_\mu\) został usunięty: ruchowe pole elektryczne
wynika teraz wyłącznie z elektrycznej składowej tego samego tensora.

Siły translacyjne nośników dipoli są wyznaczane z jednego kowariantnego
sprzężenia. Pole drugiej cząstki jest transformowane do chwilowego układu
spoczynkowego nośnika, gdzie energia ma postać
\(U_0=-\boldsymbol\mu_0\cdot\mathbf B_0\). Jej symetryczny gradient
przestrzenny daje siłę, przeliczaną następnie na laboratoryjną siłę
trójwymiarową. Zastępuje to dawne, oddzielne wzory siły dipol–dipol oraz siły
reakcji na nośnik w sprzężeniu ładunek–dipol. Siła Lorentza na ładunek nadal
wynika bezpośrednio z pełnego pola drugiego źródła.

### Dopasowania i sekcja `Experimental`

Adnotacje na wykresach rozdzielają trzy różne rzeczy: obliczenie CREM, opisowy
fit jego wyników oraz niezależny wynik eksperymentalny. Generator kinematyki
fotonów nadal korzysta z Monte Carlo, ale nie losuje czasu życia.

- czasy kolapsu CREM otrzymują opisowe dopasowanie
  \(N(t)=N_0\exp(-t/\tau_{\rm CREM})\). Nie jest ono prawem anihilacji i nie
  opisuje kształtu próbki: wszystkie trajektorie startują z \(r=a_0\), losowana
  jest wyłącznie prędkość początkowa, więc rozkład czasów kolapsu jest wąski i
  jednomodalny, a nie wykładniczy. Raportowana niepewność jest błędem
  standardowym średniej \(\sigma/\sqrt{N}\) z rzeczywistej próbki; panel podaje
  obok iloraz \(\sigma/\langle t\rangle\), żeby widać było, jak dalece kształt
  odbiega od wykładniczego. Krzywa wykładnicza jest rysowana wyłącznie
  poglądowo. Oś czasu jest skalowana do danych, a nie do \(\tau_{\rm exp}\):
  dla o-Ps obie skale dzieli około osiem rzędów wielkości i wymuszanie
  \(\tau_{\rm exp}\) w zakresie osi sprowadzało całą próbkę CREM do jednego
  piksela. Zewnętrzne krzywe porównawcze wykorzystują
  \(\tau_{p\text{-}Ps}=125.142(27)\,\mathrm{ps}\) oraz
  \(\tau_{o\text{-}Ps}=142.037(26)\,\mathrm{ns}\) z precyzyjnych
  pomiarów szybkości zaniku;
- idealna linia p-Ps jest wzorcem
  \(\delta(E-\mu)\), bez sztucznie fitowanej szerokości detektora. Benchmark
  energii to \(510.9955485(15)\,\mathrm{keV}\), obliczony z masy pozytonium
  zestawionej przez PDG;
- izotropię p-Ps opisuje
  \(C[1+a_2P_2(\cos\theta)]\), z idealnym benchmarkiem \(a_2=0\). Wpisy obu
  fotonów są skorelowane, dlatego estymator \(a_2\) korzysta z jednego
  oznaczonego fotonu na zdarzenie;
- widmo, Dalitz i kąt o-Ps używają ustalonego wzorca Ore’a–Powella bez
  swobodnych parametrów kształtu. Nie są
  zastępowane arbitralnym wielomianem ani Gaussem. Pomiar kontinuum 3γ jest
  zgodny z QED, lecz bez tabel wydajności i odpowiedzi konkretnego detektora
  nie ma uczciwego bezpośredniego overlay;
- dla wiązki fitowana jest jedna bezwymiarowa normalizacja \(C_R\) kształtu
  Rutherforda. Skumulowany wykres pokazuje projekcję tego samego fitu, a nie
  drugie dopasowanie silnie skorelowanych progów. Widmo \(\Delta E\) kanału
  cutoff nie dostaje Gaussa ani Landaua, ponieważ żadna z tych rodzin nie
  opisuje zaimplementowanego procesu;
- rozkłady reszt numerycznych mają jedynie opisowy fit Gaussa w wyświetlanej
  zmiennej. Ich sekcja `Experimental` ma wartość `n/a`, ponieważ są kontrolą
  obliczeń, a nie obserwablą detektora.

Dla idealizowanej swobodnej pary e⁺e⁻ przy aktualnej energii i cięciach nie
jest wbudowany zgodny zbiór danych eksperymentalnych. Panel mówi wtedy wprost
`no matching dataset loaded`. Krzywa Rutherforda pozostaje referencją
analityczną, nie danymi. Porównanie z eksperymentem relatywistycznym wymagałoby
amplitudy Bhabhy oraz identycznego układu odniesienia, akceptancji i poprawek
radiacyjnych.

### Druga strona: weryfikacja numeryczna

Każde okno statystyczne ma pasek z przyciskami `NEXT` i `STOP`. `NEXT`
przenosi na drugą stronę i zmienia nazwę na `BACK`; `BACK` przywraca rozkłady,
a `STOP` zamyka okno. Obie strony powstają z tych samych zdarzeń. Przełączanie
nie uruchamia ponownie obliczeń i nie zmienia równań ruchu.

Dla generatora 2γ/3γ druga strona pokazuje reszty w jednostkach precyzji
maszynowej \(\epsilon\):

\[
\frac{\sum E_\gamma-E_{Ps}}{E_{Ps}\epsilon},\qquad
\frac{c|\sum\mathbf p_\gamma|}{E_{Ps}\epsilon},\qquad
\frac{E_\gamma^2-c^2p_\gamma^2}{E_\gamma^2\epsilon},
\]

oraz resztę masy niezmienniczej odtworzonego czteropędu rodzica. Są to testy
domknięcia energii, pędu i warunku bezmasowego fotonów. Całkowity moment pędu
rozpadu nie jest testowany, ponieważ zdarzenie nie przechowuje helicytetów,
polaryzacji ani początkowej macierzy gęstości spinu; orbitalne
\(\sum\mathbf r\times\mathbf p\) ze wspólnego wierzchołka byłoby tautologią.

Dla wiązki porównywane są wyłącznie trajektorie `escaped`, których stan
początkowy i interpolowany stan końcowy leżą na tej samej sferze
\(R_{match}\). Zdarzenia cutoff i time-gate nie są mieszane z tą próbą.
Wyświetlane są

\[
\delta_E=\frac{\Delta(E_N+E_{rad}+E_{bound})}{E_{scale}},\qquad
\delta_P=\frac{|\Delta(\mathbf P_N+\mathbf P_{rad}+\mathbf P_{bound})|}{P_{scale}},
\qquad
\delta_J=\frac{|\Delta(\mathbf J_N+\mathbf J_{rad}+\mathbf J_{bound})|}{J_{scale}}.
\]

\(E_N\), \(\mathbf P_N\) i \(\mathbf J_N\) są wielkościami Noethera
przybliżonego modelu. Wielkości `bound` księgują objętościową część pola
bliskiego i interferencyjnego. Residuum LL–strumień jest raportowane osobno.

### Automatyczny zapis PDF

Każdy rzeczywisty wykres (bez paneli zawierających wyłącznie tekst) jest
automatycznie zapisywany do katalogu `distributions`. Nazwa ma kolejno numer
eksperymentu, numer ekranu, numer padu i opis, rozdzielone znakami `_`, na
przykład:

```text
distributions/1_1_1_crem_collapse_time.pdf
distributions/2_1_4_three_photon_dalitz.pdf
distributions/3_1_3_energy_loss_cross_section.pdf
distributions/4_1_1_differential_cross_section.pdf
```

Diagnostyki drugiej strony również są zapisywane jako osobne pliki. Każdy PDF
zawiera dokładnie jeden pad i jeden odpowiadający mu rozkład; plik nazwany
`diagnostic_energy_balance` nie zawiera pozostałych diagnostyk. Zbiorcze
PDF całych stron Statistical nie są tworzone. Kolejne uruchomienie atomowo
nadpisuje pliki o tych samych nazwach; program najpierw sprawdza poprawność
nowego PDF, aby nie utracić poprzedniego obrazu przy błędzie renderowania.
Eksport działa także z `--no-gui`.

#### Katalog generowanych plików

Poniżej wszystkie pliki, jakie program może zapisać. Zestaw zależy od
wybranego eksperymentu; dwa wpisy są warunkowe i zaznaczono to osobno.

**Statistical 1 i 2 — para- i orto-pozytonium** (po 6 plików; `N` = 1 lub 2)

Oba eksperymenty mają teraz identyczny zestaw paneli: różnią się wyłącznie
wymuszonym wyrównaniem dipoli przy losowaniu i jednostką czasu (ps / ns).

Panele fotonowe (energia 2γ, kąt biegunowy, widmo Ore'a-Powella, wykres
Dalitza, kąt między wiodącymi fotonami) zostały **usunięte**. Były dokładnymi
krzywymi kwantowymi i nie zawierały żadnego wyniku CREM, więc nie dawały się
z niczym porównać. Zwolnione pady zajmują porównania modelu z zamkniętymi
wzorami elektrodynamiki i z pomiarem. Samospójność generatora anihilacji jest
nadal sprawdzana w `positronium_validation` jako test `annihilation-generator`.

| Plik | Zawartość |
| --- | --- |
| `N_1_1_crem_collapse_time.pdf` | Krzywa przeżycia kolapsu CREM, estymator Kaplana-Meiera. Schodki ze znacznikami cenzury i słupkami Greenwooda; krzywa `exp(-t/τ_exp)` z pomiaru jako odniesienie skali. Oś czasu przełącza się na logarytmiczną dopiero powyżej ~1,3 dekady rozpiętości. |
| `N_1_2_collapse_time_distribution.pdf` | Histogram **zmierzonej** próby czasów kolapsu CREM (ten sam surowy `decayTimes`, co panel 1, tu bez cenzury Kaplana-Meiera) na osi liniowej czasu i logarytmicznej liczby zdarzeń — konwencja realnych widm spektroskopii czasu życia pozytonu. Pionowa kropkowana linia = τ_exp, narysowana tylko gdy mieści się w zasięgu osi (dla o-Ps zwykle nie: podpis wtedy mówi "off-scale"); to WYŁĄCZNIE odniesienie skali, nie test predykcji. |
| `N_1_3_collapse_time_vs_theory.pdf` | **Porównanie z teorią.** Zmierzony czas kolapsu wobec zamkniętego wzoru klasycznej inspirali `da/dt = −C/a²`, `C = 8ke⁴/(6πε₀c³m²)`, uśrednionego czynnikiem **dipolowym** `(1+e²/2)/(1−e²)^{5/2}` przy własnych `a` i `e` każdej trajektorii. Zero parametrów swobodnych, żaden składnik CREM nie wchodzi do odniesienia. Przerywana przekątna = zgodność dokładna. |
| `N_1_4_radiated_power_vs_larmor.pdf` | **Porównanie z teorią.** Stosunek zmierzonej mocy dysypacji orbitalnej do larmorowskiej mocy koherentnego dipola elektrycznego dla tej samej orbity oskulacyjnej. Linia ciągła przy 1 = dipol koherentny, kropkowana przy 0,5 = dwa ładunki promieniujące niezależnie. |
| `N_2_1_diagnostic_calibration_power.pdf` | Histogram mocy promieniowania uśrednionej po trajektorii, tylko dla przebiegów zakończonych na granicy. |
| `N_2_2_dipole_coupling_vs_hyperfine.pdf` | **Porównanie z pomiarem.** Rozkład klasycznej energii oddziaływania dipol-dipol przygotowanej pary, wyrażonej jako częstość, zestawiony z mierzonym rozszczepem nadsubtelnym o-Ps/p-Ps 203,3941 GHz. Panel podaje, jaki procent rozszczepu pokrywa człon klasyczny; reszta to anihilacja wirtualna i człon kontaktowy Fermiego, których model klasyczny nie zawiera. |

**Statistical 3 i 4 — wiązka e⁺e⁻** (`N` = 3 lub 4)

Eksperyment 3 to kanał krótkiego zasięgu (`shortRangeFocus`), 4 — rozpraszanie
elastyczne. Poza doborem akceptancji kątowej i `b_max` różnią się jednym
plikiem: panel strat energii powstaje **wyłącznie w eksperymencie 3**, więc 3
daje 6 plików, a 4 daje 5.

> **Ostrzeżenie o zakresie stosowalności — eksperyment 3.** Audyt kompletności
> fizycznej (sekcja *Wynik audytu kompletności fizycznej*) pokazuje, że kanał
> krótkiego zasięgu wypada poza zakres stosowalności modelu wedle jego własnej
> diagnostyki, i to we wszystkich trzech niezależnych miarach naraz: rezerwuar
> pola związanego sięga \(1{,}4\)–\(9{,}4\) energii samej orbity, człon Schotta
> jest \(3{,}4\)–\(3{,}7\) raza większy od energii wypromieniowanej, a praca
> reakcji Landaua–Lifshitza rozmija się ze strumieniem dalekiego pola
> o czynnik \(6{,}2\)–\(6{,}4\). Nie jest to efekt małego mianownika:
> bezwzględnie przy \(E_{rel}=-7{,}8\) eV wychodzi \(E_{bound}=+11{,}0\) eV
> i \(E_{Schott}=-1{,}5\) eV wobec \(0{,}45\) eV wypromieniowanych.
>
> Dla porównania kanał związany (1, 2) i szerokie rozpraszanie (4) są w tych
> samych miarach czyste, na poziomie \(10^{-4}\) energii orbity: eksperyment 4
> na produkcyjnym przebiegu N=1000 (ziarno 42, 3,5 minuty) daje 1000/1000
> trajektorii bez awarii i dopasowanie \(C_R=0{,}966\pm0{,}024\) (95% Wilson
> \([0{,}918;\,1{,}012]\)) — zgodne z czystym Coulombem w granicach błędu.
>
> **Ten sam produkcyjny przebieg N=1000 dla eksperymentu 3 (ziarno 42, 6,6
> minuty) najpierw potwierdził ostrzeżenie wprost, dosłownie**: 8 z 1000
> trajektorii stało się numerycznie nieskończone, program sam odrzucił cały
> raport przekroju czynnego jako nieważny i nie zapisał żadnego wykresu —
> nie tylko ostrzeżenie w tekście, ale twardy brak wyniku. 40,7% zdarzeń
> trafiło do kanału `captured` zamiast `escaped`, co samo w sobie jest
> sygnaturą reżimu, w którym te ostrzeżenia obowiązują.
>
> **Zbadane i naprawione.** Zdiagnozowane bezpośrednio (te same 8 trajektorii,
> instrumentacja na każdym miejscu awarii): integrator wyczerpuje głębokość
> podziału kroku (`maximumDepth`, 8 przy pierwszej próbie, 12 przy retry) w
> trakcie zbliżania, **zanim** istniejący próg `nuclearCutoff` (10 fm,
> analogiczny do granicy `Collision` w eksperymencie 5) w ogóle dostanie
> szansę zadziałać — wszystkie 8 awarii nastąpiło przy promieniu 14–27 fm,
> ciągle *powyżej* 10 fm, ale głębiej niż bariera Comptona (193,3 fm), gdzie
> CREM już wcześniej wykazał kruchość. Retry z ciaśniejszą tolerancją odzyskał
> tylko 1 z 9 — to fizyka, nie tolerancja.
>
> Naprawa: w miejscu awarii, **zamiast poddać się**, program liczy analitycznie
> peryapsis Keplera z zachowanej energii i momentu pędu ostatniego
> rozstrzygniętego stanu (ten sam wzór, którego `crem_collapse.hpp` już używa
> do decyzji sekularnych z samego \((E,L)\)) — jeśli wynik jest już
> \(\le\)`nuclearCutoff`, trajektoria i tak była nieuchronnie skazana na wynik
> `Collision`, więc zgubienie ostatnich kilku fm rozdzielczości nie powinno
> unieważniać całego przekroju czynnego z jej powodu. Zastosowane tylko, gdy
> trajektoria wciąż leci do środka (`!passedClosestApproach`) — po minięciu
> peryapsis `minimumSeparation` jest zmierzonym faktem, nie wnioskiem.
>
> **Zmierzone po naprawie, ten sam bieg**: `failed: 8 -> 0`,
> `collision: 0 -> 8` — te same osiem zdarzeń, poprawnie sklasyfikowane
> zamiast odrzucone. Raport przekroju czynnego jest teraz ważny:
> \(\sigma(\text{reach cutoff})=26061\pm9177\) barn wobec czysto
> kulombowskiego odniesienia \(22622\) barn — zgodne w granicach błędu.
> Eksperyment 4 (gdzie żadna trajektoria nigdy nie schodzi tak głęboko) daje
> po tej zmianie bit-w-bit ten sam wynik co przedtem: \(C_R=0{,}965625\pm
> 0{,}0240074\), zero awarii. Ostrzeżenie o rezerwuarze pola związanego
> (poniżej) i tak nadal się uruchamia na tym samym przebiegu — poprawka
> dotyczy klasyfikacji awarii integratora, nie zmienia oceny bilansu
> energetycznego kanału.
>
> Wartości \(1{,}4\)–\(9{,}4\) pochodzą z pojedynczej trajektorii
> diagnostycznej. **Zespół wiązkowy wypada łagodniej**, bo dominują w nim
> dalekie przeloty: mediana \(|E_{bound}|/|E_{rel}|\) wynosi \(0{,}11\),
> a 95. percentyl \(0{,}35\), wobec \(7{,}1\cdot10^{-8}\) i
> \(1{,}8\cdot10^{-7}\) dla eksperymentu 4 — nadal sześć rzędów różnicy.
> Program **sam to teraz wykrywa i ostrzega w czasie działania**, progiem na
> 95. percentylu, a nie zakazem wpisanym na numer eksperymentu: ostrzeżenie
> zadziała dla dowolnej pary i dowolnych ustawień, które wejdą w ten sam reżim,
> i milczy, gdy nie wejdą.
>
> W praktyce znaczy to, że **z eksperymentu 3 nie wolno czytać wielkości
> energetycznych**. Dotyczy to wprost panelu strat energii poniżej, a także
> interpretacji panelu bilansu energii, którego niezależne residuum
> LL-vs-strumień jest dla tego kanału rzędu jedności, a nie promila. Rozkłady
> kątowe \(d\sigma/d\Omega\) i \(\sigma(\theta\ge\theta_{min})\) opierają się
> na geometrii rozproszenia, nie na bilansie energii, więc to ostrzeżenie ich
> bezpośrednio nie unieważnia — ale i one pochodzą z trajektorii liczonych
> w tym samym reżimie.

| Plik | Zawartość |
| --- | --- |
| `N_1_1_differential_cross_section.pdf` | Różniczkowy przekrój `dσ/dΩ` z błędami dwumianowymi, na tle Rutherforda i dopasowania `C_R × Rutherford`. |
| `N_1_2_cumulative_cross_section.pdf` | Skumulowany przekrój `σ(θ ≥ θ_min)`. Bez niezależnego dopasowania — progi są skorelowane, więc rzutowane jest `C_R` z panelu różniczkowego. |
| `3_1_3_energy_loss_cross_section.pdf` | **Tylko eksperyment 3.** Różniczkowy przekrój po `ΔE = K_CM − E_out` w oknie fiducjalnym. **Nie czytać jako pomiaru energetycznego** — patrz ostrzeżenie o zakresie stosowalności powyżej: dla tego kanału bilans energii modelu się rozjeżdża. |
| `N_2_1_diagnostic_energy_balance.pdf` | Residuum **tożsamości** bilansu energii. To nie jest test zachowania: `E_bound` jest zdefiniowane jako reszta domykająca sumę. Panel podaje osobno niezależne residuum fizyczne LL-vs-strumień — i to właśnie ono jest dla eksperymentu 3 rzędu jedności zamiast promila, patrz ostrzeżenie powyżej. |
| `N_2_2_diagnostic_momentum_balance.pdf` | To samo dla pędu, w skali `log₁₀`. |
| `N_2_3_diagnostic_angular_momentum_balance.pdf` | To samo dla momentu pędu. |

**Statistical 5 — oddziaływania** (7 plików)

| Plik | Zawartość |
| --- | --- |
| `5_1_1_outcome_summary.pdf` | Klasyfikacja zakończeń: zderzenie, rozproszenie, para-Ps, orto-Ps, nierozstrzygnięte, awaria numeryczna. |
| `5_1_2_collision_energy.pdf` | Wyliczona niezmiennicza energia w układzie środka masy, z zaznaczeniem podpróbki związanej. |
| `5_1_3_impact_parameter.pdf` | Losowany parametr zderzenia (rozkład półnormalny), z podpróbką związaną. |
| `5_1_4_dipole_alignment.pdf` | Wyrównanie dipoli w stanach związanych, z progiem para/orto przy `+0,5`. Stosunek 1:3 wynika tu z **geometrii** progu na izotropowej sferze, a nie ze statystyki spinowej — zbieżność liczbowa jest przypadkowa. |
| `5_1_5_collapse_time_distribution_para.pdf` | Histogram zmierzonych czasów kolapsu CREM schwytanych par sklasyfikowanych jako para-Ps (zwykle nieliczne: tylko te, których okno obserwacji pokryło pełny obieg Keplera). Oś prowadzona danymi, τ_exp jako pionowa linia tylko gdy się mieści w zasięgu. |
| `5_1_6_collapse_time_distribution_ortho.pdf` | To samo dla klasy orto. |
| `5_2_1_diagnostic_summary.pdf` | Panel tekstowy z medianami residuów bilansu i zastrzeżeniem o definicji `E_bound`. |

**Tryb wizualny** (1 plik na eksperyment, `N` = 1…4)

| Plik | Zawartość |
| --- | --- |
| `N_1_1_visual_simulation.pdf` | Zrzut całego ekranu animacji: scena, tabela i przyciski. |

Po zakończeniu animacji trybu wizualnego cały ekran (scena, tabela i przyciski)
jest zapisywany jako `N_1_1_visual_simulation.pdf`, gdzie `N` jest numerem
eksperymentu. Kliknięcie `EXIT` zapisuje w tym samym pliku aktualny stan kanwy
przed zakończeniem programu, również gdy animacja nie dobiegła jeszcze końca.

### Katalog źródeł

Program nie tworzy archiwów `.root`. Wyniki przeznaczone do oglądania są
zapisywane jako PDF, a centralny katalog referencji zasilający adnotacje jest
po każdym przebiegu statystycznym zapisywany do
`ScientificalReferences.txt`. Plik tekstowy zawiera wartości pierwotne,
niepewności statystyczne/systematyczne, wartości pochodne, wzory, założenia,
cytowania, DOI i adresy źródeł. Kolejne uruchomienie nadpisuje go atomowo.

## Struktura źródeł

Nagłówki w `modules/` są **modułami tekstowymi**, nie osobnymi jednostkami
kompilacji: każdy jest włączany dokładnie raz, wewnątrz anonimowej przestrzeni
nazw `positronium.cpp`, i korzysta z tego, co zdefiniowano przed nim. Dzięki
temu cały program pozostaje jedną jednostką kompilacji, a podział służy
czytelności i rozdzieleniu odpowiedzialności.

| Moduł | Zawartość |
| --- | --- |
| `vector3.hpp`, `state.hpp`, `dipole_tensor.hpp` | typy podstawowe |
| `physical_constants.hpp` | stałe fizyczne wraz z wyprowadzeniami |
| `particle_species.hpp` | tablica gatunków (masa, ładunek, `g`), para integrowana przez przebieg i skale, które z niej wynikają: masa zredukowana, promień Bohra pary, energia wiązania, promień regularyzacji dipola, granica zderzenia oraz wyszukiwanie gatunku dla `--pair` |
| `two_body_kinematics.hpp` | stabilna relatywistyczna kinematyka dwuciałowa: czteropędy, energia niezmiennicza, boost Lorentza i stan wejściowy na sferze dopasowania |
| `electrodynamics.hpp` | **prawa sił**: pola opóźnione, Darwin, sprzężenie dipolowe, reakcja promieniowania, strumień dalekiego pola |
| `crem_engine.hpp` | **numeryka**: rekonstrukcja historii przyczynowej i adaptacyjny integrator (sonda błędu, podział kroku, retencja historii) |
| `crem_trajectory.hpp` | warstwa trajektorii: `Frame`, wspólna pętla całkowania, sampler warunków początkowych i klasyfikacja zjawiska |
| `crem_collapse.hpp` | estymator kolapsu: całkowanie sekularne z uśrednianiem po orbitach oraz zamknięte odniesienia elektrodynamiczne |
| `statistics_archive.hpp` | katalog wartości zmierzonych i teoretycznych |
| `root_export.hpp` | atomowy zapis PDF |
| `maxwell_validation*.hpp` | zestaw testów budowany do `positronium_validation` |

Trzy moduły `crem_*` nie zawierają **żadnego kodu ROOT** — cała prezentacja
pozostaje w `positronium.cpp`. Bieżąca walidacja ma 33 nazwane bramki; trzy
pierwsze sprawdzają kinematykę dwuciałową, w tym zamianę ról nierównych mas,
boost i boost odwrotny oraz warunek \(|v|<c\). Osobna bramka kontrolera
adaptacyjnego wymusza odrzucenie kroku, który na `maximumDepth` nadal
przekracza tolerancję, oraz sprawdza atomowy rollback stanu i historii.

## Metoda numeryczna

Równania są całkowane jawnym schematem kick–drift–kick. Siły
Coulomba, dipolowe, lagrangianowe ładunek–dipol i Darwina wykonują dwa półkroki pędu,
a precesja dipoli jest dzielona na dwa półkroki. Pomiędzy
półkrokami wykonywany jest krok położenia. Krok czasu jest adaptacyjny. Tryb
Visual proponuje krok nie większy niż \(5\times10^{-18}\,\mathrm s\) i 1/128
chwilowego okresu Coulombowskiego, a następnie dzieli go według lokalnej
estymaty błędu. Animacja przechowuje około 240 klatek i kończy się wcześniej
po przekroczeniu granicy modelu. Próby wiązkowe nie przechowują klatek. Statystyka zaniku stanów
związanych korzysta z niezależnego generatora zdarzeń i nie całkuje
klasycznej orbity aż do cutoffu.
Przekroczenie progu zderzenia jest lokalizowane przez interpolację wewnątrz
ostatniego kroku. Klatka znajdująca się już poniżej granicy ważności modelu nie
jest przekazywana do animacji ani diagnostyki.

Funkcja wykonująca jeden pełny krok nazywa się
`integrateElectrodynamicStep`: nazwa podkreśla, że mutuje stan i wykonuje
integrację, precesję dipoli, reakcję promieniowania oraz akumulację energii
promieniowania. Samo obliczanie zachowawczych oddziaływań pozostaje oddzielone
w `allExternalForces`.

Ponieważ siły Darwina i ładunek–dipol zależą od prędkości, schemat wymaga
osobnego badania zbieżności i nie jest ogólnie integratorem symplektycznym dla
całego lagrangianu. Obie cząstki są rysowane jako kulki tej samej wielkości.

## Wielkości prezentowane na ekranie

Widok ROOT pokazuje trajektorie i wektory momentów magnetycznych. Strzałki
`Spin ↑/↓` w tabeli są jedynie skrótem graficznym: wskazują znak składowej
\(z\) klasycznego momentu magnetycznego, a nie pomiar kwantowego spinu.

Nad tabelą widoczne są: wynik klasyfikacji, \(E_{rel}\), \(L_{orb}\),
\(r_{min}\) i diagnostyczny `Cutoff time`. Tabela bieżąca zawiera czas, odległość, energię
mechaniczną oraz skumulowaną energię promieniowania. Wykres położenia używa
jednostek \(a_0\), odległość jest podawana w pm, czas w ps, a energie w eV.

## Uruchomienie i odtwarzalność

Wymagany jest CERN ROOT z programem `root-config` dostępnym w `PATH`.

```bash
make build
./positronium
```

Samo `make` kompiluje program i od razu go uruchamia.

Pierwsze menu wybiera sposób uruchomienia:

```text
1 -> Visual simulation
2 -> Statistical analysis
```

W trybie wizualnym wybierany jest najpierw sposób prezentacji, a następnie
zjawisko:

```text
1 -> Line
2 -> Dot
```

```text
1 -> Para-positronium
2 -> Ortho-positronium
3 -> Direct collision
4 -> Scattering
```

W trybie statystycznym te same numery wybierają odpowiednio kolaps CREM dla
konfiguracji p-Ps i o-Ps wraz z kinematyką fotonów, wiązkowy kanał krótkiego
zasięgu i elastyczne rozpraszanie
wiązki e⁺e⁻. Interaktywne menu przy każdej pozycji podaje też zmierzony,
nie zgadnięty, orientacyjny czas przy N=1000: eksperymenty 1/2 (pełne
mechaniczne całkowanie każdej trajektorii do granicy) ~18 minut każdy,
eksperyment 3 ~7 minut, 4 ~3,5 minuty, 5 ~5 minut — na 4 wątkach, dla e⁺e⁻ bez
pola zewnętrznego; skaluje się w przybliżeniu liniowo z `--runs`. Tryb, wybór
eksperymentu i ziarno można podać bez interakcji:

```bash
./positronium --mode visual --visual-style line --phenomenon 2 --seed 42
./positronium --mode visual --visual-style dot --phenomenon 2 --seed 42
./positronium --mode statistical --phenomenon 1 --runs 1000 --seed 42
./positronium --mode statistical --phenomenon 4 --runs 1000 --seed 42 \
    --beam-energy-ev 20 --theta-min-deg 5 --angle-bins 10
```

Parametry wiązki można dodatkowo kontrolować przez `--bmax-pm` i
`--matching-radius-pm`; gdy ich nie podano, program dobiera je z energii i
akceptancji kątowej. Sama energia wiązki `--beam-energy-ev` jest domyślnie
monochromatyczna (jak wymaga porównania z formułą Rutherforda przy jednym
\(K_{CM}\)); `--beam-energy-sigma-ev` (domyślnie `0`) opcjonalnie włącza
próbkowanie \(K_{CM}\) każdego zdarzenia z rozkładu Gaussa o tym
odchyleniu standardowym, modelując skończoną rozdzielczość energetyczną
realnej wiązki — patrz audyt warunków startowych wyżej.

Pięć opcji sterują samą fizyką i kosztem eksperymentów związanych:

| Opcja | Domyślnie | Znaczenie |
| --- | --- | --- |
| `--zpf`, `--zpf-band` | `0` (wyłączone) | **Eksperyment, nie część modelu.** Klasyczne pole punktu zerowego elektrodynamiki stochastycznej: losowe fale płaskie o widmie \(\rho(\omega)=\hbar\omega^3/2\pi^2c^3\), 64 mody o równej energii, orientacje i fazy z ziarna `--seed`. `--zpf` skaluje **amplitudę** (1 = poziom fizyczny, moc pochłaniana rośnie jak kwadrat), `--zpf-band lo,hi` ustala pasmo w jednostkach częstości orbitalnej pary (domyślnie `0.3,3`). To jest fluktuacyjna połowa pary fluktuacja–dyssypacja; dyssypacyjną, czyli reakcję promieniowania, model ma od zawsze. Wchodzi w te same trzy miejsca co pole jednorodne, ale próbkowane osobno dla każdej cząstki, bo zależy od położenia i czasu. **Nie odtwarza stanu podstawowego SED — patrz niżej.** |
| `--external-field` | brak (pytanie na starcie) | Jednorodne zewnętrzne pole magnetyczne w mikroteslach; `0` wyłącza. Orientacja jest losowana izotropowo z ziarna `--seed`, więc odtwarza się razem z resztą przebiegu, i jest wypisywana na starcie. Gdy opcji nie podano, a przebieg jest interaktywny, program pyta o to **przed wszystkimi pozostałymi pytaniami** i oferuje 50 µT (skala pola ziemskiego). Przebieg wsadowy z podanym `--mode` i `--phenomenon` nigdy nie pyta i domyślnie nie ma pola. Pole wchodzi w sumę sił chwilowych, w sumę sił retardowanych oraz w pole lokalne widziane przez obie cząstki, przez co obejmuje precesję Thomasa-BMT. Przy 50 µT tempo cyklotronowe \(eB/m\) wynosi 8,8·10⁶ rad/s wobec tempa orbitalnego rzędu 3·10¹⁵ rad/s, więc orbita pozostaje nietknięta, a widocznym kanałem jest precesja dipoli — około 3·10⁻⁴ rad w ciągu 35 ps kolapsu. |
| `--pair` | `electron,positron` | Para cząstek, którą całkuje przebieg, podana jako `pierwsza,druga`. Dostępne gatunki: `electron`, `positron`, `muon`, `antimuon`, `proton`, `antiproton`. Para musi być przyciągająca i nieść przeciwne ładunki elementarne, inaczej opcja jest odrzucana. Wybrana para jest wypisywana na starcie wraz z masą zredukowaną, promieniem Bohra pary i energią wiązania. Honoruje ją także `./positronium_validation`. |
| `--radiation-reaction` | `stochastic` | Model reakcji promieniowania ładunku: `disabled`, `coherent` (Abraham-Lorentz na dipolu elektrycznym pary), `individual` (Landau-Lifszyc zredukowanego rzędu, osobno dla każdej cząstki), `automatic` (mieszanka obu) albo `stochastic` (domyślny od tego miejsca w historii projektu — kwantowane, Poissonowskie kopnięcia fotonowe zamiast ciągłego hamowania, patrz niżej). **Każda liczba czasu kolapsu cytowana wcześniej w tym README (rzędu 36-40 ps dla e⁺e⁻) została zmierzona pod `individual`, nie pod obecnym domyślnym `stochastic`** — żeby je odtworzyć, trzeba dziś podać `--radiation-reaction individual` jawnie. Przy `disabled` żaden kanał nie odbiera energii orbitalnej, więc klasyczna inspirala nie zachodzi i eksperymenty 1/2 zgłaszają brak zaniku. Wybrany model jest wypisywany na starcie. |
| `--beam-energy-sigma-ev` | `0` (wyłączone) | **Tylko eksperymenty 3, 4.** Odchylenie standardowe rozkładu Gaussa, z którego próbkowana jest energia środka masy \(K_{CM}\) każdego zdarzenia wiązki, wokół `--beam-energy-ev`; `0` zachowuje dotychczasową, monochromatyczną wiązkę bit w bit. Próbkowanie odrzuca wyniki \(\le 0\) (do 1000 prób, jak `sampleKinetic` eksperymentu 5), a widmo strat energii liczy się względem faktycznie wylosowanej energii zdarzenia, nie ustalonej średniej. Modeluje skończoną rozdzielczość energetyczną realnej wiązki kosztem rozmycia porównania z formułą Rutherforda, która jest zdefiniowana przy jednym \(K_{CM}\). |

Zasięg `--pair` nie jest jednakowy dla wszystkich eksperymentów.
Eksperymenty 1 i 2 mierzą klasyczną inspiralę CREM — ta część jest ogólna dla
par i liczy się poprawnie dla każdej z nich; dla p+e⁻ daje średni czas kolapsu
17,8 ps przy podręcznikowej wartości klasycznej około 16 ps. Raportowana jest
jednak obok mierzonych czasów życia pozytonium (widmo czasu anihilacji,
odniesienie \(\tau_{\rm exp}\), etykiety para/orto), a te są danymi e⁺e⁻:
proton i elektron nie anihilują, tylko tworzą stabilny wodór, a mionium
i protonium anihilują kanałami, których ten model nie niesie. Dlatego
eksperymenty 1 i 2 są odrzucane dla par innych niż `electron,positron`.
Eksperymenty 3, 4 i 5 nie mają odniesienia anihilacyjnego w części pomiarowej
i działają dla dowolnej pary przyciągającej — z zastrzeżeniem, że eksperyment 5
nadal etykietuje wychwycone stany jako „Para-/Ortho-Positronium" i zapisuje
panele `collapse_time_distribution_*` z odniesieniem τ_exp z danych
pozytonium, co dla innej pary jest mylące.

### Klasyczna gęstość kontaktowa wzdłuż spirali

Kwantowe tempo anihilacji jest rządzone gęstością w zerze,
\(\Gamma_{2\gamma}=4\pi r_e^2c\,|\psi(0)|^2\) przy
\(|\psi(0)|^2=1/(\pi a_{Ps}^3)=2{,}685\cdot10^{29}\,\mathrm m^{-3}\), co daje
\(124{,}5\) ps wobec zmierzonych \(125{,}14\) ps. Naturalne jest pytanie, czy
model klasyczny dostarcza odpowiednika tej gęstości. **Zmierzone: dostarcza go
co do rzędu, ale nie tam, gdzie trzeba.**

Histogram separacji ważony krokiem czasu, zebrany wzdłuż sześciu pełnych
spiral, daje gęstość objętościową \(\rho(r)\) uśrednioną po czasie —
klasyczny odpowiednik \(|\psi(r)|^2\). Zajmuje ona 30 powłok logarytmicznych
od \(5{,}96\) pm do \(168\) pm, z profilem bliskim \(\rho\propto r^{-1{,}8}\).

Na wewnętrznej krawędzi \(\rho=6{,}36\cdot10^{30}\,\mathrm m^{-3}\), czyli
**23,7 raza** \(|\psi(0)|^2\). Rząd wielkości się więc zgadza — spirala
faktycznie zagęszcza parę do właściwej skali.

Poniżej krawędzi \(\rho\) wynosi jednak **dokładnie zero**, i to jest różnica
jakościowa, nie ilościowa. Każda orbita keplerowska ma barierę odśrodkową
i twardy punkt zwrotny; zespół ma krawędź w najmniejszym osiągniętym peryapsis.
Funkcja falowa stanu \(1s\) jest w zerze skończona. Żadna ekstrapolacja nie
łączy twardego zera ze skończoną wartością.

Skale mówią resztę i warto je uporządkować, bo dwie różne rzeczy bywają mylone.
Trajektorie zatrzymuje **granica zderzenia** \(0{,}01a_0=529\) fm, czyli
zadeklarowana rozdzielczość przestrzenna modelu, równa rozmiarowi chmury
ładunku. Bariera dipol-dipol \(193\) fm leży \(2{,}74\) raza **poniżej** tego
punktu i nigdy nie jest osiągana. Krawędź gęstości \(5{,}96\) pm leży wyżej od
obu, bo bieg kolapsu jest ucinany na peryapsis \(0{,}1a_0\).

Rozróżnienie ma znaczenie, bo pierwsza z tych granic jest **parametrem**,
a druga nie. Granicę zderzenia da się obniżyć, zmniejszając
`chargeCloudRestRadius`; zyskuje się wtedy co najwyżej czynnik \(2{,}74\), do
bariery. Barierę wyznacza natomiast tożsamość

\[r^\*=\frac{g}{2}\cdot\frac{\bar\lambda_C}{2}=193{,}30\ \mathrm{fm},
\qquad \frac{r^\*}{r_e}=\frac{g/2}{2\alpha}=68{,}60,\]

czyli **połowa zredukowanej długości Comptona**, a jej odległość od promienia
anihilacji to \(1/(2\alpha)=68{,}52\). To nie jest wybór regularyzacji ani nic
nastawialnego — to stała struktury subtelnej. Bariera stoi dokładnie tam, gdzie
lokalizacja elektronu wymaga pędów rzędu \(mc\), czyli na granicy, poniżej
której klasyczna elektrodynamika punktowej cząstki przestaje obowiązywać;
anihilacja zachodzi \(137\) razy głębiej w obszarze kwantowym.

Sprawdzone, jak głęboko model da się zepchnąć. Bieg kolapsu kończy się na
peryapsis \(10\times\)`chargeCloudRestRadius`; obniżenie tej krotności do
\(1\) działa i schodzi dziesięciokrotnie głębiej, do \(596\) fm, przy 3 z 4
ukończonych trajektoriach i czasie kolapsu \(39{,}5\) ps wobec \(35{,}6\) ps.
Przy krotności \(0{,}5\), czyli \(265\) fm, **wszystkie cztery trajektorie
kończą się awarią numeryczną** — i to jest istotne, bo \(265\) fm leży wciąż
**powyżej** bariery \(193\) fm. Ścianą jest więc integrator, a nie fizyka
bariery; do samej bariery model nie dochodzi.

Awarie przy krotności \(0{,}5\) zostały zbadane osobno, bo nasuwały podejrzenie,
że `double` przestaje wystarczać. **Nie przestaje.** Estymator liczy stratę
energii na obieg jako różnicę `measured = raw - background`, więc kasowanie się
cyfr jest tam realne, ale zmierzone: przez cały bieg zabiera najwyżej
\(1{,}13\) cyfry znaczącej, a przez większość biegu \(0{,}2\)–\(0{,}5\).
Wobec \(15{,}95\) cyfr `double` zostaje margines około \(14{,}8\) cyfr;
wyczerpanie go wymagałoby stosunku `raw/measured` rzędu \(10^{16}\), podczas
gdy największy zaobserwowany to \(13{,}5\). Szersze typy nic by tu nie dały.

Zbadane też, jak głęboko da się zejść po obniżeniu **obu** sprzężonych granic
naraz, czyli krotności końcowego zbliżenia i zasięgu wewnętrznego całkowania.
Odpowiedź: **nie da się poniżej mniej więcej \(529\) fm.** Przy zasięgu
\(53\) fm, \(5{,}3\) fm i \(1\) fm wzorzec awarii jest **identyczny** —
te same dwie trajektorie z trzech giną tak samo, niezależnie od tego, jak
głęboko się prosi. Trajektorie umierają więc w tym samym miejscu, a żądanie
\(1\) fm nie przybliża do niego ani o krok. Skala \(1\) fm leży zresztą
poniżej promienia regularyzacji dipola \(68\) fm, poniżej bariery
Comptona \(193\) fm i poniżej samego \(r_e=2{,}82\) fm, więc model
produkowałby tam liczby w reżimie, w którym klasyczna elektrodynamika
punktowej cząstki na pewno nie obowiązuje.

Na samej podłodze, przy zasięgu \(529\) fm, ukończenie spada do 2 z 3
trajektorii, a gęstość na krawędzi \(596\) fm wynosi
\(8{,}2\cdot10^{31}\,\mathrm m^{-3}\), czyli \(307\times|\psi(0)|^2\).

Podłoga wyglądała na **ostrą** — \(520\) fm dawało 2 z 3 trajektorii,
\(500\) fm 1 z 3, a \(470\) fm i niżej żadnej — ale ten skan zrobiono na
jednym ziarnie i wniosek był przez to zbyt mocny. Powtórzony na ośmiu różnych
ziarnach **cel \(193\) fm osiąga 2 z 8 trajektorii**. Nie ma więc ściany,
tylko malejący z głębokością **odsetek dochodzących**: bariera Comptona
\(193\) fm jest osiągalna dla mniej więcej ćwiartki trajektorii, a które
konkretnie, zależy od mimośrodu wylosowanej orbity.

Głębsze zejście **nie zmienia mierzonego czasu kolapsu**, co potwierdza
założenie o ucinaniu ostatniego odcinka: dla ziarna 7 wychodzi \(40{,}488\) ps
przy \(5{,}29\) pm i \(40{,}492\) ps przy \(193\) fm, dla ziarna 23
odpowiednio \(19{,}804\) i \(19{,}806\) ps — w obu wypadkach różnica
\(+0{,}01\%\), przy zejściu \(27\) razy głębiej.

Jej charakter zbadano trzema gałkami. Zwiększenie głębokości podziału kroku
z 12 na 20 **nie zmienia nic** — wzorzec awarii pozostaje ten sam.
Zacieśnienie tolerancji z \(10^{-5}\) do \(10^{-7}\) usuwa awarie przy tym
samym budżecie, ale trajektorie są wtedy cenzurowane zegarem zamiast kończyć;
przy budżecie podniesionym czterokrotnie awaria wraca. Podłoga jest więc
**wspólnym ograniczeniem dokładności i kosztu** integratora adaptacyjnego, a nie
ścianą strukturalną ani fizyczną — da się ją przesuwać, ale coraz drożej i o
coraz mniej.

**Nieaktualne po naprawie siatki historii retardowanej** (`d164f69`, patrz
"Podłoga zasięgu CREM: bariera Comptona zamiast granicy zderzenia" wyżej).
Cała powyższa podłoga \(529\) fm okazała się artefaktem sprzężenia siatki
historii z krokiem, nie ograniczeniem dokładności/kosztu integratora
adaptacyjnego, za jakie ją tu wzięto. Po naprawie awarie na sześciu ziarnach
skupiają się w paśmie \(47\)–\(495\) fm (mediana \(155\) fm), rząd wielkości
głębiej i wokół bariery \(193\) fm, nie \(529\) fm. `crem_collapse.hpp` celuje
teraz w `comptonBarrierRadius` zamiast w `collisionBoundaryRadius`.

Rozważana w związku z tym wymiana silnika na schemat symplektyczny **nie ma
podstaw** i pomiar to przesądził. Obecny krok jest już strukturą
kick-drift-kick, a dla członu zależnego wyłącznie od położenia taka struktura
jest symplektyczna. Zmierzony błąd energii zachowawczej przy wyłączonej reakcji
promieniowania **oscyluje i nie narasta**: amplituda w kolejnych ćwiartkach
przebiegu wynosi \(4{,}84\), \(4{,}66\), \(4{,}62\) i
\(4{,}83\cdot10^{-5}\), a wartość przechodzi przez zero. Dryfu sekularnego
nie ma, więc przepisanie na splitting nic by tu nie kupiło.

Co pozostaje prawdziwe, to **stosunek skal**: ta ograniczona oscylacja
\(4{,}8\cdot10^{-5}\) jest osiemnastokrotnie większa od mierzonego sygnału
fizycznego, czyli straty energii na obieg \(2{,}6\cdot10^{-6}\). To jest
powód, dla którego estymator odejmuje przebieg tła — kasuje zależny od fazy,
ograniczony błąd schematu, a nie dryf. Nasuwało to receptę **wyższego rzędu**, a nie
symplektyczności. Została zaimplementowana i **zmierzona jako gorsza**.

Złożenie Yoshidy czwartego rzędu jest dostępne pod `--integrator-order 4`
i domyślnie wyłączone. Porównane z krokiem bazowym przy tej samej tolerancji
wypada źle na obu osiach naraz: przy \(10^{-7}\) daje błąd
\(1{,}62\cdot10^{-5}\) w \(16{,}9\) s wobec \(5{,}04\cdot10^{-6}\)
w \(5{,}68\) s dla rzędu drugiego, czyli jest \(3{,}2\) raza mniej dokładne
i \(3{,}0\) raza wolniejsze. Powodów jest kilka: krok środkowy biegnie wstecz
w czasie, wszystkie trzy podkroki czytają historię zamrożoną na początku kroku,
a kontroler adaptacyjny i tak wyrównuje błąd do tolerancji, więc wyższy rząd
kupowałby najwyżej szybkość — i nie kupuje nawet jej.

Pomiar ujawnił przy tym próg, którego żadna z dwóch hipotez nie trafiała. Błąd
nie schodzi poniżej około \(5\cdot10^{-6}\) niezależnie od tolerancji: przy
\(10^{-7}\) wynosi \(5{,}04\cdot10^{-6}\), a przy \(10^{-9}\) **rośnie**
do \(6{,}78\cdot10^{-6}\) przy czterokrotnie większej liczbie kroków.

Wyglądało to na akumulację zaokrągleń i zostało tak sprawdzone: `Vec3`
przestawiono na `long double`, co potwierdzono przez \(48\)-bajtowy rozmiar
struktury zamiast \(24\). **Zaokrąglenia są wykluczone.** Przebieg daje
wtedy dokładnie tę samą liczbę kroków — \(25\,434\) i \(106\,902\) — oraz
błąd zgodny co do czterech cyfr, \(5{,}043\cdot10^{-6}\) i
\(6{,}778\cdot10^{-6}\). Podwojenie precyzji mantysy nie zmienia niczego, więc
próg jest deterministyczny, a nie szumowy.

Czym jest, pozostaje nierozstrzygnięte. Wychylenie energii mechanicznej jest co
do bitu lustrzane wobec rezerwuaru pola związanego — suma wychodzi dokładnie
zero — ale jest to **ta sama tożsamość**, o której mowa wyżej, więc nie odróżnia
fizycznej wymiany z polem bliskim od deterministycznego błędu schematu
wchłoniętego przez resztę. Praktyczny wniosek jest natomiast jednoznaczny: próg
leży niecałe dwa razy powyżej sygnału \(2{,}6\cdot10^{-6}\) na obieg i **nie
usuwa go ani rząd, ani tolerancja, ani precyzja**, więc odejmowanie przebiegu
tła jest w tym schemacie konieczne.

Prawdziwe przyczyny są dwie i żadna nie jest precyzją. Po pierwsze **granice są
sprzężone**: wewnętrzne całkowanie ma `terminalSeparation` równe
`chargeCloudRestRadius`, czyli \(529\) fm, więc żądanie peryapsis \(265\) fm
prosi o coś, czego niższy poziom nie potrafi dostarczyć — obniżenie samej
krotności końcowego zbliżenia nie wystarcza. Po drugie tam, gdzie całkowanie
faktycznie rozbiega, przekroczenie wynosi \(1{,}66\cdot10^{14}\) energii
wiązania na obieg. To rozbieg zredukowanej reakcji Landaua-Lifshitza przy
zacieśniającej się orbicie, a nie utrata cyfr — czternastu rzędów żadna
precyzja nie odzyska. Wyłapuje go istniejąca bramka odrzucająca stratę powyżej
50% energii wiązania na obieg.

Ważniejszy jest kształt profilu w głębi. Gęstość **nie rośnie monotonicznie do
środka**: osiąga maksimum \(1{,}58\cdot10^{32}\,\mathrm m^{-3}\) przy
\(750\) fm, czyli \(589\times|\psi(0)|^2\), po czym **spada** do
\(7{,}2\cdot10^{31}\) na krawędzi \(596\) fm. Powód jest elementarny: przy
największym zbliżeniu para porusza się najszybciej, więc spędza tam najmniej
czasu — ułamek czasu w najgłębszej powłoce wynosi \(2{,}2\cdot10^{-5}\).

To domyka sprawę mocniej niż samo twarde zero. Klasyczna gęstość ma
**maksimum przy skończonym promieniu i maleje do wewnątrz**, podczas gdy
\(|\psi(0)|^2\) jest gęstością dokładnie w zerze i jest tam największa.
Ekstrapolacja \(r^{-1{,}8}\), którą sugerował profil zewnętrzny, w głębi po
prostu nie obowiązuje.

Wniosek: **czasu kolapsu CREM nie da się połączyć z czasem życia para-Ps przez
gęstość kontaktową.** Nawet po obniżeniu granicy zderzenia do samej bariery
model zatrzyma się \(68{,}5\,r_e\) od skali procesu, a tego czynnika nie da się
zmniejszyć bez porzucenia klasycznego opisu — jest nim \(1/(2\alpha)\).

Dodatkowo rozkłady mają inny kształt: anihilacja daje rozkład wykładniczy,
dla którego \(\sigma/\overline{t}=1\), a zmierzony rozkład czasu kolapsu ma
\(\sigma/\overline{t}=0{,}38\)–\(0{,}46\). Model nie rozróżnia też para od
orto (\(1{,}00\times\) wobec \(1135\times\) w eksperymencie), więc każde
dopasowanie do jednej z tych liczb rozmija się z drugą o trzy rzędy.

Wynik eksperymentu z polem punktu zerowego (`--zpf`) jest **negatywny i warto
go znać przed uruchomieniem**. Pasmo **podąża za orbitą**: jest zdefiniowane
w krotnościach oskulacyjnej częstości orbitalnej pary i przeliczane przy każdym
wywołaniu siły, więc zostaje w rezonansie przez całą spiralę. Wymaga to niesienia
fazy w stanie (`State::zeroPointPhase`, całka z częstości orbitalnej), bo faza
zapisana jako \(-\omega(t)t\) skakałaby przy każdej zmianie pasma i wstrzykiwała
energię z niczego. Amplituda rośnie jak \(\omega_{orb}^2\), co wynika z widma.

To pozwala oddzielić efekt rezonansowy od efektu obcięcia i **rozstrzyga sprawę
na niekorzyść mechanizmu**. Przy paśmie podążającym, wąskim i pozostającym
w rezonansie, wpływ na czas kolapsu jest znikomy: mediana 36,12 ps przy paśmie
\([0{,}3;3]\) i 36,123 ps przy \([0{,}5;2]\), wobec 35,998 ps bez pola — czyli
około \(+0{,}3\%\), a obie wartości zgadzają się ze sobą do czterech cyfr. Wynik
jest zbieżny po liczbie modów: 36,69 / 36,12 / 36,68 ps przy 32 / 64 / 256
modach, więc nie jest artefaktem próbkowania.

Duże wydłużenia pojawiają się dopiero po rozszerzeniu górnej krawędzi daleko
poza rezonans — 83,9 ps przy \([0{,}3;10]\) i brak zaniku przy \([0{,}3;30]\) —
czyli pochodzą od modów **daleko od rezonansu**, o dużej amplitudzie, a nie od
absorpcji rezonansowej, na której opiera się równowaga SED. To jest odwrotność
tego, czego mechanizm wymaga. Mechanizm działa: energia wpływa do orbity,
orbita się rozszerza, a czas kolapsu rośnie monotonicznie z górną krawędzią
pasma — 43,0 ps bez pola, 43,2 przy `0.3,3`, 45,2 przy `0.3,10`, 57,3 przy
`0.3,30` i 152,5 przy `0.3,100`. Przy paśmie USTALONYM na starcie orbita ucieka polu: zmierzony okres spada z 0,327 fs do 0,0031 fs, czyli o czynnik 105.

Przy `0.3,300` estymator zgłasza brak zaniku — ale to **nie jest równowaga**.
Diagnostyka pokazuje promień końcowy 219,5 pm wobec 123,7 pm bez pola i zakres
101,6–225,8 pm: orbita jest przepompowywana i się rozszerza, co utrzymane
prowadziłoby do jonizacji. Przejście jest więc z „zapada się" wprost w
„rozszerza się", bez stanu stacjonarnego pomiędzy.

Fizycznie widmo pola punktu zerowego nie ma górnego odcięcia, więc pasmo jest
wyłącznie obcięciem numerycznym, a wynik od niego silnie zależy — to znaczy, że
ta implementacja nie daje odpowiedzi fizycznej. Poprawne odtworzenie równowagi
SED wymaga samouzgodnionej odpowiedzi cząstki na pełne widmo, a nie dołożenia
losowego pola do gotowych równań ruchu; równowaga jest delikatnym skasowaniem,
nie efektem rzędu wiodącego.

**Sprawdzone wprost, czy to obcięcie po prostu za mało modów, a nie sam
mechanizm.** Powyższy wynik przy `0.3,300` użył 64 modów. Zmierzone na tym
samym paśmie i ziarnie (42) przy 64/256/1024 modach: promień końcowy
232,7/115,1/161,4 pm, zakres [103,0–233,1]/[104,8–150,6]/[104,1–176,0] pm —
**żadnej zbieżności z liczbą modów**, bo każda liczba modów losuje zupełnie
nową, niepowiązaną realizację pola (ten sam seed nie zagęszcza poprzedniej
próby, tylko losuje inną). Powtórzone na dwóch dalszych ziarnach przy 256
modach (najtańszy koszt spośród testowanych): ziarno 7 daje promień końcowy
165,6 pm, zakres [105,7–165,6] pm; ziarno 1 daje 98,2 pm, zakres
[75,3–134,4] pm — ten ostatni przypadek wygląda na mniejszą ucieczkę, ale to
złudzenie: przewidywane peryapsis tej konkretnej próbkowanej orbity to
74,7 pm (wobec ~105 pm dla pozostałych ziaren), więc dolna granica to
naturalne przejście przez peryapsis Keplera, nie efekt pola.

**Wszystkie pięć prób (3 liczby modów × ziarno 42, 2 dodatkowe ziarna przy
256 modach) kończą się identycznie: `trajectory: FAIL`** — promień
przekracza apoapsis Keplera startowej orbity o wygodny margines w każdym
przypadku, bez wyjątku. Szczegółowy kształt (monotoniczna ucieczka jak przy
ziarnach 42/7, czy nie-monotoniczny wychył jak przy ziarnie 1) zmienia się
z każdą realizacją, ale kierunek — **żadnej stabilizacji, żadnego
ograniczonego stanu równowagi** — jest identyczny za każdym razem. To
domyka sprawę mocniej niż pojedynczy pomiar: gdyby wynik zależał tylko od
niedopróbkowania wysokich częstości przy 64 modach, zwiększenie rozdzielczości
16-krotnie powinno było dać zbieżność do jednej liczby. Nie dało — bo
problem nie jest w rozdzielczości próbkowania, tylko w samej strukturze
(skończona suma modów, obcięta arbitralnie w paśmie), dokładnie jak
zdiagnozowano wyżej.

Nie przetestowano wprost pasma sięgającego częstości Comptona
`ω_C=m_ec²/ħ≈7,76·10²⁰` rad/s — przy startowej częstości orbitalnej
`ω_orb≈2,07·10¹⁶` rad/s dla p-Ps to krotność `≈37 600`, ponad sto razy
szersza niż już przetestowane `0.3,300`, a koszt jednego biegu przy `0.3,300`
i 1024 modach to już rząd godziny. Ekstrapolacja z ustalonego trendu
(szersze pasmo = więcej energii z modów nierezonansowych = mocniejsza
ucieczka, nigdy stabilizacja) czyni bieg przy pełnym odcięciu Comptona
bardzo drogim testem obarczonym niemal pewnym wynikiem — nie uzasadnia to
kosztu. Wątek SED/`--zpf` jako sposobu na wydłużenie symulowanego czasu
kolapsu uznaję za zamknięty.

Dwie własności modeli reakcji, które warto znać przed użyciem:

- **`individual` zawiera człon wzajemny i odtwarza koherentne tempo
  dipolowe.** Rozkładając moc koherentną
  `|q₁a₁+q₂a₂|² = |q₁a₁|² + |q₂a₂|² + 2q₁q₂(a₁·a₂)`, samosiła per-cząstkowa
  daje tylko dwa pierwsze człony; dla e⁺e⁻ człon interferencyjny równa się ich
  sumie, więc jego brak kosztował dokładnie czynnik 2 w stracie energii
  orbitalnej. Po dołożeniu `F_i = (q_i q_j/(6πε₀c³))·ȧ_j` panel
  `N_1_4_radiated_power_vs_larmor.pdf` mierzy `⟨P_CREM/P_Larmor⟩ = 0,998`,
  panel `N_1_3_collapse_time_vs_theory.pdf` daje `⟨t_CREM/t_klasyczne⟩ = 1,02`,
  a residuum bilansu `|dE_LL-vs-strumień|/E_rad` spada z 0,501 do **0,003**.
  Czas kolapsu skrócił się przy tym dokładnie dwukrotnie.
- **`coherent` wykazuje rozbieg Abrahama-Lorentza i nie nadaje się do pomiaru
  sekularnego.** Siła jest budowana z trzeciej pochodnej momentu dipolowego,
  co dopuszcza rozwiązania samoprzyspieszające: część zmierzonych orbit
  *zyskuje* energię zamiast tracić, a sporadycznie stencil dzielony przez
  `2h³` zwraca wartość niefizyczną. Estymator odrzuca teraz pomiar implikujący
  stratę powyżej 50% energii wiązania na jedną orbitę i zgłasza go jako awarię
  numeryczną, zamiast wcielać go do stanu oskulacyjnego.
| `--crem-wallclock-budget-s` | `20` | Skończony dodatni budżet zegarowy na jedno zdarzenie w eksperymentach 1/2. Po jego wyczerpaniu trajektoria jest cenzurowana prawostronnie. Wartość domyślna daje niską frakcję ukończenia przy `N = 1000`; do wniosków ilościowych trzeba ją podnieść tak, aby kompletacja zbliżyła się do 100%. |

`--no-gui` wykonuje serię i wypisuje podsumowanie bez
otwierania okna ROOT, ale nadal renderuje i zapisuje pliki PDF, co jest
przydatne w obliczeniach wsadowych:

```bash
./positronium --mode statistical --phenomenon 2 --runs 1000 \
    --seed 42 --no-gui
```

Tryb diagnostyczny nie otwiera okna i wypisuje zakres odległości, bilans
energii oraz dryf \(\mathbf P_N\) i \(\mathbf J_N\):

```bash
./positronium --diagnose --phenomenon 4 --seed 42
```

Dla referencyjnego ziarna 42 wszystkie cztery scenariusze przechodzą kontrolę
kierunku trajektorii. Po włączeniu Darwina względny dryf bilansu wynosi około
0,03% dla rozproszenia, 0,32% dla bezpośredniego zderzenia oraz 4,6% dla
końcowej, ciasnej fazy obu trajektorii związanych. Te wartości są testem
numerycznym konkretnej konfiguracji, a nie oszacowaniem niepewności fizycznej
modelu. Wyniki przeznaczone do analizy ilościowej wymagają osobnego badania
zbieżności względem kroku czasu.

Przyciski `STOP`/`START` sterują animacją, a `EXIT` zamyka program.

## Ograniczenia — efekty nieuwzględniane

- mechanika kwantowa w dynamice orbitalnej, funkcja falowa, zasada Pauliego,
  splątanie i dynamiczne tworzenie stanów para-/ortopozytonium;
- kwantowy czas anihilacji nie jest wyprowadzany z dynamiki klasycznej;
  Statistical 1/2 raportuje odrębny, operacyjny czas kolapsu CREM mierzony
  przez całkowanie mechaniczne z uśrednianiem po orbitach, natomiast
  zewnętrzne czasy życia służą wyłącznie do porównania;
- odpowiedź detektora, oddziaływania pozytonium w materiale, pick-off,
  quenching, poprawki radiacyjne i rzadsze kanały rozpadu;
- relatywistyczne amplitudy QED niezwiązanej pary, w tym rozpraszanie Bhabhy
  i przekrój anihilacji e⁺e⁻;
- kompletna dynamika pola Maxwella w torze produkcyjnym. Same **wzajemnie
  opóźnione pola Liénarda–Wiecherta są uwzględnione** i to one napędzają
  trajektorię (`retardedExternalForces`, patrz sekcja *Retardowane pola
  Liénarda–Wiecherta*); nieobecne jest rozwiązywanie równań Maxwella na
  siatce, które kompiluje się wyłącznie do `positronium_validation` jako
  kontrola krzyżowa i nie bierze udziału w produkcyjnym całkowaniu;
- ograniczenie do rzędu \(v^2/c^2\) **nie obowiązuje** w torze produkcyjnym:
  siła między ładunkami liczona jest z pełnego pola Liénarda–Wiecherta,
  dokładnego we wszystkich rzędach \(v/c\) wraz z retardacją. Przybliżenie
  Darwina, które jest rzędu \(v^2/c^2\), występuje wyłącznie w raportowanej
  diagnostyce energii oraz jako zalążek rekonstrukcji historii przyczynowej,
  natychmiast zastępowany dwiema iteracjami Picarda na siłach retardowanych;
- odrzut orbitalny od promieniowania dipola magnetycznego. Koherentna moc M1
  jest liczona z \(|\ddot{\mathbf m}_1+\ddot{\mathbf m}_2|^2\), a jej energia
  księgowana w `dipoleConstraintEnergy`, czyli poza
  orbitą. Walidator sprawdza oba skrajne przypadki interferencji: dwa zgodne
  źródła dają \(4P_1\), a dwa przeciwne — zerową moc. Odrzut M1 pozostaje
  poza ruchem orbitalnym; pominięcie jest świadome;
- kwantowa struktura subtelna i kwantowe elementy macierzowe sprzężenia
  spin–orbita. Klasyczne odpowiedniki **są uwzględnione**: relatywistyczna
  precesja BMT z mierzonym \(g\), precesja Thomasa i transformacja tensora
  dipolowego między układami odniesienia; człon \(\boldsymbol\beta\times\mathbf E\)
  w polu efektywnym BMT to właśnie ruchowe pole magnetyczne, czyli klasyczne
  sprzężenie spin–orbita;
- skończony rozmiar cząstek i struktura krótkiego zasięgu — zastępuje je próg
  \(10^{-14}\,\mathrm m\);
- zderzenia z materią, grawitacja i wpływ ośrodka (jednorodne zewnętrzne pole
  magnetyczne jest natomiast dostępne — patrz `--external-field` — ale
  zewnętrzne pole elektryczne i pola niejednorodne nadal nie);
- energia orbitalnego pola magnetycznego w prezentowanym bilansie energii;
- wyższe multipole promieniowania jako kanały reakcji: moc dipola
  magnetycznego i kwadrupola elektrycznego jest liczona, ale wyłącznie jako
  wejście bramki `dominance` ważącej model koherentny E1 przeciw
  Landauowi-Lifszycowi przy `--radiation-reaction automatic`. Do bilansu
  energii nie wchodzi żadna z nich, a dalsze człony rozwinięcia nie są
  liczone w ogóle.

**Audyt fizyki (przegląd trzech niezależnych torów: pola retardowane, reakcja
promieniowania, precesja spinu/Darwin/dipol-dipol) znalazł i naprawił jeden
błąd numeryczny**: `electricQuadrupoleRadiatedPower()` miał współczynnik SI
`1/(180πε₀c⁵)` zamiast poprawnego `1/(720πε₀c⁵)` (brakujący czynnik `4π`,
jaki niosą wszystkie sąsiednie stałe typu kulombowskiego w tym pliku) —
potwierdzone trzema niezależnymi przeliczeniami (Gaussa→SI tą samą metodą co
już zwalidowane człony E1/M1). Wpływ produkcyjny znikomy: ta moc wchodzi
wyłącznie do nasyconej bramki `dominance` (próg 10-20, zmierzone wartości
~3,5·10¹³ dla e⁺e⁻ i ~3,3·10³ dla p+e⁻ — 4× korekta niczego nie przełącza),
nigdy do bilansu energii. Zweryfikowane po poprawce: `positronium_validation`
33/33 bez zmian.

**Ten sam audyt znalazł też nierozstrzygniętą jeszcze rozbieżność w precesji
spinu**, nienaprawioną celowo (wymaga osobnej derywacji, nie zgadywania,
która strona jest błędna). Kod ma dwie implementacje: `thomasBmtEffectiveField`
(dokładnie odtwarza wzór Jacksona, ale nie jest wołana z produkcyjnego kroku)
i `advanceCovariantBmt` (kowariantne równanie tensorowe RK4, faktyczna ścieżka
produkcyjna). Zmierzone bezpośrednio (niezależna reimplementacja w Pythonie,
zbieżna po kroku czasowym): przy `g=2` obie zgadzają się dokładnie dla
każdego β, ale przy `g≠2` rozbieżność rośnie z anomalią i prędkością —
5,8·10⁻⁸ dla e⁺e⁻ przy orbitalnym β≈0,007 (nieistotne), ale 3,0%/13,8%/45,4%
dla protonu (`g=5,5857`) przy β=0,3/0,6/0,9. Istniejący test `covariant BMT`
w `positronium_validation` sprawdza `advanceCovariantBmt` tylko względem
samego siebie pod boostem (samospójność), nigdy względem niezależnie
poprawnej `thomasBmtEffectiveField` — 33/33 nic tu nie gwarantuje. Dla
domyślnej pary e⁺e⁻ przy prędkościach pozytonium to nieistotne; dla ciężkich,
wysoko-anomalnych par pod `--radiation-reaction automatic` przy
relatywistycznych prędkościach to nadal otwarte pytanie.

**Dochodzenie kontynuowane.** Sprawdzone bezpośrednio w kodzie: konwersja
`properDipole`→`state.firstDipole` (`synchronizeCovariantDipoles`,
[positronium.cpp:926-941](positronium.cpp)) idzie przez boost **tensorowy**
(`lorentzBoostDipole`), nie przez samą część przestrzenną czterowektora
(`properDipoleFromFourVector`) — więc porównanie tempa precesji przez odczyt
tensorowy (ten, którego użyto wyżej) jest metodologicznie poprawne i
produkcyjnie-relewantne; sprawdzona odwrotna konwencja daje jeszcze gorszą
zgodność, nawet przy `g=2`. Rozbieżność jest więc realna, nie artefaktem
metody porównania. Prawdopodobne źródło: to ten sam mechanizm co
`covarianceRepresentabilityGap` ("dipole repr gap", `<5·10⁻³`, zmierzone
`0,00386`) — rotacja Wignera przy złożeniu boostu z przyspieszeniem cząstki
— tylko silniej widoczna w dynamice (tempo precesji) niż w statycznym teście
jednego wektora. Dodany diagnostyczny (bez progu pass/fail — przedwczesne,
zanim znane jest źródło) test w `positronium_validation`: `BMT vs eff field`,
mierzący dokładnie tę rozbieżność dla aktywnej pary przy β=0,1 oraz dla
syntetycznej sondy proton/β=0,9, żeby przyszła zmiana kodu nie pogorszyła
tego po cichu. Zmierzone bazowe wartości: `5,8·10⁻⁶` (e⁺e⁻, β=0,1) /
`0,4538` (syntetyczna sonda) — druga liczba jest stała niezależnie od
`--pair`, pierwsza rośnie z anomalią pary (`0,0032` dla p+e⁻ przy tym samym
β=0,1). Naprawa nadal wymaga osobnej derywacji relacji między czterowektorem
spinu a tensorem dipolowym pod przyspieszeniem przy wszystkich rzędach v/c —
nie podjęta tutaj.

**Naprawione.** Wyprowadzone od zera (redukcja kowariantnego równania BMT do
3D, przy prędkości zamrożonej dokładnie tak, jak `applyDipolePrecession`
faktycznie wywołuje precesję — dwa pół-kroki, każdy przy stałym `v`):
zredukowane równanie **dokładnie** odtwarza `thomasBmtEffectiveField`
(zgodność do zera maszynowego, `relerr=0`), potwierdzając, że oba wzory są
fizycznie poprawne z osobna. Rozbieżność w faktycznej implementacji
`advanceCovariantBmt` miała inne źródło: więz `a·u=0` (spin czterowektor
prostopadły do czteroprędkości) jest zachowywany przez ciągłe równanie
kowariantne **tylko wtedy, gdy `u` współewoluuje** razem ze spinem
(`du/dτ=(q/m)F·u`) — a `advanceCovariantBmt` z konstrukcji zamraża `v` (i
`u`) na czas całego wywołania, dokładnie tak jak produkcja go zawsze
wywoływała. Przy zamrożonym `u` więz **nie jest** zachowywany przez samo
równanie (`d(a·u)/dτ=-(q/m)Kc²≠0` dla ogólnego `K=a·(F·u)/c²`), a końcowa
projekcja/renormalizacja (potrzebna, żeby czterowektor zostawić poprawnym)
wstrzykiwała poprawkę tego samego rzędu co sam sygnał fizyczny — dokładnie
proporcjonalną do członu anomalnego `(g/2−1)`, stąd zero przy `g=2` i wzrost
z anomalią/prędkością gdzie indziej.

Naprawa: nowa `advanceThomasBmtDipole`/`properDipolePrecessionRate`
([modules/electrodynamics.hpp](modules/electrodynamics.hpp)) omija
czterowektor całkowicie. Zamiast tego: (1) rozwiń zmierzony dipol
laboratoryjny (`lorentzBoostDipole` — ten sam tensorowy boost, którego
`synchronizeCovariantDipoles` już używa do wyprowadzenia `state.firstDipole`
wszędzie indziej w kodzie) z bieżącego `properDipole`; (2) policz jego
tempo wprost ze wzoru Jacksona, `d(μ_lab)/dt=(q/m)μ_lab×B_eff`
(`thomasBmtEffectiveField`) — czysta rotacja, więc dokładnie zachowuje normę
`μ_lab` bez żadnej renormalizacji; (3) odwzoruj to tempo z powrotem na
`properDipole` przez wyprowadzony w zamkniętej formie odwrotny boost —
boost tensorowy skaluje czynnikiem γ **wyłącznie** składową prostopadłą do
`v`, składowa równoległa przechodzi bez zmian, więc odwrotność to po prostu
podzielenie części prostopadłej przez γ. Zweryfikowane niezależnie w
Pythonie do precyzji maszynowej (`relerr~10⁻¹⁶`) na całym zakresie kroku
czasowego, dla tego samego przypadku wysokiej anomalii/prędkości, który
ujawnił błąd starej formuły.

RK4 po tym torze nie zachowuje dokładnie normy `properDipole` (obcięcie
O(dt⁵), inaczej niż dla `μ_lab`, który jest chroniony samą strukturą
rotacji) — zmierzone bezpośrednio: dryf normy dipola na pełnej trajektorii
e⁺e⁻ skoczył z `4·10⁻¹⁴%` do `2,5·10⁻¹⁰%` bez renormalizacji, przekraczając
próg samotestu `1·10⁻¹²`. Naprawione dodaniem tej samej jawnej
renormalizacji `properDipole`, jakiej zawsze używała stara
`advanceCovariantBmt` — po niej dryf wraca do `2,9·10⁻¹⁴%`, w granicach
poprzedniego poziomu.

Zweryfikowane: `positronium_validation` 33/33 dla e⁺e⁻ i p+e⁻ (test
`role-routing`, jedyny bezpośrednio wołający starą funkcję, przełączony na
nową). Diagnostyka `BMT vs eff field` (dodana wcześniej w tym samym audycie)
spadła z `5,8·10⁻⁶`/`0,4538` do `2,2·10⁻¹⁶`/`1,7·10⁻¹⁶` (e⁺e⁻) i z `0,0032`
do `0`/`1,7·10⁻¹⁶` (p+e⁻) — czyste zaokrąglenie maszynowe, dokładnie
przewidziane przez niezależną weryfikację w Pythonie. `--diagnose` na pięciu
ziarnach (e⁺e⁻ i p+e⁻) daje `trajectory: PASS` z dryfem normy dipola
porównywalnym do stanu sprzed całej tej serii poprawek. `advanceCovariantBmt`
pozostaje w kodzie (oznaczona `[[maybe_unused]]` w zwykłej kompilacji) —
używana już tylko przez własne testy samospójności pod boostem w
`maxwell_validation.hpp`, nie przez produkcję.

**Ostatnie znalezisko tego samego audytu, teraz zmierzone.** Pole
ładunek-ładunek (`lienardWiechertField`) interpoluje historię retardowaną
przez Hermite'a (`interpolatedCharge`, rząd zbieżności ~2, zmierzone
wyżej), ale pole dipolowe (`retardedElectricDipoleField`/
`retardedMagneticDipoleField`, przez `historicalDipoleKinematics`→
`historicalSource`) interpoluje pozycję/prędkość źródła **liniowo** — bez
dopasowania pochodnych. Nie było dotąd dedykowanego testu geometrycznej
dokładności tej ścieżki na zakrzywionej trajektorii (istniejący test
dipolowy sprawdzał tylko stencile pochodnych w dokładnych węzłach historii,
nie samą interpolację między nimi).

Dodany test (`sourceInterpolationErrors`, ta sama syntetyczna trajektoria
kołowa i ten sam rozstaw, co istniejący test Hermite'a, więc liczby są
bezpośrednio porównywalne): błąd `historicalSource` wychodzi **dokładnie**
`0,0012497396` — bit w bit ta sama liczba, co już zmierzony błąd liniowy
`lienardWiechertField`-owej ścieżki (`interpolationCoarse[0]`/
`interpolationFine[0]`, drukowane jako `history linear/H`), zbieżny rzędu
`1,999` (czysto drugiego rzędu, jak oczekiwane dla interpolacji liniowej).
To **nie jest nowe, ukryte źródło błędu** — to dokładnie ten sam,
już scharakteryzowany błąd interpolacji liniowej, tylko liczony inną
ścieżką kodu (`historicalSource` zamiast `linearlyInterpolatedCharge`), i
z góry ograniczony (Hermite dałby ~3× mniej, `0,000417`, ale sam błąd
liniowy jest już mały na tej skali czasu i typowe kroki produkcyjne są
drobniejsze niż testowany rozstaw — `CFL dt` rzędu `10⁻²²` s wobec
testowanego `10⁻¹⁸`/`0,5·10⁻¹⁸` s). Wynik: asymetria jest realna, ale
nieszkodliwa w praktyce — udokumentowana i mierzona trwale (`dipole src
linear` w `positronium_validation`), bez progu pass/fail (jak inne
diagnostyki w tym audycie), żeby przyszła zmiana nie pogorszyła jej po
cichu.

## Literatura pomocnicza

1. NIST, *2022 CODATA Recommended Values of the Fundamental Physical
   Constants*, [JPCRD 53, 030201 (2024)](https://physics.nist.gov/cuu/pdf/JPCRD2022CODATA.pdf).
2. H. Spohn, *The critical manifold of the Lorentz–Dirac equation*,
   [Europhys. Lett. 50, 287 (2000)](https://doi.org/10.1209/epl/i2000-00233-3).
3. C. Bild, D.-A. Deckert i H. Ruhl, *Radiation reaction in classical
   electrodynamics*, [Phys. Rev. D 99, 096001 (2019)](https://doi.org/10.1103/PhysRevD.99.096001).
4. C. G. Darwin, *The Dynamical Motions of Charged Particles*,
   [Philosophical Magazine 39, 537–551 (1920)](https://doi.org/10.1080/14786440508636066).
5. Y. Aharonov i A. Casher, *Topological Quantum Effects for Neutral
   Particles*, [Phys. Rev. Lett. 53, 319–321 (1984)](https://doi.org/10.1103/PhysRevLett.53.319).
6. Particle Data Group, *Positronium: Atomic and Nuclear Properties* (2025),
   [PDG tables](https://pdg.lbl.gov/2025/AtomicNuclearProperties/positronium.html).
7. D. B. Cassidy, *Experimental progress in positronium laser physics*,
   [Rev. Mod. Phys. 95, 021002 (2023)](https://doi.org/10.1103/RevModPhys.95.021002).
8. A. Ore i J. L. Powell, *Three Photon Annihilation of an Electron-Positron
   Pair*, [Phys. Rev. 75, 1696 (1949)](https://doi.org/10.1103/PhysRev.75.1696),
   wraz z algorytmem generatora opisanym w
   [Geant4 Physics Reference Manual](https://geant4.web.cern.ch/documentation/dev/prm_html/PhysicsReferenceManual/decay/OrthoPositronium.html).
9. A. H. Al-Ramadhan i D. W. Gidley, *New precision measurement of the decay
   rate of singlet positronium*,
   [Phys. Rev. Lett. 72, 1632 (1994)](https://doi.org/10.1103/PhysRevLett.72.1632).
10. R. S. Vallery, P. W. Zitzewitz i D. W. Gidley, *Resolution of the
    Orthopositronium-Lifetime Puzzle*,
    [Phys. Rev. Lett. 90, 203402 (2003)](https://doi.org/10.1103/PhysRevLett.90.203402).
11. T. C. Chang, C. M. Tang i W. L. Li, eksperymentalny test kontinuum 3γ,
    [Phys. Lett. B 157, 357 (1985)](https://doi.org/10.1016/0370-2693(85)90380-6).
12. M. E. Peskin i D. V. Schroeder, *An Introduction to Quantum Field
    Theory* (1995), wzór 5.61 (rozpraszanie t-kanałowe dwóch różnych
    fermionów) i wynik e⁺e⁻→μ⁺μ⁻ z rozdziału 5.1 — punkty odniesienia
    użyte do weryfikacji granicy bezmasowej krzywej Bhabhy/Motta.
13. H. J. Bhabha, *The Scattering of Positrons by Electrons with Exchange
    on Dirac's Theory of the Positron*,
    [Proc. R. Soc. A 154, 195 (1936)](https://doi.org/10.1098/rspa.1936.0046).

Literatura uzasadnia użyte elementy elektrodynamiki klasycznej i wejścia
generatora zaniku, ale nie waliduje fenomenologicznego utożsamienia orientacji
klasycznych dipoli ze stanami para-/ortopozytonium w animacji.
