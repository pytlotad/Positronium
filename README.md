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

Nie jest też **brakującym członem interferencyjnym**. Niezależny test dwóch
pełnych orbit wykazał rzecz przeciwną: pełne pole Liénarda–Wiecherta partnera
już niesie wzajemną część radiacyjną, więc dodanie na nim `firstMutual` albo
pełnego \(q_i\dddot{\mathbf d}\) liczy ją drugi raz. Zmierzone reszty
\((\Delta E_{mech}+E_{far}+\Delta E_S)/E_{far}\) wynoszą 0,0088 dla
Coulomb–Darwin + pełna reakcja koherentna, −0,0594 dla pola retardowanego +
tylko self LL i −0,5101 dla błędnego pola retardowanego + pełna reakcja
koherentna. Produkcja używa więc członów wzajemnych wyłącznie z
nieretardowanym polem Coulomba–Darwina; z `retardedExternalForces` dokłada
tylko reakcję własną.

Historyczny rozkład sprzed usunięcia podwójnego liczenia, wzdłuż kolapsu
(e⁺e⁻, ziarno 12345, w eV): strumień retardowany 0,422,
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

Sama część własna `individualLandauLifshitzSelfForces` jest **poprawna**. Sprawdzona
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

Walidacja zawiera dwa deterministyczne profile statystyczne. Mały profil
(`N=4096`) służy do codziennych regresji, a profil publikacyjny (`N=1000000`)
sprawdza te same momenty rozkładu kierunku fotonu E1 i rozkładu wykładniczego
procesu Poissona z ostrzejszymi tolerancjami. Oba używają stałych, wbudowanych
ziaren i kontrolują bitowo dokładne odtworzenie próbki:

```bash
make validation-small
make validation-publication
```

Przenośny profil publikacyjny nie używa `-march=native`, usuwa bezwzględną
ścieżkę źródeł z artefaktu i ustala ziarno kompilatora. Osobne binaria nie
nadpisują szybkiego buildu lokalnego:

```bash
make reproducible
make reproducible-validation
```

Wersje narzędzi użyte dla bieżącego punktu odniesienia są zapisane w
`toolchain-versions.txt`; `make toolchain-info` wypisuje wersje aktualnego
środowiska. Profile AddressSanitizer i UndefinedBehaviorSanitizer buduje się
i uruchamia poleceniem:

```bash
make sanitizers-check
```

ASan ma wyłączone wykrywanie wycieków, ponieważ proces ładuje ROOT; nadal
sprawdza błędy dostępu do pamięci. UBSan zatrzymuje test przy pierwszym
wykrytym niezdefiniowanym zachowaniu. Sanitizatory domyślnie korzystają z
`clang++` (niezależnie od produkcyjnego GCC); można to zmienić przez
`SANITIZER_CXX=...`.

Metadane bibliograficzne z niepustym DOI można porównać automatycznie z
Crossref. Kontrola sprawdza składnię DOI, kanoniczny adres `doi.org`, istnienie
rekordu oraz zgodność tytułu po normalizacji Unicode i interpunkcji. Tytuł
z Crossref jest wcześniej oczyszczany z formatowania wydawcy: znaczników HTML,
numeru artykułu starszych czasopism (`LXXIX.` dla Rutherforda 1911) oraz
zapisu greckimi literami tam, gdzie katalog używa nazwy litery (`α` = `alpha`):

```bash
make references-check
```

Test wymaga dostępu do sieci; rekordy źródeł internetowych bez DOI są jawnie
raportowane jako pominięte, a nie uznawane za zweryfikowane. Dwa rodzaje
niepowodzenia są rozróżniane kodem wyjścia: `1` oznacza rzeczywistą niezgodność
metadanych z Crossref (błędny tytuł, adres niekanoniczny, DOI niezarejestrowany),
a `2` — że kontroli nie dało się przeprowadzić, bo katalog jest nieczytelny albo
Crossref nieosiągalny. Kod `2` nigdy nie twierdzi, że bibliografia jest błędna;
mówi tylko, że pozostała niezweryfikowana.

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
\(P_\mu=\mu_0|\ddot{\boldsymbol\mu}|^2/(6\pi c^3)\). W trybie
`--radiation-reaction stochastic` (domyślnym) ta moc **nie jest** odprowadzana
w sposób ciągły — patrz „Kwantyzacja wszystkich kanałów" niżej. Odpowiadający mu
moment reakcyjny jest wyznaczany z \(\boldsymbol\mu\times
\dddot{\boldsymbol\mu}\), działa na orientację dipola, a energia jest
jednocześnie przenoszona z sektora wewnętrznego do energii promieniowania.
Bilans obejmuje także konwekcyjny pęd promieniowania dipolowego i wynoszony
moment pędu.

\[
\frac{d\boldsymbol\mu_i^{\rm proper}}{dt}
=\frac{q_i}{m_i}\boldsymbol\mu_i^{\rm proper}\times\mathbf B_{\mathrm{BMT},i},
\]

gdzie \(\boldsymbol\mu_i^{\rm proper}\) jest momentem w **chwilowym układzie
spoczynkowym** cząstki (`state.firstProperDipole`), bo to jego dotyczy
równanie Jacksona 11.170 — moment laboratoryjny `state.firstDipole` jest
z niego wyprowadzany boostem tensorowym i **nie** spełnia tego równania.
Pomylenie tych dwóch wektorów było błędem naprawionym w drugim audycie
(patrz „audyt fizyki" niżej). Dla \(a=(g-2)/2\),

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
Ponieważ w każdym półkroku \(\mathbf B_{\rm BMT}\) jest zamrożone, równanie
ma postać \(\dot{\boldsymbol\mu}=\boldsymbol\omega\times\boldsymbol\mu\)
ze stałym \(\boldsymbol\omega=-(q/m)\mathbf B_{\rm BMT}\), więc obrót jest
wykonywany **wzorem Rodriguesa w dwóch symetrycznych półkrokach** — jest to
rozwiązanie ścisłe, nie całkowanie numeryczne, dla dowolnej długości kroku.
Zachowuje \(|\boldsymbol\mu_i|=(g/2)\mu_B\) z samej konstrukcji, bez
renormalizacji.

### Bilans energii: co naprawdę mierzy 5,9%

Macierz bilansu długiego horyzontu raportowała dla ścieżki produkcyjnej
(retardowana + sama siła własna) residuum **5,9%**, i przez długi czas
figurowało to jako jakość domknięcia bilansu energii przez model. **Nie jest
to własność fizyczna modelu, tylko błąd dyskretyzacji wzmocniony przez
kondycjonowanie samego pomiaru.**

Residuum dzieli różnicę dwóch energii mechanicznych przez energię pola
dalekiego, mniejszą o cztery rzędy. Na promieniu sondy energia wiązania pary
wynosi 34,0 eV, a wypromieniowana przez dwa okresy 2,475·10⁻³ eV, więc
**każdy** względny błąd energii mechanicznej jest przed zaraportowaniem
mnożony przez

\[
\frac{|E_{\rm mech}|}{E_{\rm far}}\approx 1{,}37\cdot10^{4}.
\]

Przy tym wzmocnieniu błąd 4,3·10⁻⁶ — czyli dokładnie skala tolerancji
integratora, \(10^{-6}\) — daje odczyt 5,9%. Nie trzeba więc szukać
brakującej fizyki: sam dryf integratora w pełni tłumaczy tę liczbę.

Trzy niezależne pomiary mówią to samo:

- **Sektor reakcji jest poprawny do trzech cyfr.** Siła Landaua–Lifshitza
  oddaje 0,4997–0,4998 strumienia dalekiego wobec teoretycznych 0,5000 dla
  tej pary (przy \(\mathbf a_1=-\mathbf a_2\) moc koherentna
  \(|q_1\mathbf a_1+q_2\mathbf a_2|^2\) jest dokładnie dwukrotnością sumy
  członów własnych).
- **Zbieżne jest to, co fizyczne, a błądzi to, co numeryczne.** \(E_{\rm far}\)
  jest zbieżne do 2·10⁻⁵ przez cztery dekady tolerancji, podczas gdy różnica
  mechaniczna waha się w tym samym zakresie o 12%.
- **Zagęszczanie dyskretyzacji przepycha residuum przez zero.** Rząd 2 przy
  \(10^{-8}\) daje −1,0%, rząd 4 przy \(10^{-10}\) daje +2,9%, wobec −5,9%
  przy ustawieniu produkcyjnym. Prawdziwy wyciek nie zmienia znaku przy
  poprawianiu numeryki.

Próba naprawy przez podniesienie rzędu interpolacji historii retardowanej —
kubiczny Hermite (pozycja, prędkość; przyspieszenie odtwarzane przez dwukrotne
różniczkowanie, więc \(O(h^2)\)) na **kwintyczny** (dopasowujący także
zapisane przyspieszenia, \(O(h^4)\)) — została wykonana i **cofnięta**. Rząd
zbieżności faktycznie wzrósł z 1,999 do 3,999 na analitycznej linii świata,
ale residuum bilansu **pogorszyło się** do −17,4%, a zbieżność trajektorii
przy tolerancji \(10^{-7}\) się załamała. Powód: zapisane przyspieszenia są
wzajemnie spójne z próbkami pozycji i prędkości tylko na tyle, na ile czyni je
integrator rzędu drugiego, więc interpolant klasy \(C^2\) je przetrenowuje.
Rząd interpolacji nie jest tu wąskim gardłem.

Zamiast tego bramka `long-horizon-radiative-balance` została przeformułowana
na stwierdzenie o **zbieżności**, a nie o zachowaniu energii: obok odczytu
zgrubnego liczony jest ten sam bilans o krok gęstszy (tolerancja \(10^{-8}\)),
i wymaga się, żeby zagęszczenie nie pogorszyło wyniku i zeszło poniżej 5%.
Zmierzone: **5,94% → 0,98%** po jednym kroku zagęszczenia, koszt 1,5 s.

Test ma moc rozróżniającą — sprawdzone przez wstrzyknięcie **realnego** błędu
fizycznego (siła reakcji zawyżona o 10%): residuum zgrubne 10,9%, po
zagęszczeniu 5,98%, czyli **nie znika**, tylko zostaje na skali wstrzykniętego
błędu, i test oblewa. Błąd dyskretyzacji spada sześciokrotnie; realny wyciek
nie spada.

Raportowane jest też samo wzmocnienie (`balance amplification`), żeby liczby
z tej macierzy nie dało się już odczytać bez kontekstu.

#### Przemiatanie promienia orbity

Żeby żaden z tych wniosków nie był przypadkiem jednego promienia, ten sam
bilans jest liczony w trzech punktach rozpiętych na czynnik 16, każdy w obu
dyskretyzacjach:

| \(R/a_{\rm pary}\) | wzmocnienie | zgrubnie | zagęszczone |
|---|---|---|---|
| 0,05 | 1 718 | 7,24% | 2,61% |
| 0,2 | 13 741 | 5,94% | 0,98% |
| 0,8 | 109 912 | 23,6% | 6,77% |

Dwie rzeczy naraz.

**Residuum maleje po zagęszczeniu na każdym promieniu**, nie tylko na tym, z
którego cytowana jest liczba nagłówkowa. Zauważmy też, że samo residuum
zgrubne **nie skaluje się z \(\beta\)** — idzie 7,2%, 5,9%, 23,6% — tylko
podąża za wzmocnieniem, dokładnie jak dla artefaktu kondycjonowania.

**Wzmocnienie podąża za prawem \(R^{3/2}\)** — i to jest test fizyczny, nie
księgowy. Ponieważ \(E_{\rm mech}\propto1/R\), a z Larmora
\(E_{\rm far}\) na orbitę \(\propto R^{-5/2}\), iloraz musi rosnąć jak
\(R^{3/2}\), czyli 8,0 na czterokrotność promienia. Zmierzone: **7,997 i
7,998**. Sprawdzenie to obowiązuje z pasmem 2%, więc wyłapuje złamanie prawa
potęgowego inspiralu, którego żadna sonda jednopunktowa tu nie widzi.

Moc rozróżniająca sprawdzona przez wstrzyknięcie fałszywej zależności
promieniowanej energii od promienia: ilorazy spadają do 7,766 i 7,171
(odchyłki 2,9% i 10,4%) i test oblewa. Uczciwe zastrzeżenie — ta konkretna
wada trafia również w `larmor-normalization`; unikalna wartość przemiatania
dotyczy wad **poprawnie skalibrowanych na jednym promieniu, a źle
skalujących się**.

Koszt: +2,8 s.

### Kwantyzacja wszystkich kanałów promieniowania

W trybie `stochastic` (domyślnym, `gRadiationReactionModel`) **żaden** kanał
promieniowania nie odprowadza energii w sposób ciągły. Cała moc klasyczna —
sektor ładunkowy E1 **oraz** magnetyczny dipol M1 — jest sumowana i bankowana
jako jeden strumień hazardu Poissona, wypłacany dyskretnymi kwantami
\(\hbar\omega\):

\[
\frac{d\Lambda}{dt}=\frac{P_{E1}+P_{M1}+P_{E2}}{\hbar\omega_{\rm orb}} .
\]

Moce multipolowe są w tym rzędzie ortogonalne, więc jest to suma, a nie
podwójne liczenie.

Wcześniej skwantowany był wyłącznie E1, a E2 **nie promieniowało w ogóle**,
w żadnym trybie: jego moc zasilała jedynie bramkę dominacji modelu
`automatic` i nigdzie nie wchodziła do bilansu. M1 miał **dwa** ujścia ciągłe, w
dodatku bez żadnej bramki na model reakcji: moment siły reakcyjnej
(`applyDipoleRadiationTorque`) oraz drenaż wewnętrznego rezerwuaru
(`dipoleConstraintEnergy`). Oba są teraz wyłączane w trybie skwantowanym —
inaczej energia M1 wychodziłaby dwa razy, raz ciągle i raz jako fotony.

Wspólna częstość \(\omega_{\rm orb}\) to świadomy koszt jednego kanału: M1
jest w rzeczywistości emitowany blisko tempa precesji spinu, o kilka rzędów
**niższego** niż \(\omega_{\rm orb}\), więc jego kwanty wychodzą rzadsze i
większe, niż dałoby tamto tempo. Nie zmienia to jednak **średniej** energii
wypromieniowanej: tempo hazardu wynosi \(P/\hbar\omega\), a każdy foton unosi
\(\hbar\omega\), więc oczekiwana energia na jednostkę czasu równa się \(P\)
dla **dowolnego** \(\omega\). Od wyboru zależy tylko ziarnistość i statystyka
zliczeń.

Konsekwencja, którą trzeba wypowiedzieć wprost: energia M1 opuszcza teraz układ
przez orbitalny odrzut fotonu, a nie przez tłumienie sektora spinowego. Przy
zmierzonych udziałach M1 (poniżej) jest to praktycznie nieobserwowalne, ale jest
to realna zmiana drogi przepływu energii.

`radiatedEnergy` pozostaje **ciągłą** całką Poyntinga i pełni rolę diagnostyki
porównawczej, nie sumy wyemitowanych kwantów — tak jak dla E1 od początku.

#### Anihilacja związana z promieniem terminalnym

Generator anihilacji w tym projekcie jest receptą kwantową **celowo
niezależną** od modelu klasycznego: zakłada parę **w spoczynku** i wkłada w
fotony pełne \(2m_ec^2\), przez co linia \(2\gamma\) jest stałą kompilacji,
a widmo \(3\gamma\) to Ore-Powell dla pary spoczywającej. Ani sufit
kinematyczny, ani podłoga go nie wyzwalają.

Obok tego liczona jest teraz **druga, równoległa odpowiedź**: co zostaje
parze, którą ten model faktycznie scałkował, w chwili gdy zatrzymuje ją
reguła stopu. Orbita jest wtedy związana, więc jej energia niezmiennicza leży
**poniżej** sumy mas spoczynkowych o nagromadzone wiązanie:

\[
W=(m_1+m_2)c^2+\mu\varepsilon,\qquad \varepsilon<0,
\]

i to \(W\), a nie \(2m_ec^2\), mają do podziału fotony stanu końcowego.
Para dzieli je na \(2\gamma\) po \(W/2\) (wymusza to samo zachowanie pędu),
orto na \(3\gamma\) z widma Ore-Powella, którego jedyną skalą jest \(W/2\) —
więc przeskalowuje się z \(W\) zamiast być przypięte do \(m_ec^2\).

Zmierzone (seed 4242, para):

| wielkość | wartość |
|---|---|
| promień terminalny | \(258\) fm |
| energia wiązania | \(2{,}79\) keV |
| \(W\) | \(1019{,}21\) keV wobec \(1022{,}00\) keV → **\(-0{,}27\%\)** |
| \(W/2\) (linia \(2\gamma\), koniec widma \(3\gamma\)) | \(509{,}60\) keV wobec \(510{,}999\) keV |

Przesunięcie linii \(511\) keV jest wielkością **mierzoną doświadczalnie**,
więc jest to konkretna, sprawdzalna przepowiednia — wynikająca jednak z
klasycznej energii wiązania, nie z modelu anihilacji, którego CREM nie ma.

*Dwie liczby, które łatwo pomylić.* Wielkością, którą przesuwa wiązanie, jest
\(W\) (a więc i \(W/2\)) — to samo dla obu kanałów. **Wiodący foton to co
innego**: dla orto leży średnio \(12\%\) poniżej \(W/2\) z samej przestrzeni
fazowej trzech fotonów (zmierzone: \(0{,}878\pm\) rozrzut \(0{,}73\)–\(1{,}00\)),
co z wiązaniem nie ma nic wspólnego. Pierwsza wersja tego raportu porównywała
wiodący foton z \(511\) keV dla obu kanałów, przez co orto wyglądało na
przesunięte \(-8{,}7\%\), czyli trzydziestokrotnie bardziej niż para —
mieszając dwa niezależne efekty. Raportowane są teraz osobno.

*Kontrola.* Suma energii fotonów równa się \(W\) do dwunastu cyfr w obu
kanałach, a każdy foton spełnia \(E\le W/2\), czyli warunek fizycznej
przestrzeni fazowej. Wypełniane tylko dla trajektorii, które **realnie
dochodzą** do promienia terminalnego; przebiegi ocenzurowane i nieudane nie
mają stanu, z którego można anihilować.

*Człon dipol-dipol: duży, ale o zerowej wartości oczekiwanej.* Element
`elements.specificEnergy` jest **oskulacyjnym elementem keplerowskim**, nie
ogólną energią: osiem miejsc czyta go jako \(a=-K/2\varepsilon\),
\(e^2=1+2\varepsilon\ell^2/K^2\) i okres. Wsunięcie w niego członu
\(1/r^3\) nakarmiłoby każdy z tych wzorów wielkością, dla której nie zostały
wyprowadzone, a zepsucie byłoby największe **dokładnie przy promieniu
terminalnym**, gdzie sektor dipolowy stanowi \(56\%\) kulombowskiego.
Element zostaje więc keplerowski, a energia dipolowa jest dodawana do
niezmiennika, który dzielą fotony:

\[
W=(m_1+m_2)c^2+\mu\varepsilon+U_{\rm dd}.
\]

Nie jest to korekta **dynamiki**: trajektoria już czuje tę siłę
(`mutualForces` niesie `regularizedDipoleForce`), a mechaniczna księga już
niesie tę energię (`conservativeParticleEnergy` dodaje `dipolePotential`).
Pominięta była wyłącznie **etykieta energetyczna** księgowania sekularnego.

Stosunek dipol/Coulomb rośnie jak \(1/r^2\): \(3{,}3\cdot10^{-6}\) przy
\(a_{Ps}\), ale \(0{,}56\) przy promieniu terminalnym i \(0{,}998\) na
barierze. Warto odnotować, że \(r^*=193{,}3035\) fm, gdzie oba człony się
zrównują, jest **co do cyfry** równe barierze Comptona.

*Azymut trzeba uśrednić, nie zgadnąć.* `osculatingPeriapsisState` resetuje
płaszczyznę orbity do kanonicznej x-y i kładzie separację wzdłuż \(+\hat x\),
mówiąc wprost, że orientacja płaszczyzny „jest nieistotna dla estymaty czasu
kolapsu". Dla czasu kolapsu owszem — ale **nie dla czynnika kątowego**
\(\boldsymbol\mu_1\!\cdot\!\boldsymbol\mu_2-3(\boldsymbol\mu_1\!\cdot\!\hat{\mathbf n})(\boldsymbol\mu_2\!\cdot\!\hat{\mathbf n})\),
który zależy od kierunku \(\hat{\mathbf n}\). Pierwsza wersja tej sekcji
korzystała z tamtej konwencji i podawała liczby, które sprawdzone przeciw
poprawnej średniej miały **inny znak** (na losowych konfiguracjach
\(-0{,}268\) wobec \(+0{,}484\)).

Prawdziwa anomalia przy terminacji jest w reprezentacji oskulacyjnej
nieznana — kod mówi to sam przy losowaniu azymutu emisji fotonu. Śledzona
jest natomiast normalna płaszczyzny, więc azymut jest **uśredniany**:

\[
\langle n_in_j\rangle=\tfrac12(\delta_{ij}-\hat L_i\hat L_j)
\ \Longrightarrow\
\langle\text{czynnik}\rangle=-\tfrac12(\boldsymbol\mu_1\!\cdot\!\boldsymbol\mu_2)
+\tfrac32(\boldsymbol\mu_1\!\cdot\!\hat{\mathbf L})(\boldsymbol\mu_2\!\cdot\!\hat{\mathbf L}),
\]

zweryfikowane przeciw całkowaniu brute-force po azymucie do \(10^{-16}\).
Nie wchodzi ważenie \(1/r^3\): \(r\) jest ustalone na perycentrum
terminalnym, więc jest to średnia po azymucie przy jednym promieniu, nie
średnia po orbicie.

*Co dokładnie jest uśredniane — bo to łatwo pomylić.* Uśredniana jest
orientacja **sztywnej pary momentów względem \(\hat{\mathbf L}\)**, przy
**ustalonym** wzajemnym kącie \(\cos=\hat{\boldsymbol\mu}_1\!\cdot\!\hat{\boldsymbol\mu}_2\).
Ten kąt jest właśnie tym, co odróżnia kanały, i nie jest uśredniany. Tożsamość
niżej sprawdzono osobno dla \(\cos=-0{,}9;\,-0{,}4;\,0;\,+0{,}4;\,+0{,}9\),
więc zerowanie zachodzi **wewnątrz każdej klasy**, a nie przez wymieszanie ich
ze sobą.

*Uwaga o tym, który kanał ma zgodne momenty.* Wcześniejsza wersja tego opisu
mówiła „para: antyrównoległe, orto: równoległe" — to dotyczy **spinów**. Dla
**momentów magnetycznych** jest odwrotnie, bo przeciwne ładunki odwracają
korelację spin–moment (\(\boldsymbol\mu=(gq/2m)\mathbf S\)). Kod wymusza
\(\cos\ge0{,}5\) dla para i \(\cos<0{,}5\) dla orto, co README dokumentuje
osobno jako Sondę 4. Zmierzone przy terminacji potwierdza konwencję: para
\(+0{,}53\ldots+0{,}89\), orto \(-0{,}77\ldots+0{,}47\). Warunek orto jest
przy tym **jednostronny** — to nie „antyrównoległe", lecz „wszystko poza
silnie zgodnym", od prostopadłych po przeciwne.

*Czy założenie izotropii przeżywa dynamikę.* Izotropia orientacji względem
\(\hat{\mathbf L}\) wynika z warunków początkowych, ale trzeba sprawdzić, czy
ewolucja jej nie łamie — i czy sam \(\cos\) dożywa promienia terminalnego.
Nasuwa się podejrzenie, że nie: przeciwne ładunki dają przeciwne \(q/m\), więc
składowe prostopadłe mogłyby kontrarotować i rozmyć wzajemny kąt. **Zmierzone,
\(\cos\) jest zachowany**:

```
para:  0,885 -> 0,893    0,552 -> 0,552    0,539 -> 0,535    0,798 -> 0,648
orto: -0,733 -> -0,738  -0,391 -> -0,391  -0,765 -> -0,765  -0,117 -> +0,028
```

Największa zmiana to \(\approx0{,}15\), większość poniżej \(0{,}01\). Powód
jest ten sam, co odwrócenie z Sondy 4: skoro przeciwne ładunki odwracają
relację spin–moment, to \(q/m\) i orientacja momentu odwracają się **razem**,
więc oba momenty precesują w tę samą stronę wokół wspólnego
\(\mathbf B_{\rm BMT}\).

*Dlaczego momenty się nie ustawiają — poprawiony argument.* Wcześniejsza
wersja tego akapitu tłumaczyła to tak: \(\mathbf B_{\rm BMT}\) jest
zdominowane przez człon motoryczny \(\mathbf v\times\mathbf E/c^2\), skierowany
wzdłuż \(\hat{\mathbf L}\), więc precesja zachowuje
\(\boldsymbol\mu\!\cdot\!\hat{\mathbf L}\). **To jest prawdą tylko przy dużych
separacjach.** Zmierzony stosunek pola dipolowego partnera do motorycznego
wynosi \(0{,}25\) przy \(a_{Ps}\), zrównuje się przy \(6630\) fm i sięga
\(5{,}07\) przy promieniu terminalnym — czyli tam, gdzie ten człon w ogóle
zaczyna być duży, precesją rządzi już pole partnera, a nie \(\hat{\mathbf L}\).
Argument nie obowiązuje dokładnie tam, gdzie był potrzebny.

Nie obowiązuje też odpowiedź „nie zdąży": okres precesji dipolowej przy
\(a_{Ps}\) wynosi \(45{,}6\) ps wobec \(\approx133\) ps kolapsu, czyli około
trzech pełnych obrotów. **Czasu jest dość.**

Nie broni się też odpowiedź „nie ma czym stracić energii", którą wpisano tu
w pierwszym podejściu do tej poprawki. Kanał **istnieje**: M1 wchodzi do
kwantu na równi z E1 i E2,

```cpp
const double quantizedPower=
    stepRadiation.leadingElectricDipolePower
    +stepRadiation.magneticDipoleFlux.energy
    +stepRadiation.electricQuadrupolePower;
```

a scałkowany po całym kolapsie kanał M1 unosi \(\approx7\) eV wobec
\(3837\) eV z E1 — udział \(1{,}89\cdot10^{-3}\) uśredniony po orientacjach
momentów, mało, ale nie zero.
Momenty **mogą** tracić energię w kwantach.

*Ile faktycznie się obracają — liczba rozstrzygająca.* Wyprowadźmy to z
momentu siły reakcji M1, tego samego, który integrator stosuje jako
`mu += gamma*T*dt`. Zmierzona chwilowa szybkość obrotu wzdłuż spirali:

| promień | \(\Omega_{\rm reakcji}\) [rad/s] | \(\Omega_{\rm reakcji}/\Omega_{\rm BMT}\) | \(1/\Omega_{\rm reakcji}\) [s] |
|---|---|---|---|
| \(a_{Ps}\) | \(0{,}85\) | \(2{,}1\cdot10^{-12}\) | \(1{,}18\) |
| \(a_{Ps}/10\) | \(5{,}2\cdot10^{5}\) | \(2{,}2\cdot10^{-9}\) | \(1{,}93\cdot10^{-6}\) |
| \(6630\) fm | \(2{,}6\cdot10^{8}\) | \(2{,}8\cdot10^{-7}\) | \(3{,}85\cdot10^{-9}\) |
| \(258\) fm (terminalny) | \(1{,}0\cdot10^{16}\) | \(1{,}2\cdot10^{-3}\) | \(9{,}90\cdot10^{-17}\) |
| \(193\) fm (bariera) | \(2{,}9\cdot10^{16}\) | \(2{,}1\cdot10^{-3}\) | \(3{,}48\cdot10^{-17}\) |

Przy promieniu terminalnym czas ustawienia to \(10^{-16}\) s — o czternaście
rzędów krócej niż kolaps. Gdyby patrzeć tylko na tę kolumnę, momenty
ustawiałyby się natychmiast, i zarzut byłby w pełni słuszny.

Rozstrzyga jednak **scałkowany** obrót, bo tam, gdzie \(\Omega_{\rm reakcji}\)
eksploduje, para prawie nie przebywa. Całkując \(\int\Omega_{\rm
reakcji}\,dt\) po spirali \(dt=(dE/da)\,da/P_{E1}\) od \(a_{Ps}\) do
bariery:

```
całkowity czas spirali = 3,16e-11 s
```

Czas życia spirali jest zdominowany przez obszar zewnętrzny —
\(3{,}16\cdot10^{-11}\) s upływa już przy \(a\approx2\cdot10^{-11}\) m i dalej
się praktycznie nie zmienia, podczas gdy szybki obrót zaczyna się cztery rzędy
niżej.

> **Sprostowanie, dwukrotne — i drugie znosi pierwsze.** W tym miejscu stała
> najpierw liczba \(5{,}24\cdot10^{-3}\) rad \(=0{,}30°\), potem —
> po uśrednieniu \(\int\Omega_{\rm reakcji}\,dt\) po \(96\) losowych
> parach orientacji — \(2{,}33\cdot10^{-2}\) rad \(=1{,}33°\). **Obie są
> nieaktualne, bo zły był sam estymator, nie jego uśrednienie.**
>
> \(\int\Omega_{\rm reakcji}\,dt\) sumuje **moduł** chwilowej szybkości
> obrotu, czyli milcząco zakłada dryf sekularny. Moment siły reakcji ma
> jednak kierunek \(\propto\boldsymbol\mu\times\dddot{\mathbf m}_{\rm tot}\),
> a \(\dddot{\mathbf m}_{\rm tot}\) obraca się razem z orbitą — obrót w
> kolejnych częściach okresu znosi się, zostaje **oscylacja** o amplitudzie
> \(\sim\Omega_{\rm reakcji}/\omega_{\rm orb}\), a nie narastanie
> \(\Omega_{\rm reakcji}\cdot t\).

*Pomiar bezpośredni, zamiast całki.* Przebieg `stochastic` i przebieg
`disabled` z tego samego stanu początkowego różnią się **wyłącznie** momentem
siły M1 — dopóki nie padnie foton, trajektorie są identyczne co do bitu — więc
\(|\Delta\hat{\boldsymbol\mu}|\) izoluje go czysto. Przy \(a_{Ps}/200\):

| orbity | \(1\) | \(2\) | \(4\) | \(8\) | \(16\) |
|---|---|---|---|---|---|
| \(\|\Delta\hat{\boldsymbol\mu}\|\) | \(2{,}16\cdot10^{-6}\) | \(2{,}14\cdot10^{-6}\) | \(2{,}17\cdot10^{-6}\) | \(1{,}84\cdot10^{-6}\) | \(2{,}72\cdot10^{-6}\) |

**Ograniczone, bez trendu liniowego** — szesnastokrotne wydłużenie czasu nie
zwiększa efektu. To samo przy \(a_{Ps}/50\): \(3{,}02\cdot10^{-9}\) po
dwóch orbitach i \(2{,}77\cdot10^{-9}\) po ośmiu. Rzędy zgadzają się z
\(\Omega_{\rm reakcji}/\omega_{\rm orb}\) (\(\approx8\cdot10^{-10}\) przy
\(a_{Ps}/50\)), a nie z \(\Omega_{\rm reakcji}\cdot t\), która byłaby o
rzędy większa.

Skala pominięcia wynosi więc \(\Omega_{\rm reakcji}/\omega_{\rm orb}\):
\(4\cdot10^{-17}\) rad przy \(a_{Ps}\), \(2\cdot10^{-6}\) przy
\(a_{Ps}/200\), a w najgorszym punkcie — przy promieniu terminalnym —
\(\approx6\cdot10^{-5}\) rad, czyli \(0{,}003°\). Nie \(1{,}33°\).

> Uwaga metodologiczna: liczby zbieżności siatki
> (\(4{,}78/4{,}44/4{,}36/4{,}10/4{,}12\cdot10^{-3}\) dla
> \(N=36/100/300/800/2000\)) były poprawne, i to właśnie czyniło je mylącymi
> — zbieżna całka źle postawionej wielkości wygląda dokładnie tak samo jak
> zbieżna całka dobrze postawionej.

### Luka implementacyjna wokół kanału M1 — rozbiór

W trybie skwantowanym bramkowane są **dwie** rzeczy, przy tej samej fladze i w
tym samym miejscu: moment siły reakcji
(`if(!quantizedRadiation) applyDipoleRadiationTorque(...)`) oraz drenaż
rezerwuaru (`if(!quantizedRadiation) trial.dipoleConstraintEnergy-=...`).
Wypłata fotonu, `applyStochasticDipolePhoton`, przelicza wyłącznie czteropędy —
nie dotyka `firstProperDipole` ani `secondProperDipole`. Pierwotny opis tej
luki („energia wychodzi przez kinematykę, a orientacja nie jest obciążana")
był zbyt zgrubny i po rozbiorze na składniki okazał się mylący. Poniżej cztery
aspekty, zmierzone osobno.

**(a) Bilans energii — usterka realna, ale w innym miejscu, niż ją najpierw
umiejscowiłem; poprawiona w kodzie.**

Usterka: obie ciągłe studnie sektora M1 — moment siły reakcji i drenaż
`dipoleConstraintEnergy` — były bramkowane **wyłącznie** dla
`stochasticElectricDipole`. Tymczasem sektor ładunkowy od zawsze wyklucza
**dwa** modele, `disabled` i `stochasticElectricDipole`, z komentarzem
mówiącym wprost, że w `disabled` „nic nie ma ciągnąć orbity". A `disabled`
to właśnie model, którym `crem_collapse.hpp` całkuje przebieg **tła**,
odejmowany potem od przebiegu rzeczywistego.

Skutki, zmierzone na jednym okresie przy \(a_{Ps}/50\) (przed poprawką):

```
                       tryb ciągły      tło (disabled)    tryb skwantowany
dipoleConstraintEnergy  -4,81317e-27      -4,81336e-27       0,00000e+00
```

Tło drenowało ten sam rezerwuar, co przebieg, od którego je odejmowano, więc
w modelach ciągłych studnia M1 kasowała samą siebie — przez odjęcie
przechodziło \(4{,}0\cdot10^{-5}\) jej wartości. W trybie produkcyjnym
(skwantowanym) było gorzej: tło niosło studnię, której przebieg rzeczywisty
w ogóle nie miał, więc między fotonami oba przebiegi **nie były identyczne**,
choć powinny być.

Poprawka: bramkowanie **drenażu** objęło `disabled`, dokładnie jak w sektorze
ładunkowym. Po niej tło ma \(0\), a przez odejmowanie przechodzi \(100\%\)
energii M1.

(Bramkowanie momentu siły **rozdzielono** od bramkowania drenażu — to
poprawka (b) niżej. Przebieg skwantowany i jego tło nie są więc między
fotonami identyczne co do bitu: przebieg rzeczywisty niesie moment reakcji
M1, którego tło z definicji nie ma. Tak ma być — to sygnał fizyczny, a nie
artefakt dyskretyzacji, i odejmowanie tła ma go zachować, nie skasować.)

*Ile to zmienia w wynikach — niewiele, i warto powiedzieć dlaczego.*
Prześledzenie ścieżki energii pokazuje, że rezerwuar i tak nie zasilał
kolapsu:

- pole `specificEnergy` z `measuredDelta` jest liczone, ale **celowo nieczytane**
  (komentarz w `crem_collapse.hpp`: „ANGULAR MOMENTUM only, from here on");
  tło jest odejmowane wyłącznie od momentu pędu;
- sekularny ubytek energii to
  `deltaEnergyPerOrbit = -finalState.orbitalRadiatedEnergy/reducedMass`, a
  `orbitalRadiatedEnergy` jest z definicji strumieniem **bez M1**
  (`radiatedEnergyIncrement-dipoleRadiatedEnergy`);
- w ścieżce produkcyjnej nawet to nie jest używane: `lossPerOrbit` bierze
  `expectedLossPerOrbit`, czyli analityczną moc Larmora **czystego E1**.

Zostaje więc jedno wejście M1 do dynamiki: hazard wewnątrz **jednej**
mierzonej orbity na checkpoint, gdzie `quantizedPower` sumuje E1, M1 i E2.
Skip analityczny pokrywa do \(200000\) orbit na checkpoint i jest w całości
E1-owy, a mierzona orbita niesie \(\sim5\cdot10^{-6}\) hazardu. Przy udziale
M1 rzędu \(1{,}9\cdot10^{-3}\) daje to wpływ M1 na kolaps rzędu
\(10^{-8}\) — poprawka jest więc naprawą spójności, nie liczb.

> Przy okazji poprawiony nieaktualny komentarz w `state.hpp`, który nazywał
> `orbitalRadiatedEnergy` „diagnostic-only: nothing currently reads it in
> production". Czyta je `crem_collapse.hpp` i jest to **jedyny** produkcyjny
> kanał energii sekularnej. Zdanie „M1 does not recoil the orbit at all" jest
> prawdziwe o księdze sekularnej i fałszywe o hazardzie skwantowanym.

**(b) Reakcja na orientację — jedyny realnie pominięty składnik; domknięty
w kodzie, ale mały.** Moment siły reakcji był odrzucany i był to jedyny
składnik, którego skwantowanie nie odtwarzało żadną drogą. Jego skala to
jednak \(\Omega_{\rm reakcji}/\omega_{\rm orb}\), czyli najwyżej
\(\approx6\cdot10^{-5}\) rad przy promieniu terminalnym i \(10^{-17}\)
rad przy \(a_{Ps}\) — obrót **oscylacyjny**, nie narastający (pomiar i
sprostowanie wyżej).

Przyczyną było sklejenie dwóch **różnych** decyzji pod jedną flagą. Drenaż
energii trzeba wyłączyć w trybie skwantowanym, bo tę samą energię unosi
foton — inaczej odjęto by ją dwa razy. Ale reakcja **orientacyjna** nie ma
takiego duplikatu: `applyStochasticDipolePhoton` przelicza czteropędy i nic
więcej, foton niesie wyłącznie odrzut liniowy i zostawia moment, który go
wypromieniował, skierowany dokładnie tam, gdzie był. Bramkowanie momentu siły
razem z drenażem nie zapobiegało więc podwójnemu liczeniu — po prostu
kasowało jedyną reakcję, do której ścieżka skwantowana nie miała innej drogi.

Poprawka rozdziela flagi:

```cpp
const bool quantizedRadiation=            // drenaż energii
    model==stochasticElectricDipole || model==disabled;
const bool applyDipoleReactionTorque=     // reakcja orientacyjna
    model!=disabled;
```

Energia nie jest przez to liczona podwójnie: kosztem energetycznym samego
obrotu jest zmiana energii oddziaływania dipol-dipol, którą
`conservativeParticleEnergy` już niesie przez
`regularizedDipoleInteractionEnergy`. Foton rozlicza energię wypromieniowaną
do strefy falowej, moment siły — orientację. `disabled` nie dostaje ani
jednego, ani drugiego, zgodnie z (a).

**(c) Czy kanał M1 stać na to, co promieniuje.** Skumulowana energia M1 wzdłuż
spirali nigdy nie przekracza energii orientacyjnej dipoli dostępnej przy tym
samym promieniu — stosunek rośnie monotonicznie, ale kończy na
\(3{,}9\cdot10^{-3}\) (\(3{,}70\) eV skumulowanego M1 wobec \(950\) eV
\(U_{dd}\) przy barierze). Żadne z dwóch kont nie jest przekroczone, więc
pytanie „skąd ta energia pochodzi" nie ma w tym modelu wymuszonej odpowiedzi.

**(d) Moment pędu — pomijalny.** Udział M1 w strumieniu wypromieniowanego
momentu pędu: \(2{,}2\cdot10^{-15}\) przy \(a_{Ps}\),
\(3{,}5\cdot10^{-11}\) przy \(a_{Ps}/10\), \(9{,}3\cdot10^{-5}\) przy
promieniu terminalnym. Po zważeniu czasem przebywania — nieistotny.

*Czy luka rozdziela kanały.* Destruktywna interferencja M1 znosi się
**dokładnie** tylko przy \(\cos=-1\); realne zdarzenia orto zajmują całe
\(\cos<0{,}5\). Uśredniony po orientacjach udział M1 w strumieniu
skwantowanym wynosi \(1{,}875\cdot10^{-3}\pm3{,}7\cdot10^{-4}\) dla para i
\(1{,}903\cdot10^{-3}\pm2{,}3\cdot10^{-4}\) dla orto — różnica
\(2{,}8\cdot10^{-5}\) przy niepewności \(4{,}4\cdot10^{-4}\), czyli
**zgodna z zerem**. Obraz „M1 tylko dla para, zero dla orto" jest prawdziwy
wyłącznie w punkcie \(\cos=-1\) i nie przenosi się na statystykę. Luka nie
wprowadza mierzalnej asymetrii kanałów.

Zmierzone wprost, rozkład \(\boldsymbol\mu\!\cdot\!\hat{\mathbf L}\) na
\(144\) składowych z \(72\) trajektorii:

| | średnia | SD |
|---|---|---|
| na starcie | \(-0{,}081\) | \(0{,}601\) |
| przy terminacji | \(-0{,}090\) | \(0{,}596\) |
| izotropowo | \(0\) | \(1/\sqrt3=0{,}577\) |

Mediana zmiany \(\boldsymbol\mu\!\cdot\!\hat{\mathbf L}\) przez cały kolaps to
\(0{,}019\); frakcja blisko \(\hat{\mathbf L}\) \(12\%\) wobec \(10\%\)
izotropowych. Momenty precesują, ale kończą tam, gdzie zaczęły, z rozkładem
nieodróżnialnym od izotropowego. Reakcja M1 wnosi do tego co najwyżej
\(\approx6\cdot10^{-5}\) rad (oscylacja, nie dryf — patrz pomiar wyżej),
czyli **trzysta razy mniej** niż zmierzone \(0{,}019\). Werdykt o izotropii
ma więc zapas trzech rzędów wielkości.

> Poprzednie dwie wersje tego akapitu podawały tu \(5{,}2\cdot10^{-3}\), a
> potem \(2{,}33\cdot10^{-2}\) rad, i na tej drugiej podstawie **wycofałem**
> stwierdzenie o wygodnym zapasie, pisząc, że domknięcie luki podwoiłoby
> obrót momentów. To wycofanie samo było błędne i zostaje wycofane: zapas
> jest wygodny, i to o dwa rzędy bardziej, niż głosiła pierwsza wersja.
Założenie izotropii jest podtrzymywane przez dynamikę — nie dlatego, że
precesja zachodzi wokół \(\hat{\mathbf L}\), i nie dlatego, że nie ma kanału
strat, tylko dlatego, że **kanał strat robi się szybki dopiero tam, gdzie para
już nie przebywa**.

Warto dodać, że gdyby dyssypacja istniała, i tak nie dałaby odpychania w
orto: układ dążyłby do **minimum** energii dipolowej, a konfiguracja
minimalna jest przyciągająca w obu kanałach.

*I stąd wniosek, który obala wcześniejszą wersję tej sekcji.* Po uśrednieniu
czynnik zależy już tylko od \(\boldsymbol\mu_1\!\cdot\!\boldsymbol\mu_2\) i
\((\boldsymbol\mu_1\!\cdot\!\hat{\mathbf L})(\boldsymbol\mu_2\!\cdot\!\hat{\mathbf L})\).
Dla izotropowo losowanych momentów o ustalonym
\(\cos=\hat{\boldsymbol\mu}_1\!\cdot\!\hat{\boldsymbol\mu}_2\) zachodzi
\(\langle(\hat{\boldsymbol\mu}_1\!\cdot\!\hat{\mathbf L})(\hat{\boldsymbol\mu}_2\!\cdot\!\hat{\mathbf L})\rangle=\cos/3\)
(zmierzone na \(2\cdot10^6\) próbkach: \(-0{,}3001\) wobec \(-0{,}3000\)
przy \(\cos=-0{,}9\), i tak dalej), więc

\[
\langle\text{czynnik}\rangle=-\tfrac{\cos}2+\tfrac32\cdot\tfrac{\cos}3=0
\quad\textbf{tożsamościowo,}
\]

**dla każdego \(\cos\), czyli dla obu kanałów.** Człon dipolowy ma zerową
wartość oczekiwaną i **nie rozdziela para od orto**. Potwierdzone pomiarem
(4 ziarna × 15 trajektorii): para \(-0{,}061\pm0{,}219\) keV, orto
\(+0{,}157\pm0{,}147\) keV, różnica \(-0{,}217\pm0{,}263\) keV — wszystko
zgodne z zerem, znaki mieszane wewnątrz obu kanałów.

Wcześniejsze zdanie tego dokumentu, że „kanały rozdzielają się przy
promieniu terminalnym", było więc **nieprawdziwe**: ten model nie rozróżnia
kanałów w energii anihilacji nigdzie.

*Co zostaje.* Człon jest przy promieniu terminalnym **duży** i wnosi realny
**rozrzut**: rzędu \(\pm1{,}5\) keV na pojedynczą trajektorię przy zerowej
średniej. Przepowiednia, którą wspiera, to zatem linia **poszerzona**, nie
przesunięta — szerokość rzędu \(0{,}15\%\) z \(511\) keV, wobec
\(0{,}27\%\), o które wiązanie przesuwa jej środek.

*Zastrzeżenie skali.* Klasyczna energia wiązania zjadłaby całe \(W\) przy
\(a=r_e/4=0{,}70\) fm. Podłoga modelu leży \(274\times\) powyżej tej
granicy, więc problem nie jest osiągany — ale to niezależne potwierdzenie, że
bariera Comptona kończy trajektorię na długo przed tym, jak opis klasyczny
zjadłby masę spoczynkową pary.

#### Domknięcie sufitu: emisja w oknie zamiast chwilowego kopnięcia

Sufit opisany niżej wykluczał emisję w całej dziedzinie modelu, a obie
oczywiste naprawy (zatrzymanie hazardu, przeskalowanie kwantu) zostały
zmierzone i odrzucone. Domyka go zmiana **rodzaju** emisji, nie jej rozmiaru.

Zdarzenie nie przekazuje już fotonu w jednej chwili: **otwiera okno** długie
na jeden okres orbitalny, przez które kwant jest wypłacany w przyrostach. W
żadnej chwili nie ubywa więcej niż różniczka, a przez okno **pozycje
ewoluują**, więc energia pochodzi z orbity — ze studni potencjału — a nie z
zamrożonej chwilowej energii kinetycznej. To ta sama zmiana rodzaju, na
której opiera się ścieżka sekularna, wykonana na trajektorii rozdzielczej.
Jeden okres, bo fotonu o częstości \(\omega\) nie da się złożyć z krótszego
ciągu falowego.

*Dlaczego jest dokąd oddać energię.* Przy \(n=1\) orbita idzie
\(105{,}8\to35{,}3\) pm, wciąż \(184\times\) powyżej bariery Comptona; przy
\(n=0{,}25\) nowe \(a'=0{,}735\) pm, \(3{,}8\times\) nad barierą.

*Moment pędu wychodzi za darmo.* Przyrosty idą przez to samo
`applyStochasticDipolePhoton`, które skaluje **długość** pędu względnego,
zachowując kierunek. Daje to \(dE/dL=v^2/(v_tr)\), czyli dokładnie
\(\omega\) dla orbity kołowej — ten stosunek unosi promieniowanie
obrotowego dipola E1 — i mieści się w \(0{,}25\%\) od niego przy zmierzonej
medianie mimośrodu emisji \(0{,}05\).

*Okno reguluje się samo.* Tempo nominalne \(\hbar\omega/(2\pi/\omega)\)
rośnie jak \(\omega^2\), a energia kinetyczna, z której musi być zapłacone,
tylko jak energia wiązania — więc głęboko w kolapsie tempo nominalne
wyprzedza rezerwuar. Zmierzone przed dodaniem ogranicznika: odrzucane
przyrosty sięgały średnio \(2{,}5\times\), a maksymalnie \(13{,}7\times\)
dostępnej energii. Ograniczenie przyrostu do jednej dziesiątej **niezmiennika
CM** \(W-(m_1+m_2)c^2\) (nie laboratoryjnej energii kinetycznej — ta różnica
sama zostawiała \(1{,}3\%\) odmów) po prostu **wydłuża tam okno**: suma
wypłat to nadal dokładnie jeden kwant, a odmowa staje się niemożliwa.
Zweryfikowane przy sztucznie zawyżonym hazardzie
(\(10^5\times\)): **3155 wypłat, zero odmów**.

*Koszt, wypowiedziany wprost.* Trajektoria przestaje być zachowana dokładnie
między fotonami — w oknie działa siła ciągła. Współczynnik wypełnienia
ogranicza to do \(1{,}5\cdot10^{-6}\): jeden foton na \(669\,000\) obiegów,
okno długie na jeden obieg. **\(99{,}99985\%\) trajektorii zachowuje dokładną
zachowawczość**, na której ten tryb stoi.

Nakładanie się okien **akumuluje**, nie nadpisuje. Przy prawdziwym hazardzie
to zdarzenie o prawdopodobieństwie \(1{,}5\cdot10^{-6}\) i polityka nie ma
znaczenia, ale nadpisywanie porzucałoby niezapłaconą resztę, co byłoby
wyciekiem energii, a nie wyborem modelowym.

Zachowanie produkcyjne jest **bit-identyczne**: przy fizycznym hazardzie
ścieżka mechaniczna nie wyzwala fotonu w żadnym z przebiegów, jakie ten
projekt uruchamia, więc okno nie ma się kiedy otworzyć. Domyka natomiast
mechanizm, który dotąd po cichu gubił każdy foton, jaki próbował paść.

#### Zgodność harmonik z progami energetycznymi

Dwie ścieżki emisji mają **różne** progi kinematyczne, i to jest poprawne —
różnią się rodzajem emisji, nie jakością niezmiennika.

| ścieżka | niezmiennik | warunek |
|---|---|---|
| sekularna (`crem_collapse.hpp`) | \(W=(m_1+m_2)c^2+\mu\varepsilon\) | \(W^2-2WE_\gamma>0\), czyli \(E_\gamma<W/2\approx511\) keV |
| mechaniczna (`applyStochasticDipolePhoton`) | \(W=\gamma_1m_1c^2+\gamma_2m_2c^2\) | \(W^2-2WE_\gamma\ge(m_1c^2+m_2c^2)^2\), czyli \(E_\gamma\le\) energia kinetyczna w CM |

Ścieżka mechaniczna zmienia wyłącznie prędkości — pozycje, a więc i energia
potencjalna, są w chwili emisji nietknięte. Zatem \(KE'=KE-E_\gamma\) i foton
**może** czerpać tylko z energii kinetycznej. Ścieżka sekularna operuje na
elementach oskulacyjnych, czyli na wielkości uśrednionej po orbicie, której
rezerwuarem jest energia wiązania; para związana ma \(W<(m_1+m_2)c^2\) z
definicji i emisja po prostu wiąże ją głębiej.

**Ten sufit realnie tnie.** Twierdzenie o wiriale daje energię kinetyczną
orbity kołowej równą energii wiązania, podczas gdy
\(\hbar\omega/E_{\rm wiąz}=2\sqrt{a_{Ps}/a}\) — dwa na promieniu Bohra pary
i gorzej w miarę zapadania. Ułamek **okresu**, na którym foton \(n\hbar\omega\)
w ogóle się mieści:

| \(e\) | \(n=1\) | \(n=5\) | \(n=16\) | \(n=27\) |
|---|---|---|---|---|
| 0 | **0** | 0 | 0 | 0 |
| 0,3 | **0** | 0 | 0 | 0 |
| 0,5 | 0,149 | 0 | 0 | 0 |
| 0,9 | 0,113 | 0,017 | **0** | 0 |
| 0,97 | 0,098 | 0,015 | 0,003 | 0,001 |

Orbita kołowa lub słabo mimośrodowa **nie może wyemitować podstawowej nigdzie
na swoim okresie**.

**Dwie oczywiste naprawy są obie błędne** — stąd ten opis zamiast łatki.
Zachowanie zbankowanego hazardu i ponowienie próby w późniejszym kroku
**zakleszcza się** dla \(e\lesssim0{,}3\): hazard rośnie bez ograniczeń i foton
nie pada nigdy. Rozluźnienie bramki do warunku sekularnego **nie jest
ujednoliceniem**, tylko przeniesieniem warunku między dwoma różnymi rodzajami
emisji.

Zostaje realne ograniczenie dziedziny: **prescription \(\hbar\omega_{\rm orb}\)
i model emisji jako chwilowego kopnięcia są wzajemnie niespójne na orbicie
bliskiej kołowej**, a pogodzenie ich wymaga innego modelu emisji, nie innego
progu. Produkcja tego nie dotyka: praktycznie każdy foton pada ścieżką
sekularną (hazard mechanicznej to \(\sim5{,}5\cdot10^{-5}\) na obieg, a
przebieg `--mode 2` z instrumentacją w tym miejscu zanotował **zero
wywołań**).

#### Tryb wzbudzony: `--level n` i kaskada po drabinie Bohra

Pomiar, który ten tryb uzasadnia. Na drabinie Bohra sufit kinematyczny emisji
przez chwilowe kopnięcie wynosi **dokładnie**

\[
\frac{\hbar\omega}{E_{\rm kin}}=\frac{2}{n}
\]

(z twierdzenia o wiriale energia kinetyczna orbity kołowej równa się energii
wiązania), więc foton \(\hbar\omega\) nie mieści się przy \(n=1\) ani
\(n=2\), a mieści się od \(n=3\) w górę. Zmierzone w produkcji: model nigdy
nie opuszcza \(n\le1{,}09\) — 115 próbek z przebiegów związanych daje
\(n_{\rm eff}\in[0{,}184;\,1{,}090]\), a zjawiska niezwiązane (zderzenie,
rozpraszanie) nie dają stanów związanych w ogóle: **zero próbek z
\(n\ge2\)** w którymkolwiek z czterech zjawisk.

Stąd dwa wnioski. Sufit kinematyczny nie jest przypadkiem brzegowym, tylko
obowiązuje w **całej** dotychczasowej dziedzinie modelu. I różnica poziomów
nie miałaby tam na czym działać — dlatego wariant hybrydowy dostaje sens
dopiero razem z przygotowaniem pary wyżej na drabinie.

`--level n` przygotowuje parę na \(a_n=n^2a_{\rm pary}\). Pasmo prędkości
stycznej jest podawane w jednostkach prędkości kołowej **przy tej
separacji**, więc rozrzut \(L/(n\hbar)\) pozostaje niezmieniony i przesuwa
się wyłącznie poziom. Energia fotonu podąża wtedy za **odstępem poziomów**
\(\Delta E(n\to n-1)\), dopóki \(n\ge2\), a poniżej wraca do
\(\hbar\omega\).

To rozróżnienie ma znaczenie tylko przy małych \(n\): zmierzone
\(\Delta E/\hbar\omega\) wynosi \(1{,}0152\) przy \(n=100\),
\(1{,}0523\) przy \(n=30\), \(1{,}1728\) przy \(n=10\), ale **3,0000**
przy \(n=2\) — czyli zasada korespondencji działa asymptotycznie, a przy
dnie drabiny \(\hbar\omega\) myli się trzykrotnie.

Sprawdzone przeciw wzorowi analitycznemu
\(\Delta E/\hbar\omega=n(2n-1)/(2(n-1)^2)\):

| `--emission` | `poisson` | **Model produkcyjny.** Próg, na którym pada skwantowany foton. `poisson` losuje go z \(\mathrm{Exp}(1)\) — emisja spontaniczna z szumem śrutowym, zachowanie domyślne i fizyczne. `deterministic` ustawia próg na stałe \(1\), więc foton pada dokładnie wtedy, gdy orbita straciła ciągle jeden kwant (bankowany hazard to wypromieniowana energia w kwantach) — krok w stronę deterministycznego wyznaczania parametrów kwantowych, w którym energia fotonu wynika z orbity, a czas z jej ciągłej utraty energii. Zmierzone: \(\sigma/\mu\) spada z \(0{,}894\) do \(0{,}343\), średnia pozostaje bez zmian. Nie sprowadza rozkładu do odpowiedzi kontinuum — patrz sekcja „`--emission deterministic`" wyżej. |
| `--level` | \(n\) startowe | zmierzone | wzór |
|---|---|---|---|
| 2 | 2,0994 | 2,778 | 2,778 |
| 3 | 3,1492 | 1,806 | 1,806 |

**Kaskada schodzi po jednym poziomie na foton.** Zmierzone per emisja, a nie
per checkpoint (pierwsza wersja tego akapitu myliła jedno z drugim i podawała
liczby checkpointów jako liczby fotonów):

```
--level 3:  n=3.149 -> 2.130   dn=-1.019   dE
            n=2.130 -> 1.087   dn=-1.044   dE
            n=1.087 -> 0.634   dn=-0.452   hbar*omega
            n=0.634 -> 0.270   dn=-0.364   hbar*omega
            n=0.270 -> 0.080   dn=-0.191   hbar*omega
--level 2:  n=2.099 -> 1.053   dn=-1.046   dE
            (trzy dalsze emisje juz na hbar*omega)
```

Cała trajektoria to **pięć** fotonów przy `--level 3` i **cztery** przy
`--level 2`, z czego z \(\Delta E\) korzystają odpowiednio dwa i jeden.

Para **ląduje o jeden poziom niżej**, z dokładnością do 2–4%. To nie
przypadek, tylko tożsamość: skoro \(n=\sqrt{R/|E|}\), to \(|E|=R/n^2\) z
definicji, więc

\[
E'=-\frac{R}{n^2}-R\Big(\frac1{(n-1)^2}-\frac1{n^2}\Big)=-\frac{R}{(n-1)^2},
\]

czyli dokładnie poziom \(n-1\). Resztkowe 2–4% pochodzi z czynnika
\((\text{energyRatio})^{3/2}\), którym kod skaluje energię referencyjną
między checkpointami, oraz z relatywistycznej poprawki
\(\sqrt{W^2-2WE_\gamma}\approx W-E_\gamma-E_\gamma^2/2W\).

**Czego ten tryb nie robi.** Nie kwantyzuje energii orbitalnej: krok wynosi
jeden poziom, ale **odstęp od całkowitych poziomów Bohra jest zachowywany** —
kaskada startująca z \(n=3{,}149\) idzie przez \(2{,}13\) i \(1{,}09\), a
nie przez \(3\to2\to1\). Drabina ma poprawny **rozstaw**, lecz przesunięty
**offset**, równy części ułamkowej \(n\) w chwili przygotowania. Pełna
kwantyzacja byłaby inną zmianą, która przy domyślnym \(n=1\) zabroniłaby
emisji orbitalnej w ogóle, a więc usunęłaby z modelu badany kolaps.

**Koszt — i to nie jest jedno prawo potęgowe dla obu rodzin modeli.**
Pierwsza wersja tego akapitu podawała \(n^6\) bez zastrzeżeń, co jest
poprawne wyłącznie dla modeli ciągłych: tam \(dE/dt\sim a^{-4}\) przy
\(E\sim-1/a\) całkuje się do \(t\sim a^3\), czyli \(n^6\).

Domyślny tryb stochastyczny rządzi się hazardem, a hazard czyta także
**energię** fotonu. Przy \(P\sim a^{-4}\sim n^{-8}\) i
\(\hbar\omega\sim a^{-3/2}\sim n^{-3}\) czas oczekiwania szedłby jak
\(n^5\); gałąź odstępu poziomów dzieli następnie tempo przez
\(\Delta E/\hbar\omega=n(2n-1)/(2(n-1)^2)\). Zatem

\[
t\ \propto\ n^5\cdot\frac{n(2n-1)}{2(n-1)^2}\ \longrightarrow\ n^5 .
\]

Zmierzone bezpośrednio z tempa hazardu (bez szumu Poissona, w odróżnieniu od
czasu do pierwszego fotonu, którego odchylenie standardowe równa się średniej):

| n | 1/λ (s) | zmierzone × | \(n^5\) | \(n^6\) | wzór wyżej |
|---|---|---|---|---|---|
| 1,0497 | 2,352·10⁻¹⁰ | 1,0 | 1 | 1 | 1,0 |
| 2,0995 | 2,091·10⁻⁸ | **88,9** | 32 | 64 | **88,9** |
| 3,1492 | 1,032·10⁻⁷ | **438,9** | 243 | 729 | **438,9** |
| 4,1989 | 3,655·10⁻⁷ | **1554,2** | 1024 | 4096 | **1554,2** |

Wzór odtwarza pomiar co do ostatniej drukowanej cyfry; \(P\sim n^{-8}\)
potwierdzone niezależnie (zmierzony iloraz mocy między \(n=1\) a \(n=2\)
wynosi dokładnie \(256=2^8\)). Prawdziwa potęga asymptotyczna to więc
\(n^5\), nie \(n^6\) — dla `--level 8` daje to \(4{,}0\cdot10^4\times\), a
nie \(2{,}6\cdot10^5\times\), jak podawała pierwsza wersja. Komunikat
startowy liczy właściwy mnożnik dla aktywnego modelu reakcji i przypomina o
`--crem-wallclock-budget-s`.

`--level 1` jest domyślne i **bit-identyczne** z zachowaniem sprzed tej
zmiany, co jest sprawdzane przez porównanie `--diagnose` z budową HEAD.

#### Skala hazardu w stanie podstawowym wobec czasu życia para-Ps

Wielkość zmierzona przy okazji analizy progów i warta utrwalenia, bo nie
występowała dotąd nigdzie w tej dokumentacji.

Klasyczny hazard fotonowy w stanie podstawowym, liczony **wzorem samego
projektu** (`larmorOrbitAveragedPower`, \(|\ddot d|=|q_{\rm eff}|k|q_1q_2|/(\mu a^2)\),
co dla e⁺e⁻ daje \(q_{\rm eff}=e\), \(\mu=m_e/2\)) przy \(a=a_{Ps}\):

\[
P=1{,}16732208\cdot10^{-8}\,\mathrm{W},\qquad
\lambda=\frac{P}{\hbar\omega}=5{,}355\cdot10^{9}\,\mathrm{s^{-1}} .
\]

Wyrażone w naturalnej jednostce obie liczby są **wymierne i dokładne**:

\[
\hbar\omega\big|_{a_{Ps}}=\tfrac12\alpha^2m_ec^2,\qquad
\lambda=\tfrac13\,\frac{\alpha^5m_ec^2}{\hbar},\qquad
\Gamma_{\rm para}^{\rm LO}=\tfrac12\,\frac{\alpha^5m_ec^2}{\hbar}.
\]

Zmierzone w kodzie: \(0{,}50000000\) i \(0{,}33333333\). Stąd

\[
\frac{\tau_{\rm kl}}{\tau_{\rm para}}=\frac{1/2}{1/3}=\frac32
\quad\Longrightarrow\quad
186{,}74\ \mathrm{ps}\ \text{wobec}\ 124{,}49\ \mathrm{ps}.
\]

**Czym to jest.** Klasyczne tempo emisji E1 w stanie podstawowym i wiodące
tempo anihilacji parapozytonium mają **tę samą potęgę \(\alpha\)** i różnią
się tylko liczbą wymierną rzędu jedności. Po stronie klasycznej ma to
strukturalne wyjaśnienie — tempo radiacyjne orbity Bohra idzie jak
\(\alpha^3\omega\), a \(\omega=\alpha^2m_ec^2/2\hbar\), więc
\(\alpha^3\cdot\alpha^2=\alpha^5\). Po stronie QED \(\alpha^5\) bierze się
z zupełnie innego rachunku. **Zgodność współczynnika \(O(1)\) (1/3 wobec 1/2)
nie ma w tym modelu wyjaśnienia** i nie jest wyprowadzeniem anihilacji: to
różne procesy, jeden emituje kwanty \(13{,}6\) eV, drugi dwa po \(511\) keV,
czyli \(3{,}75\cdot10^4\) razy większe.

Warto zestawić z mierzoną medianą kolapsu w trybie stochastycznym,
\(149\) ps (10 ziaren), która leży **pomiędzy** tymi dwiema liczbami.

**Dlaczego żadne przeskalowanie kwantu nie zamyka sufitu.** Nasuwa się
propozycja, żeby w stanie podstawowym czerpać z energii całkowitej, a fotony
podstawowe uznać za anihilacyjne — co usunęłoby sufit kinematyczny opisany
wyżej. Rozstrzygnięcie wymaga rozdzielenia dwóch odczytań tej propozycji,
bo wyklucza je co innego.

*Odczytanie pierwsze: tempo anihilacji wyprowadzone z mocy klasycznej.*
Przy niezmienionym \(P\) podmiana kwantu z \(\hbar\omega\) na \(m_ec^2\)
dzieli tempo przez \(2/\alpha^2\approx3{,}8\cdot10^4\), dając \(7{,}0\) µs
zamiast \(125\) ps — cztery do pięciu rzędów wielkości za długo. Wykluczone
liczbowo.

*Odczytanie drugie: anihilacja we własnym, kwantowym tempie.* Tego pierwszy
rachunek nie dotyka i nie wyklucza — ale wyklucza je co innego: para
przestaje wtedy emitować orbitalnie, żyje do anihilacji i **kolaps, który
model bada, znika**. Jest to dokładnie ten sam powód, dla którego odrzucono
pełną kwantyzację poziomów, a nie nowy argument.

*Czy pomogłoby inne przeskalowanie kwantu.* Naturalne pytanie, skoro sufit
to \(\hbar\omega/E_{\rm kin}=2/n\), a przy \(n=1\) twierdzenie o wiriale
daje \(E_{\rm kin}=R\). Zmierzone:

| kwant | \(E_\gamma\) | czas życia | wobec \(124{,}5\) ps | sufit przy \(n=1\) |
|---|---|---|---|---|
| \(\hbar\omega\) | \(13{,}61\) eV | \(186{,}74\) ps | \(\times1{,}500\) | odrzucony |
| \(\tfrac23\hbar\omega\) | \(9{,}07\) eV | \(124{,}49\) ps | \(\times1{,}000\) | odrzucony |
| \(\tfrac12\hbar\omega\) | \(6{,}80\) eV | \(93{,}37\) ps | \(\times0{,}750\) | granicznie |
| \(\tfrac13\hbar\omega\) | \(4{,}54\) eV | \(62{,}25\) ps | \(\times0{,}500\) | mieści się |
| \(m_ec^2\) | \(511\) keV | \(7{,}01\) µs | \(\times5{,}6\cdot10^4\) | mieści się |

Wiersz \(\tfrac23\hbar\omega\) nie jest przypadkiem: hazard klasyczny
wynosi \(\tfrac13\alpha^5m_ec^2/\hbar\), a \(\Gamma_{\rm para}\) —
\(\tfrac12\), więc kwant mniejszy o \(\tfrac23\) odtwarza zmierzone tempo
**dokładnie**.

I tu jest sedno: **sufit wymaga \(E_\gamma<E_{\rm kin}=\tfrac12\hbar\omega\),
a zgodność z czasem życia wymaga \(E_\gamma=\tfrac23\hbar\omega\).** Te dwa
warunki są niespełnialne jednocześnie, bo \(\tfrac23>\tfrac12\) — rozmijają
się o czynnik \(\tfrac43\). Pierwszy kwant, który sufit faktycznie
przepuszcza, daje zgodność \(\times0{,}5\) lub gorszą, czyli gorszą niż
obecne \(\times1{,}5\). **Żadne przeskalowanie kwantu, w żadną stronę, nie
zamyka sufitu bez utraty tej zgodności** — a to jest mocniejsze stwierdzenie
niż samo „nie wolno podnosić energii", bo zamyka również kierunek przeciwny.

Część intuicji o rezerwuarze jest natomiast zgodna z kodem: ścieżka
sekularna już liczy \(W=(m_1+m_2)c^2+\mu\varepsilon\), czyli ma masę
spoczynkową w niezmienniku, i dlatego nie ma tam tego problemu.

Osobno warto pamiętać o ostrzeżeniu z sekcji o dipolach: podobieństwo liczb
\(2/3\) w sektorze dipolowym do \(2\gamma/3\gamma\) jest tam nazwane
numerologicznym zbiegiem okoliczności. Ta sekcja nie twierdzi nic więcej niż
zgodność skali — anihilacja pozostaje w CREM osobnym, jawnie oznaczonym
generatorem opartym na zmierzonych czasach życia.

#### Kwadrupol E2 i dlaczego znika dla pozytonium

Kwadrupol pary wynosi \(Q_{ij}=\kappa\,(3d_id_j-d^2\delta_{ij})\) z

\[
\kappa=\frac{q_1m_2^2+q_2m_1^2}{(m_1+m_2)^2},\qquad \mathbf d=\mathbf r_1-\mathbf r_2 .
\]

Dla pary **masowo symetrycznej** \(\kappa\equiv0\) — i to dokładnie, także w
arytmetyce zmiennoprzecinkowej, bo obie masy są tym samym literałem.
Zmierzone: \(\kappa=0\) dla e⁺e⁻, μ⁺μ⁻ i p+p̄, ale \(\kappa=-0{,}9989\,e\)
dla p+e⁻, czyli praktycznie pełna siła kanału. Powód jest geometryczny:
kanał E2 żyje z rozjechania się środka masy ze środkiem ładunku, a dla pary
cząstka–antycząstka te punkty się pokrywają.

Dlatego włączenie E2 do kwantu **nic nie zmienia dla pary domyślnej** (i dla
każdej symetrycznej), a dla p+e⁻ dokłada kanał, którego wcześniej po prostu
nie było. Zmierzone udziały w hazardzie (ortopozytonium/p+e⁻, seed 12345):
E2 = \(9{,}3\cdot10^{-25}\) (szum zaokrągleń) dla e⁺e⁻ wobec
\(1{,}09\cdot10^{-4}\) dla p+e⁻ — tam jest to **największy** kanał
podwiodący, około trzy rzędy wielkości powyżej M1.

#### Gdzie kończy się drabina multipoli

Następnym rzędem po E1 jest **M2/E3, a nie M1/E2**. Orbitalny dipol
magnetyczny niesie ten sam czynnik \(\kappa\) co kwadrupol, więc znika razem
z nim dla pary symetrycznej — ale znika też dla **każdej** pary, bo
\(\mathbf m=(\kappa/2)(\mathbf d\times\dot{\mathbf d})\) jest proporcjonalny
do zachowanego orbitalnego momentu pędu, czyli stały, więc \(\ddot{\mathbf m}=0\)
i nie promieniuje. (Człon M1 w sumie powyżej to moment **wewnętrzny**, inny
obiekt.)

To, które multipole elektryczne w ogóle istnieją, rozstrzyga **reguła
parzystości** — warto ją zapisać raz, bo odpowiada na całą nieskończoną
rodzinę, a nie na jeden szczebel. Pisząc \(\mathbf r_1=(m_2/M)\mathbf d\) i
\(\mathbf r_2=-(m_1/M)\mathbf d\), moment rzędu \(l\) niesie

\[
\kappa_l=\frac{q_1m_2^l+(-1)^lq_2m_1^l}{M^l},
\]

co dla pary **masowo symetrycznej** zwija się do \([q_1+(-1)^lq_2]/2^l\):
proporcjonalne do \(q_1+q_2=0\) dla \(l\) parzystego i do \(q_1-q_2=-2e\)
dla nieparzystego. Neutralna para symetryczna nie ma więc **żadnego**
parzystego multipola elektrycznego i ma **wszystkie** nieparzyste. Zmierzone
dla e⁺e⁻: \(\kappa_1=-1{,}000\,e\), \(\kappa_2=0\),
\(\kappa_3=-0{,}250\,e\), \(\kappa_4=0\), \(\kappa_5=-0{,}0625\,e\).

**E3 nie znika zatem dla pozytonium tak, jak znika E2** — jest tam wiodącą
poprawką do E1, razem z M2. Oba leżą na poziomie \(\beta^4\) względem E1
(\(E_l\) skaluje się jak \(\beta^{2(l-1)}\), \(M_l\) jak \(\beta^{2l}\)) —
o rząd dalej niż zwykłe \(\beta^2\), bo M2 płaci jednocześnie karę
„magnetyczną" i „o jedno l wyżej". Na promieniu Bohra pary \(\beta=\alpha\)
dokładnie, więc jest to \(2{,}84\cdot10^{-9}\).

Dla porównania, próg pomiarowy samego modelu: \(6{,}9\cdot10^{-5}\) (zgodność
kwadratury dalekiego pola z analitycznym Larmorem) i w najlepszym razie
\(1{,}9\cdot10^{-8}\) (zanieczyszczenie polem bliskim produkcyjnej powierzchni
kontrolnej). **M2 leży cztery rzędy poniżej szumu wszystkiego, z czym można by
je porównać** — i kosztowałoby kolejny komplet stencyli historii.

To jest różnica wobec E2, które dla p+e⁻ wynosi \(1{,}09\cdot10^{-4}\), czyli
**powyżej** tego progu, a nie promieniowało wcale. E2 było realną luką; M2 i
E3 nie są. Ten sam werdykt obejmuje E5 na poziomie \(\beta^8=8\cdot10^{-18}\),
o które nie trzeba już pytać.

Nie brakuje go też w **bilansie energii**: kwadratura Poyntinga w strefie
dalekiej nie jest rozwinięciem obciętym — kierunkowa retardacja w poprzek pary
zachowuje E3, M2, człony toroidalne i ich interferencję dokładnie
(`electromagneticFieldFluxRates`). Obcięcie dotyczy wyłącznie **analitycznej
sumy hazardu**, której zadaniem jest chwilowe tempo; strumień tego zadania
pełnić nie może, bo jest opóźniony o promień kontrolny, \(1{,}8\cdot10^{-13}\) s,
czyli dłużej, niż trwają niektóre przebiegi.

Moc E2 jest liczona tylko wtedy, gdy ktoś ją czyta: w modelu `automatic` (dla
bramki dominacji) i w trybie skwantowanym (dla hazardu). Modele ciągłe jej nie
liczą, bo wyrzuciłyby wynik razem z kosztem stencyli historii.

Zmierzone udziały M1 w hazardzie (seed 12345): ortopozytonium \(1{,}5\cdot10^{-24}\)
(przy zewnętrznym polu \(10^9\,\mu\)T rośnie do \(1{,}7\cdot10^{-14}\), co
potwierdza, że kanał jest podpięty i reaguje na wymuszenie precesji), zderzenie
czołowe \(6{,}0\cdot10^{-5}\). M1 jest więc wszędzie kanałem marginalnym.

Sprostowanie wcześniejszego zapisu: „zderzenie czołowe promieniuje w 100%
kanałem M1" było **artefaktem przyczynowości**, nie własnością fizyczną.
`radiatedEnergy` zbiera strumień z powierzchni kontrolnej oddalonej o
\(1{,}8\cdot10^{-13}\) s świetlnych, a ten przebieg trwa \(3{,}3\cdot10^{-17}\) s
— front E1 po prostu jeszcze tam nie dotarł, podczas gdy moc M1 jest liczona
natychmiastowo. Chwilowe moce dają M1/E1 = \(6\cdot10^{-5}\).

Estymator sekularny w `crem_collapse.hpp` pozostaje E1-only i inaczej być nie
może: całkuje hazard analitycznie przez pomijane orbity, mając do dyspozycji
wyłącznie elementy oskulacyjne, w których nie ma stanu spinowego pozwalającego
odtworzyć moc M1. W praktyce nie jest to luka — stany związane, które ta
ścieżka liczy, mają M1 dokładnie zero, bo oba momenty wchodzą do koherentnej
amplitudy M1 jako \(\mathbf m_1+\mathbf m_2\) i się znoszą.

Bramka regresyjna `quantized-radiation-drain` pilnuje tego wprost: w trybie
skwantowanym ciągła siła reakcji ładunku i drenaż rezerwuaru dipolowego muszą
być **dokładnie** zerowe. Test ma moc rozróżniającą — na budowie z przywróconym
ciągłym drenażem M1 czyta \(3{,}19\cdot10^{-13}\) i oblewa. Wcześniej żaden
wymuszany test nie wykonywał w ogóle ścieżki stochastycznej, więc pomyłka w
bramkowaniu (albo jej brak, czyli podwójne liczenie energii kanału) była
niewidoczna.

Niskoprędkościowym punktem odniesienia dla sprzężenia poruszającego się
ładunku z poruszającym się dipolem drugiej cząstki jest człon lagrangianu

\[
L_{q\mu}=q\,(\mathbf v_q-\mathbf v_\mu)\cdot\mathbf A_\mu,
\qquad
\mathbf A_\mu=\frac{\mu_0}{4\pi}\,
w(r)\frac{\boldsymbol\mu\times\mathbf r}{r^3}.
\]

Prędkość **względna** nie jest tu symetryzacją narzuconą po to, by siły się
równoważyły — jest ścisła w tym rzędzie. Dipol magnetyczny poruszający się z
\(\mathbf v_\mu\) niesie motoryczny dipol **elektryczny**
\(\mathbf p=(\mathbf v_\mu\times\boldsymbol\mu)/c^2\), a jego potencjał
skalarny wnosi do lagranżjanu ładunku dokładnie

\[
-q\varphi_p=-\frac{q}{4\pi\varepsilon_0}w(r)\frac{\mathbf p\cdot\mathbf r}{r^3}
=-q\,\frac{\mu_0}{4\pi}w(r)\frac{\mathbf v_\mu\cdot(\boldsymbol\mu\times\mathbf r)}{r^3}
=-q\,\mathbf v_\mu\cdot\mathbf A_\mu,
\]

gdzie skorzystano z \((\mathbf a\times\mathbf b)\cdot\mathbf c
=\mathbf a\cdot(\mathbf b\times\mathbf c)\) oraz
\(1/(4\pi\varepsilon_0c^2)=\mu_0/4\pi\). Zatem
\(q\mathbf v_q\!\cdot\!\mathbf A_\mu-q\varphi_p\) **jest**
\(q(\mathbf v_q-\mathbf v_\mu)\!\cdot\!\mathbf A_\mu\): to zwykłe
sprzężenie jednocząstkowe \(q\mathbf v_q\!\cdot\!\mathbf A\) z już
uwzględnionym własnym polem elektrycznym poruszającego się dipola, a nie jego
przybliżenie. Równość akcji i reakcji **wynika** stąd (lagranżjan zależy tylko
od \(\mathbf r\) i \(\mathbf v_q-\mathbf v_\mu\), więc niezmienniczość
translacyjna wymusza przeciwne wariacje), a nie jest powodem wyboru tej formy.
Ponieważ motoryczny dipol elektryczny jest tu już zawarty, członu tego nie
wolno sumować z `retardedElectricDipoleField`, opisującym tę samą fizykę przez
`state.firstElectricDipole`; i nie jest — `allExternalForces` używa jednego,
`retardedExternalForces` drugiego, a obie sumy sił są wariantami, nigdy nie
dodawanymi razem.

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

  **Domodelowanie tego brakującego orbitalnego momentu pędu — trzy próby,
  wszystkie odrzucone, udokumentowane, żeby następna próba nie musiała
  odkrywać tych samych ślepych zaułków od nowa.**

  *Próba 1 — dosłowne \(r\times p\) z próbkowanej anomalii prawdziwej.*
  Losowanie \(\nu\) jednostajnie w anomalii średniej (blok pominiętych
  orbit obejmuje ich nieznaną liczbę, więc "jednostajnie w czasie" jest
  jedynym nieuprzedzonym założeniem), rekonstrukcja \(r\) przez
  standardowe równanie Keplera \(r(E)=a(1-e\cos E)\), iloczyn wektorowy z
  już wylosowanym pędem fotonu. **Odrzucona przed wpięciem do stanu**:
  \(|r\times p_\gamma|\sim a\cdot(\hbar\omega_\gamma/c)=\hbar n(v_{orbita}/c)\)
  — standardowe tłumienie multipolowe wkładu orbitalnego względem skali
  \(\hbar\). Dla pozytonium \(v/c\sim10^{-3}\), czyli o **trzy rzędy
  wielkości za mało**, żeby domknąć lukę rzędu \(\hbar\), niezależnie od
  tego, jak dokładnie próbkuje się \(\nu\).

  *Próba 2 — przeskalowanie istniejącego mechanizmu ×2, skalibrowane do
  \(k(e)\) w średniej zespołowej.* Wymagana średnia usunięta wartość
  osiowa na foton, żeby odtworzyć niezależnie zweryfikowaną klasyczną
  formułę \(k(e)=-(1-e^2)/(2+e^2)\) (zgodność \(10^{-4}\) z mierzoną
  trajektorią), w zamkniętej postaci: \(\text{target}(e)=2|k(e)|
  \sqrt{1-e^2}/S(e)\). Policzone na pełnej tabeli \(S(e)\): stosunek
  target/(sam spin, \(\hbar/2\)) wynosi \(2{,}000\) przy \(e=0\), dryfując
  tylko do \(1{,}867\) przy \(e=0{,}97\) — niemal uniwersalna stała \(2\),
  nie poprawka zależna od \(e\). To nie numerologia: to podręcznikowy
  wynik dla jednostajnie wirującego klasycznego dipola — strumień momentu
  pędu dzieli się po połowie między "spin" i "orbitalny" wkład fali
  kołowo spolaryzowanej.

  **Zaimplementowana, zmierzona, wycofana**: na \(10\) niezależnych
  ziarnach \(|L|\) po pierwszym fotonie wyszło \(2{,}3\)–\(4{,}5\times\)
  celu klasycznego — GORZEJ niż niepoprawiona wersja (dla ziarna \(42\):
  z \(1{,}84\times\) na \(2{,}97\times\)). Przyczyna znaleziona i
  udowodniona ściśle (prawo cosinusów dla odejmowania wektorów):
  \[
  |L_{po}| = |L_{przed}|\sqrt{1+x^2-2xy},\qquad
  x=\frac{R}{|L_{przed}|},\quad y=\cos\theta\cdot h\in[-1,1].
  \]
  \(R=2\hbar/\mu\) jest, z tej samej konwencji Bohra/SED, na której stoi
  własny warunek początkowy tego modelu (\(L_{początkowe}=\hbar/\mu\)),
  niemal dokładnie DWUKROTNOŚCIĄ własnej charakterystycznej skali momentu
  pędu modelu — więc \(x\) wychodzi \(1{,}8\)–\(2{,}2\) na KAŻDEJ
  sprawdzonej trajektorii, nie przez przypadek, tylko z konstrukcji.
  Przy \(x\geq2\) minimum \(\sqrt{1+x^2-2xy}\) po WSZYSTKICH możliwych
  \(y\) wynosi \(|x-1|\geq1\): zmniejszenie wielkości jest geometrycznie
  NIEMOŻLIWE dla każdego pojedynczego losowania fotonu, nie tylko mało
  prawdopodobne. Żadne jednorodne przeskalowanie pojedynczego kopnięcia
  rzędu \(\hbar\) nie może tego ominąć, bo moment pędu pozytonium z
  konstrukcji modelu żyje w skali \(\hbar\) przez całą trajektorię, nie
  tylko blisko pierwszego fotonu.

  *Próba 3 — rozwiązanie ścisłe (nie w średniej), dla magnitudy \(R\)
  wzdłuż JUŻ wylosowanego kierunku fotonu, trafiającej dokładnie w
  `classicalAngularMomentumMagnitude` dla TEGO jednego zdarzenia.*
  \(|\mathbf L_{przed}+R\hat{\mathbf n}|^2=\text{cel}^2\) to kwadrat
  względem \(R\), z rozwiązaniem rzeczywistym tylko gdy
  \(\text{cel}\geq|L_{przed}|\sin\theta\) — geometryczna podłoga
  osiągalna przez ruch wzdłuż USTALONEGO kierunku. Sprawdzone na tych
  samych \(10\) ziarnach, każdy z własnym wylosowanym \(\theta\):
  **\(10/10\) NIEOSIĄGALNYCH** (cel zawsze poniżej podłogi, bo \(k(e)\)
  typowo żąda spadku wielkości o \(\sim45\)–\(50\%\) na jednym fotonie —
  te same duże, dyskretne skoki, dla których i energia skacze o rząd
  wielkości — a tylko \(\approx19\%\) kierunków ważonych
  \((1+\cos^2\theta)\) leży wystarczająco blisko osi, żeby taki spadek w
  ogóle był osiągalny przy jakiejkolwiek magnitudzie).

  **To wynik rozstrzygający.** Nie chodzi o dobór właściwego \(R\) —
  ścisłego czy przybliżonego. Wektor "skalar razy ustalony kierunek" ma
  JEDEN stopień swobody, a musi jednocześnie spełnić dwa na ogół
  niezgodne wymagania: trafić we właściwą KOŃCOWĄ wielkość (z \(k(e)\))
  i leżeć w kierunku wyznaczonym przez fizykę energii/polaryzacji fotonu,
  niezwiązaną z wymaganą redukcją \(L\). Domknięcie tej luki naprawdę
  wymagałoby dodatkowego stopnia swobody w kopnięciu — nowej fizyki, nie
  lepiej dobranej stałej. W tym momencie kod wracał do pierwotnej,
  niepoprawionej postaci (sam spin, \(\hbar\)) — **ale ten dodatkowy
  stopień swobody znalazł się, patrz próba 4 niżej.**

  **Ta sama analiza zastosowana do AKTYWNEGO mechanizmu (nie tylko do
  odrzuconej wersji ×2) ujawnia, że i on siedzi na tej samej granicy.**
  \(x=\hbar/(\mu L_{przed})\) mieści się w paśmie \(0{,}90\)–\(1{,}08\) na
  \(10\) niezależnych ziarnach (łagodniejsza wersja tego samego "limitu
  kwantowego", bo \(\hbar\), nie \(2\hbar\)) — a przy \(x\approx1\) i
  \(\langle y\rangle=1/2\) prawo cosinusów daje \(\langle L_{po}/
  L_{przed}\rangle=\sqrt{1+1-1}=1{,}0\): **praktyczny brak sekularnego
  spadku \(|L|\) na pojedynczym fotonie, tylko losowe błądzenie**.
  Zmierzone na tych samych \(10\) ziarnach: średnia \(L_{po}/L_{przed}=
  0{,}9724\) (zgodność z przewidywaniem \(1{,}0\) w granicach fluktuacji
  próby), \(3/10\) ziaren dało WZROST \(|L|\) zamiast spadku, żadne nie
  zbliżyło się do celu \(k(e)\) (`L_po/L_klasyczne` zawsze
  \(1{,}2\)–\(2{,}9\times\)).

  *Czy błądzenie się samokoryguje na dłuższej trajektorii? Sprawdzone
  wprost, nie założone.* Ziarno \(50\): trzy kolejne fotony w tym samym
  checkpoincie —
  \[
  \begin{array}{lccc}
   & L_{po}/L_{przed} & L_{po}/L_{klasyczne} & e_0^2\\
  \text{foton 1} & 1{,}673 & 3{,}34 & 0{,}0286\\
  \text{foton 2} & 0{,}997 & 4{,}40 & 0{,}0000\\
  \text{foton 3} & 0{,}667 & 3{,}51 & 0{,}0000
  \end{array}
  \]
  Przestrzelenie NIE zanika przy kolejnych losowaniach (\(3{,}34\to
  4{,}40\to3{,}51\times\)) — błądzi wokół wartości kilkukrotnie za dużej,
  bez tendencji powrotu do \(k(e)\).

  **Znaleziony konkretny skutek uboczny: sztuczna cyrkularyzacja.**
  \(e^2=\max(0,\,1+2EL^2/k_{Coul}^2)\) (\(E<0\)) — gdy \(L\) zostaje za
  duże (norma, nie wyjątek, bo sam spin z konstrukcji łapie tylko połowę
  \(k(e)\)), \(|2EL^2|\) przewyższa \(k_{Coul}^2\) i wyrażenie ucina się do
  DOKŁADNIE zera: orbita staje się sztucznie kołowa, nie z fizyki, tylko z
  tego, że \(L\) nie skurczyło się wystarczająco. Sprawdzone systematycznie
  na \(29\) parach kolejnych fotonów (\(24\) różne ziarna):
  \[
  \begin{array}{lcc}
   & e_0^2{=}0 \text{ przy nast. fotonie} & e_0^2{\neq}0\\
  \text{przestrzelone } (L_{po}/L_{klas}{>}1) & 28 & 0\\
  \text{niedostrzelone } (L_{po}/L_{klas}{\leq}1) & 0 & 1
  \end{array}
  \]
  **\(29/29\) zgodnych z przewidywaniem, zero wyjątków** — korelacja bez
  rozmytego pogranicza, bo mechanizm energii/hazardu systematycznie
  zwiększa \(|E|\) niezależnie od \(L\), więc gdy \(L\) zostanie za duże,
  \(e^2\) nie tylko lekko spada poniżej zera, tylko głęboko — obcięcie
  uruchamia się solidnie, nie na granicy. Skoro przestrzelenie jest
  regułą (nie wyjątkiem), **niemal każda trajektoria, która przeżyje do
  drugiego fotonu, ma ten i kolejne fotony sztucznie wymuszone na
  orbicie kołowej** — a wyzerowany `eccentricitySquaredHere` wchodzi z
  powrotem jako baza zarówno dla `classicalEccentricitySquared`
  (diagnostyka), jak i dla próbkowania harmoniki
  (`eccentricOrbitHarmonicNumber(e{=}0,\cdot)` daje harmonikę \(1\) z
  pewnością) — możliwa pętla samopodtrzymująca, w której raz wyzerowany
  mimośród nigdy nie ma szansy wrócić, tłumiąc strukturę harmoniczną,
  którą `CREM_HARMONIC` (punkt K) miał wprowadzać.

  *Próba 4 — rozdzielić wielkość od kierunku, zamiast prosić jeden
  wektor o oba naraz. UDANA.* Diagnoza wspólna dla prób 2 i 3: winowajcą
  nigdy nie była WIELKOŚĆ poprawki, tylko to, że jeden wektor "skalar razy
  ustalony kierunek" miał jeden stopień swobody na dwa zwykle sprzeczne
  zadania. Rozwiązanie: rozdzielić je między dwa już istniejące, osobno
  zweryfikowane źródła fizyki —
  \[
  \text{wielkość }|L| \leftarrow \texttt{classicalAngularMomentumMagnitude}
  \quad\text{(dokładna, ODE-scałkowana } k(e)\text{, już zweryfikowana}),
  \]
  \[
  \text{kierunek } \hat{\mathbf L} \leftarrow \text{ten sam wektor spinu}
  \text{ co dotąd, użyty TYLKO do wyznaczenia nowego kierunku}.
  \]
  Czy to podwójne liczenie, którego zakaz z próby 1–3 dotyczył? Nie —
  ten zakaz mówił o DODAWANIU \(k(e)\) NA WIERZCH pełnego, niezależnego
  wektora, który TAKŻE ustawia wielkość (kolizja). Tu \(k(e)\)
  **zastępuje** wielkość, dokładnie tak jak `photonEnergy` zastępuje
  klasyczne tempo mocy, a spin dostaje rolę (kierunek), której formuły
  Petersa–Mathewsa z definicji nigdy nie opisywały (są tylko o
  wielkości) — więc żadna fizyka nie jest liczona dwa razy.

  **Zaimplementowane i zweryfikowane, nie tylko wyprowadzone.** Po drodze
  znaleziony i naprawiony błąd jednostek (`classicalAngularMomentumMagnitude`
  jest już w jednostkach specyficznych — pierwsza wersja dzieliła przez
  `reducedMass` drugi raz, dając wynik zawyżony o \(\sim10^{30}\); ujawnione
  natychmiast przez rozjazd z sąsiedniej kolumny diagnostycznej w tym samym
  wierszu logu). Po poprawce, na tych samych \(10\) ziarnach:
  `L_spec(magnitude=k(e))` zgadza się z `L_spec(classical k)` **dokładnie**,
  w każdym pojedynczym zdarzeniu, nie tylko średnio. Sztuczna cyrkularyzacja
  znika: pierwszy i drugi foton dają teraz płynnie malejący, NIEZEROWY
  mimośród (np. ziarno \(42\): \(0{,}063\to0{,}010\)), zero pojawia się
  dopiero przy ostatnim, najgłębszym zdarzeniu tuż przed kolapsem — tam,
  gdzie prawdziwa fizyka też przewiduje niemal-kołową orbitę, nie jako
  artefakt. Trajektorie idą teraz rutynowo \(3\)+ fotonów głęboko (wcześniej
  kończyły się po \(1\)), sięgając energii \(10^{14}\)–\(10^{15}\) razy
  większych niż start.

  Walidacja: `positronium_validation` \(33/33\). Partia PARA \(N=20\):
  \(16/20\) do granicy, \(4\) ucięte budżetem czasu, **\(0\) awarii
  numerycznych**, mediana \(117{,}5\) ps. Partia ORTO \(N=20\): \(14/20\)
  do granicy, \(6\) ucięte, **\(0\) awarii numerycznych**, mediana
  \(162{,}1\) ps — ten sam rząd wielkości co para, zgodnie z oczekiwaniem.
  Sprawdzone bezpośrednio pod kątem NaN/Inf w stanie fizycznym na
  dodatkowych \(14\) trajektoriach — jedyne znalezione "nan" to
  nieszkodliwy artefakt statystyczny (odchylenie standardowe z próby
  \(N{=}1\)), nie korupcja fizyki. Sprawdzony ręcznie przypadek brzegowy:
  `energyAfterKick` w mianowniku mogłoby teoretycznie zejść do zera, gdyby
  `cmEnergyKick` było ujemne i większe co do wartości bezwzględnej od
  `photonEnergy` (zaobserwowane raz: `cmEnergyKick=-2,33e-23J` w ziarnie
  \(50\)) — w praktyce, na wszystkich przetestowanych trajektoriach, nigdy
  nie doprowadziło to do degeneracji.

  **Zewnętrzna zgodność, nie tylko wewnętrzna spójność**: mediana czasu
  kolapsu po tej poprawce (\(117\)–\(162\) ps w zależności od kanału)
  wraca w okolice wartości cytowanej niżej z historii tego pliku sprzed
  zastąpienia formuły `k` mechanizmem czysto spinowym
  (\(147{,}8\)–\(151{,}6\) ps) — niezależne potwierdzenie, że naprawa
  odtwarza fizykę, która była już wcześniej (osobno) zmierzona jako
  poprawna, a nie tylko wewnętrznie zgodny nowy wynik.

  Zmierzone empirycznie (dawniej, o promocji samego mechanizmu spinowego
  ponad zamrożone `k` — kontekst historyczny, patrz próba 4 wyżej dla
  bieżącego stanu): partia produkcyjna (o-Ps, \(80\) trajektorii, dwa
  ziarna) dała \(1\) awarię numeryczną wobec \(0\) awarii, ale \(3\) ucięć
  budżetem czasowym dla tej samej wielkości partii pod poprzednim kodem.
  Czas kolapsu ledwo drgnął (mediana orto \(130\)–\(155\) ps w kilku
  ziarnach, wobec \(147{,}8\)–\(151{,}6\) ps zmierzonych wcześniej z
  formułą `k`) — czas kolapsu wyznacza głównie całka energii/hazardu z
  punktów A–B, której ten mechanizm nie dotyka; mimośród (teraz naprawdę
  osiągający wartości jak \(e^2\approx0{,}9\)) wpływa na to, KTÓRY warunek
  wyjścia trajektoria trafi i jak szybko. Ta jedyna awaria doczekała się
  własnego, pełnego śledztwa — patrz punkt H.

**E2. Inne charakterystyki fotonu — czy poza spinem/momentem pędu coś jeszcze
wymaga zbadania.** Przegląd: energia (\(\hbar\omega\), `CREM_HARMONIC`),
pęd liniowy (odrzut, człon krzyżowy \(v_{cm}\cdot p\)), moment pędu (punkt
E), kierunek emisji (\((1+\cos^2\theta)\)), polaryzacja kołowa (Stokes V/I)
i harmonika — wszystkie już zweryfikowane w tej sekcji. Jeden kandydat
znaleziony i sprawdzony: **czy brak transformacji Lorentza między ramką
spoczynkową pary a laboratorium przy kolejnych fotonach kaskady psuje
fizykę.**

*Skala zmierzona wprost, nie oszacowana.* \(355\) zdarzeń fotonowych,
\(\sim200\) ziaren: \(v_{cm}/c\) rośnie systematycznie z głębokością
kaskady (mediana \(0{,}002\%\to0{,}014\%\to0{,}20\%\) dla fotonu \(1\to
2\to3\)), z sufitem \(\sim0{,}35\%\) (globalne maksimum, ziarno \(232\)) —
rzadkie głębokie kaskady (foton \(4\)+, \(\approx1\%\) trajektorii) nie
rosną dalej bez ograniczeń, bo sama fizyka (bariera Comptona) kończy
trajektorię, zanim zdąży się skumulować więcej kopnięć. \(0{,}35\%\) to
\(\sim5\) rzędów wielkości więcej niż sprzężenie dipol-dipol (\(17\) sond
w tej samej sekcji) — pierwsze wrażenie: to musi być istotne.

*Sprostowanie po dokładniejszym wyprowadzeniu — pierwszy alarm był
nietrafiony.* Wielkości orbitalne (`elements.specificEnergy`,
`specificAngularMomentum`, `period`) opisują ruch WZGLĘDNY, niezmienniczy
względem doładowania całego układu (nawet nierelatywistycznie, na czym
stoi mechaniczny rdzeń CREM) — nie zależą od \(v_{cm}\) w ogóle.
`photonEnergy`/`photonDirection` są więc liczone poprawnie w chwilowej
ramce spoczynkowej pary (S'), dokładnie tam, gdzie cała reszta już
zweryfikowanej fizyki (punkty A–E) żyje. Odrzut jest liczony w S' i
dodawany newtonowsko do zgromadzonego \(v_{cm}\) — poprawka od użycia
właściwego (relatywistycznego) składania prędkości zamiast newtonowskiego
jest rzędu \((v_{cm}/c)^2\sim1{,}2\times10^{-5}\) przy skrajnym
\(v/c=0{,}35\%\), NIE rzędu \(v/c\) — nieistotne. Sprawdzone wprost w
kodzie: `photonEnergy` z tego mechanizmu nigdy nie trafia do żadnego
histogramu/wykresu jako wielkość "obserwowana" — panele fotonowe pokazywane
użytkownikowi pochodzą z zupełnie innego, niezależnego generatora (prawdziwa
anihilacja \(2\gamma/3\gamma\), jawnie odseparowana, `positronium.cpp:1574`).
**Wniosek: żadne prawo zachowania ani istniejący wynik nie jest tu
naruszone** — pierwszy alarm mylił skalę przesunięcia \(v_{cm}\) względem
oryginalnego \(t{=}0\) z realnym błędem w księgowości fizyki.

*Prawdziwa, węższa potrzeba: nowa diagnostyka, nie poprawka.* Jedyne,
czego kod nie potrafił: podać, co zmierzyłby odległy, nieruchomy
obserwator dla energii/częstotliwości/kąta emisji głęboko-kaskadowego
fotonu — bo taka wielkość nigdy nie była liczona. Zaimplementowane jako
**nowa strona wykresów** (`--radiation-reaction stochastic`,
`distributions/1_3_*.pdf`), bez dotykania wewnętrznej fizyki:

```
E_lab = γ E' (1 + β cosθ')                    (Doppler)
θ_lab = kąt między p_lab a osią odrzutu β̂     (aberracja, pełny wzór wektorowy)
```

gdzie \(E',\theta'\) to już istniejące, zweryfikowane `photonEnergy`,
`photonDirection` w S', a \(\beta=v_{cm}/c\) brane PRZED kopnięciem tego
konkretnego fotonu (`LabFramePhoton`, `modules/crem_collapse.hpp`).
Zweryfikowane liczbowo (ziarno \(232\)): przy \(\beta=2\times10^{-4}\),
\(\Delta E/E\approx-3\times10^{-5}\) — właściwy rząd (\(\sim\beta\)),
zgodny znak z kątem. Trzy nowe panele
(`1_3_b_1_lab_frame_photon_energy.pdf`, `..._frequency.pdf`,
`..._angle.pdf`): sprawdzone wizualnie — sensowny rozkład energii
(\(N{=}43\), \(\langle E\rangle=826\) eV, od dziesiątek do
\(\sim3800\) eV, widoczna struktura kaskadowa), kąt od osi odrzutu
(\(N{=}28={=}43{-}15\) — dokładnie tyle, bo pierwszy foton każdej
trajektorii nie ma zdefiniowanej osi odrzutu, zgodnie z projektem).
Model ciągły (`--radiation-reaction coherent`, gdzie `labFramePhotons`
jest puste) renderuje bezpiecznie pusty wykres, bez awarii.
`positronium_validation` \(33/33\), bez regresji.

**E3. Audyt zasad zachowania na konkretnym zdarzeniu para-Ps.** Pełna
precyzja (\(17\) cyfr), stan przed/po pierwszym fotonie trajektorii
(ziarno \(50\), para-Ps), wszystkie trzy prawa sprawdzone niezależnie i
liczbowo, nie założone.

*Energia — zachowana do szumu numerycznego.* Przewidywana zmiana \(E\)
(z \(photonEnergy+cmEnergyKick\)) wobec rzeczywistej: różnica względna
\(1{,}4\times10^{-9}\).

*Pęd liniowy — zachowany do szumu numerycznego.* Przewidywane \(v_{cm}\)
po odrzucie (z pędu fotonu w S') wobec rzeczywistego: różnica względna
\(1{,}4\times10^{-9}\).

*Moment pędu — realna, ilościowa rozbieżność, świadoma konsekwencja
próby 4 (punkt E), teraz zmierzona po raz pierwszy.*
\[
|\mathbf L_{przed}-\mathbf L_{po}| = 1{,}4527\times10^{-4},\qquad
|\hbar\cdot h/\mu| = 2{,}3154\times10^{-4}\ (\text{stała}),
\]
stosunek \(0{,}627\), kąt między obydwoma wektorami \(78{,}5°\)
(przykładowa wartość TEGO zdarzenia — sprawdzone niżej na szerszej
próbie: to nie stała, wielkość waha się \(0{,}52\)–\(0{,}98\) i
\(7°\)–\(65°\) zależnie od zdarzenia). Usunięty z orbity moment pędu
NIE jest już równy spinowi pojedynczego fotonu — ani co do wielkości,
ani co do kierunku. Bilans wciąż się zamyka (\(L_{po}\) jest z definicji
tym, co zostaje), ale nazwa "moment pędu, jaki niesie foton" przestaje
dosłownie odpowiadać rzeczywistemu spinowi fotonu — dokładnie taki
kompromis, jaki próba 4 świadomie przyjęła (wielkość z \(k(e)\), kierunek
z wektora próbnego), tu po raz pierwszy skwantyfikowany na realnym
zdarzeniu.

*Zastrzeżenie metodologiczne: punkty 1–2 są częściowo tautologią.*
Sprawdzenie energii/pędu liniowego weryfikuje, że wzór księgowy się
zamyka — a on zamknie się ZAWSZE, bo `elements.specificEnergy` i
`centreOfMassVelocity` są aktualizowane właśnie TYMI wzorami. Nie
sprawdza, czy `photonEnergy` sam w sobie jest fizycznie poprawną
wartością dla PÓŹNIEJSZYCH fotonów w kaskadzie. To osobne pytanie —
zbadane niżej, i okazało się poważniejsze, niż sugerowałoby "znane,
udokumentowane przybliżenie brzegowe".

*Ten sam audyt powtórzony dla orto-Ps (ziarno \(50\)) — energia/pęd
zachowane z tą samą precyzją (\(1{,}361\times10^{-9}\) względnie, wobec
\(1{,}4\times10^{-9}\) dla para-Ps).* Sama trajektoria para vs orto (ten
sam seed): \(E_{przed}\) różni się o \(1{,}895\times10^{-8}\) względnie,
`photonEnergy` o \(1{,}551\times10^{-8}\) — dokładnie tego samego rzędu,
co już ustalone sprzężenie dipol-dipol (sondy 3–17, punkt L), rosnące z
głębokością trajektorii (sonda 16). Naprawa `hFraction` (punkt E4)
działa jednakowo w obu kanałach, bo mechanizm nigdzie nie odwołuje się
do flagi para/orto.

*Sprawdzone dokładniej, na żądanie: czy "\(63\%\)/\(78°\)" z ziarna
\(50\) było wartością TYPOWĄ, nie tylko zgodną między kanałami.*
Odpowiedź: **nie było** — to jeden punkt z szerokiego rozkładu, nie
charakterystyczna skala. \(14\) niezależnych zdarzeń fotonowych
(\(10\) ziaren, oba kanały, pełna precyzja):
\[
\begin{array}{cc|ccc|ccc}
\text{ziarno} & \text{zdarz.} & \text{stos. P} & \text{stos. O} & \Delta & \text{kąt P} & \text{kąt O} & \Delta\\
1 & 0 & 0{,}9795 & 0{,}9795 & {-}1{,}0{\times}10^{-10} & 7{,}35° & 7{,}35° & 0{,}000°\\
1 & 1 & 0{,}5244 & 0{,}5244 & 2{,}9{\times}10^{-10} & 61{,}28° & 61{,}28° & 0{,}000°\\
2 & 0 & 0{,}7342 & 0{,}7342 & 1{,}1{\times}10^{-10} & 64{,}60° & 64{,}60° & 0{,}000°\\
2 & 1 & 0{,}6897 & 0{,}6897 & {-}3{,}5{\times}10^{-8} & 13{,}88° & 13{,}88° & 0{,}000°\\
3 & 0 & 0{,}8638 & 0{,}8638 & {-}7{,}1{\times}10^{-11} & 54{,}32° & 54{,}32° & 0{,}000°\\
5 & 0 & 0{,}8793 & 0{,}8793 & 1{,}1{\times}10^{-10} & 25{,}22° & 25{,}22° & 0{,}000°\\
5 & 1 & 0{,}6620 & 0{,}6620 & 1{,}3{\times}10^{-7} & 10{,}08° & 10{,}08° & 0{,}000°\\
7 & 0 & 0{,}8835 & 0{,}8835 & 1{,}6{\times}10^{-12} & 28{,}64° & 28{,}64° & 0{,}000°\\
7 & 1 & 0{,}7244 & 0{,}7244 & {-}5{,}1{\times}10^{-9} & 8{,}34° & 8{,}34° & 0{,}000°\\
10 & 0 & 0{,}8949 & 0{,}8949 & {-}4{,}1{\times}10^{-12} & 25{,}37° & 25{,}37° & 0{,}000°\\
10 & 1 & 0{,}5165 & 0{,}5165 & {-}8{,}4{\times}10^{-10} & 52{,}23° & 52{,}23° & 0{,}000°\\
12 & 0 & 0{,}8590 & 0{,}8590 & 8{,}6{\times}10^{-11} & 24{,}50° & 24{,}50° & 0{,}000°\\
12 & 1 & 0{,}6182 & 0{,}6182 & 1{,}3{\times}10^{-8} & 9{,}20° & 9{,}20° & 0{,}000°\\
42 & 0 & 0{,}8903 & 0{,}8903 & 5{,}9{\times}10^{-11} & 31{,}02° & 31{,}02° & 0{,}000°
\end{array}
\]
**Stosunek \(|L_{usunięte}|/|L_{spin}|\) rozciąga się od \(0{,}517\) do
\(0{,}980\)** (nie stała \(\approx0{,}63\)), a **kąt od \(7°\) do
\(65°\)** (nie stała \(\approx78°\)) — silnie zależne od mimośrodu i
geometrii KAŻDEGO zdarzenia z osobna. Ale **niezależnie od tego, jak
bardzo dane zdarzenie odbiega od "typowej" wartości, para i orto zawsze
dają dokładnie ten sam wynik dla tego samego zdarzenia** — różnica kąta
\(0{,}000°\) w każdym z \(14\) przypadków, różnica stosunku na poziomie
szumu (\(10^{-8}\)–\(10^{-12}\)), nie ograniczona do jednego
przykładowego ziarna. To mocniejsze i bardziej precyzyjne potwierdzenie
niż pojedyncze zdarzenie z ziarna \(50\): rozbieżność momentu pędu
(próba 4) jest realną, zmienną cechą KAŻDEGO zdarzenia, ale jej
niezależność od kanału para/orto jest uniwersalna, nie przypadkiem
jednego pomiaru.

*Co właściwie napędza ten szeroki zakres — sprawdzone wprost, nie
sprzężenie dipol-dipol.* Rozbite na dwie zmienne wejściowe każdego
zdarzenia:
\[
\begin{array}{cc|ccc}
\text{ziarno} & \text{zdarz.} & \text{stos.} & \text{kąt} & \cos\theta,\,e^2\\
1 & 0 & 0{,}9795 & 7{,}35° & -0{,}7683,\ 0{,}0390\\
12 & 1 & 0{,}6182 & 9{,}20° & -0{,}9523,\ 0{,}0022\\
7 & 1 & 0{,}7244 & 8{,}34° & 0{,}9431,\ 0{,}0015\\
5 & 1 & 0{,}6620 & 10{,}08° & 0{,}9354,\ 0{,}0000\\
10 & 1 & 0{,}5165 & 52{,}23° & -0{,}2514,\ 0{,}0006\\
3 & 0 & 0{,}8638 & 54{,}32° & -0{,}2138,\ 0{,}0577\\
1 & 1 & 0{,}5244 & 61{,}28° & 0{,}0735,\ 0{,}0060\\
2 & 0 & 0{,}7342 & 64{,}60° & 0{,}0722,\ 0{,}0098
\end{array}
\]
**Kąt koreluje bezpośrednio z \(|\cos\theta|\)** (kąt emisji fotonu
względem osi orbity): blisko osi (\(|\cos\theta|\to1\)) kąt
"usunięte-vs-spin" jest mały (\(7°\)–\(10°\)); blisko płaszczyzny
równikowej (\(|\cos\theta|\to0\)) kąt jest duży (\(52°\)–\(65°\)) — czysta
geometria: gdy spin fotonu jest niemal równoległy do \(L\), przeskalowanie
wielkości (\(k(e)\)) ledwo odchyla wynik od samego spinu; gdy niemal
prostopadły, to samo przeskalowanie daje duże odchylenie kątowe.
**Stosunek koreluje z mimośrodem \(e^2\)** (razem z \(\cos\theta\)) —
wyższe \(e^2\) daje na ogół wyższy stosunek, bo \(k(e)\) silniej
zmniejsza cel wielkości przy większym mimośrodzie. Obie te zmienne
(\(\cos\theta\), \(e^2\)) pochodzą z tego samego strumienia PRNG i tej
samej klasycznej mechaniki orbity — identycznej w para i orto (stąd
zgodność \(0{,}000°\) powyżej). **Sprzężenie dipol-dipol nie ma tu
żadnego udziału** — to dwa całkowicie odrębne źródła zmienności, na
zupełnie różnych skalach: geometria pojedynczego zdarzenia (\(\cos\theta,
e^2\)) napędza cały zakres \(0{,}52\)–\(0{,}98\)/\(7°\)–\(65°\)
MIĘDZY różnymi zdarzeniami; sprzężenie dipol-dipol odpowiada wyłącznie
za maleńką resztkową różnicę MIĘDZY kanałami para/orto dla TEGO SAMEGO
zdarzenia.

*Sprawdzone dokładniej, na żądanie, wobec aktualnego stanu modelu (po
próbie 4 i naprawie `hFraction`): sama skala tej resztkowej różnicy —
poprzednie \(10^{-8}\)–\(10^{-9}\) było zbyt wąskie.* Policzone na
tych samych \(14\) zdarzeniach, rozdzielone na pierwszy foton trajektorii
(zdarzenie \(0\)) i drugi w tej samej kaskadzie (zdarzenie \(1\)):
\[
\begin{array}{lc}
\text{zdarzenie \#0 (pierwszy foton), 8 przypadków} & 7{,}9\times10^{-10}\text{ -- }1{,}0\times10^{-8}\\
\text{zdarzenie \#1 (drugi foton tej samej kaskady), 6 przypadków} & 6{,}8\times10^{-8}\text{ -- }1{,}17\times10^{-6}
\end{array}
\]
**Pełny zmierzony zakres: \(7{,}9\times10^{-10}\) do \(1{,}17\times10^{-6}\)**
— drugie zdarzenie w tej samej kaskadzie jest systematycznie \(2\)–\(3\)
rzędy wielkości większe niż pierwsze. Zgodne z już ustaloną fizyką (sonda
16: różnica rośnie z liczbą upłynionych orbit/checkpointów przed danym
zdarzeniem) — ale teraz, po naprawie `hFraction` (punkt E4), kaskady
poprawnie sięgają głębiej, więc kolejne fotony dziedziczą pełną,
nieskompresowaną akumulację różnicy sprzed nich, i efekt jest wyraźniej
widoczny na tych samych, płytkich trajektoriach niż w chwili, gdy sonda
16 była mierzona. Poprawny opis: **\(10^{-10}\)–\(10^{-6}\), rosnące z
każdym kolejnym fotonem tej samej kaskady** — nie stały rząd wielkości
\(10^{-8}\)–\(10^{-9}\).

**E4. Skala błędu przybliżenia `hFraction`, zmierzona i naprawiona.** Kod już
miał komentarz ostrzegający: gdy próg pochłaniany jest głównie z hazardu
"przeniesionego" z poprzednich checkpointów (`hFraction` nasyca się do
\(1\)), foton jest przypisywany do POCZĄTKU bieżącego skoku zamiast do
właściwego miejsca — "przybliżenie brzegowe, nie dokładne". Nigdy nie
zmierzone ilościowo. W kaskadzie trzech kolejnych fotonów (ziarno \(50\),
ta sama trajektoria co E3) wszystkie trzy dostały **identyczną** energię
(\(3{,}9314042382173936\times10^{-18}\) J co do wszystkich \(17\) cyfr),
mimo że \(|E|\) urosło między nimi o rząd wielkości. Przyczyna
zlokalizowana: `hFraction` nasyca się do \(1\) już przy pierwszym
fotonie (próg \(1{,}6633\) był wielokrotnie większy niż `skipHazard`
\(0{,}0978\) tego jednego checkpointu), więc `sAtPhoton` przypina się do
\(\text{jumpParameter}\) — maksimum, jakie przybliżenie w ogóle
dopuszcza — dla wszystkich trzech, niezależnie od tego, jak daleko orbita
faktycznie zdążyła się zapaść między nimi.

*Skala, skalibrowana wprost z danych trajektorii (nie założona).*
Wyprowadzone \(\omega(E)=(2|E|)^{1{,}5}/\text{attractionParameter}\) z
`photonEnergyReference` checkpointu \(17\) (sprawdzone na checkpoincie
\(16\): zgodność do \(6\) cyfr), zastosowane do PRAWDZIWEGO stanu \(E\)
tuż przed fotonem \(2\) i \(3\):
\[
\begin{array}{lccc}
 & \text{u\.zyta } photonEnergy & \text{"poprawna" (z prawdziwego } E) & \text{niedoszacowanie}\\
\text{foton 2} & 3{,}93\times10^{-18}\text{ J} & 2{,}27\times10^{-17}\text{ J} & 83\%\\
\text{foton 3} & 3{,}93\times10^{-18}\text{ J} & 5{,}29\times10^{-17}\text{ J} & 93\%
\end{array}
\]
**Nie kilka procent — czynnik \(5{,}8\times\) i \(13{,}5\times\) za
mało.** Przybliżenie jest ograniczone z góry przez
\((1-\text{jumpParameter})^{-1}\) (\(\approx1{,}43\times\) przy
\(\text{jumpParameter}\approx0{,}3\)) niezależnie od tego, jak daleko
poza ten jeden checkpoint prawdziwy stan już się oddalił — więc błąd
**rośnie bez ograniczenia z każdym kolejnym fotonem tej samej
kaskady**, właśnie w reżimie (szybkie, wielofotonowe kaskady), który
próba 4 uczyniła powszechnym. Energia radiowana w takich kaskadach jest
więc systematycznie ZANIŻANA względem tego, co prawdziwa, ewoluująca
orbita powinna emitować — nieco spowalniając tempo, w jakim \(|E|\)
faktycznie się pogłębia w takich seriach, względem poprawnej fizyki.

**NAPRAWIONE i zweryfikowane, nie tylko zmierzone.** Diagnoza: dla
PIERWSZEGO fotonu w kaskadzie `period`/`eccentricityHere`/
`photonEnergyReference` (liczone raz na checkpoint, przed pętlą `while`)
SĄ bieżącym stanem, bo nic jeszcze go nie ruszyło — błąd dotyczy
wyłącznie fotonu DRUGIEGO i kolejnych w tej samej pętli, dla których
`elements.specificEnergy`/`specificAngularMomentum` już zostały
zaktualizowane przez poprzedni foton, a te trzy wielkości — nie.
Poprawka: flaga `cascadeStateAlreadyMoved`, ustawiana `true` zaraz po
pierwszej aktualizacji stanu w danej pętli; dla każdego kolejnego fotonu
`period` (przez `regularizedPeriod`), mimośród i
`photonEnergyReference` są przeliczane NA NOWO z bieżącego, już
zaktualizowanego stanu, z `energyRatio=1` (wartość dokładna nie wymaga
ekstrapolacji obwiedni). `orbitsToSkip`/`jumpParameter`/`skipHazard`
(opisujące, ile orbit obejmuje CAŁA całka hazardu tego checkpointu) NIE
są przeliczane — to osobna księgowość od tego, jakiej referencji
energii/mimośrodu używa POJEDYNCZY foton przy własnym losowaniu.

*Zweryfikowane na tej samej trajektorii (ziarno \(50\)).* Foton \(2\)
po poprawce: \(2{,}27463\times10^{-17}\) J — zgadza się co do
wyświetlanej precyzji z ręcznie wyprowadzoną "poprawną" wartością
(\(2{,}27\times10^{-17}\) J) sprzed naprawy. Foton \(3\): teraz
\(2{,}82\times10^{-16}\) J (jeszcze większy niż wcześniejsze ręczne
oszacowanie \(5{,}29\times10^{-17}\) J, bo referencja dla fotonu \(3\)
teraz poprawnie startuje z JUŻ poprawionego stanu po fotonie \(2\) —
błąd kaskadowo się kumulował, więc i poprawka kaskadowo się wzmacnia).
Energia/pęd liniowy nadal zachowane dokładnie (wzór księgowy niezmieniony,
zmieniła się tylko wejściowa wartość `photonEnergy`). Partie \(N=20\)
para/orto: **\(0/20\) awarii numerycznych w obu**, \(0\) ucięć budżetem
czasowym w obu (wcześniej \(4/20\) i \(6/20\) — trajektorie teraz kończą
się szybciej, bo poprawnie promieniują więcej energii w kaskadach,
zamiast być sztucznie spowalniane). `positronium_validation` \(33/33\).

**F. Co model świadomie zostawia otwarte.** (1) Orbitalny moment pędu
fotonu względem pary jako dosłowny \(\mathbf r\times\mathbf p_\gamma\) —
wciąż wymagałby anomalii prawdziwej, niedostępnej w reprezentacji samych
elementów oskulacyjnych. Praktyczny SKUTEK tego braku (wielkość \(L\)
systematycznie za duża, sztuczna cyrkularyzacja) jest już naprawiony
inaczej — rozdzieleniem wielkości i kierunku, nie rekonstrukcją \(r\)
(punkt E, próba 4). (2) Ścieżka mechaniczna
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

**I. Dalszy audyt spójności, na zadane pytanie "co jeszcze może być
niespójne" — jedno pytanie doprecyzowujące rozstrzygnięte wprost, jedno
znalezisko naprawione.**

*Czy dokładniej byłoby to opisać jako dipol elektryczny plus kwadrupol
magnetyczny?* Sprawdzone wprost, nie na pamięć: dla pary o równych masach
i przeciwnych ładunkach (e⁺e⁻, ale też μ⁺μ⁻) orbitalny moment dipolowy
magnetyczny \(\boldsymbol\mu_{orb}=\tfrac12\sum_iq_i(\mathbf r_i\times
\mathbf v_i)\) i orbitalny moment kwadrupolowy elektryczny
\(Q_{orb}=\sum_iq_i(3\mathbf r_i\mathbf r_i-r_i^2\mathbb 1)\) znikają
**dokładnie**, tożsamościowo — współczynnik przy obu sprowadza się do
\(q_1+q_2=0\) (kwadrupol) lub kombinacji znikającej dla \(m_1=m_2\)
(dipol magnetyczny), potwierdzone bezpośrednim rachunkiem. To oznacza, że
zwykle zakładana "następna poprawka po E1" — M1+E2, rzędu \(\beta^2\)
względem E1 — **znika tu tożsamościowo przez symetrię**, a prawdziwa
najniższa nietrywialna poprawka to rzeczywiście kwadrupol magnetyczny M2
(razem z oktupolem elektrycznym E3, tego samego rzędu) — więc tak,
**opisanie tego jako "dipol elektryczny plus kwadrupol magnetyczny" jest
fizycznie dokładniejsze** niż zwyczajowe "E1 plus poprawki rzędu
\(\beta^2\)". Oszacowane: M2/E3 są rzędu \(\beta^4\) względem E1 (płaci
podwójną karę — "magnetyczny" ORAZ "o jeden rząd \(l\) wyżej"), czyli
\(\sim3\cdot10^{-9}\) przy \(a_{Ps}\), rosnące do rzędu jedności dopiero
tuż przy barierze Comptona, gdzie reszta modelu już i tak przestaje
obowiązywać. Sam ten brak — nienaprawiony, bo prawdopodobnie pomijalny w
całym praktycznym zakresie modelu — nie był jednak celem tego punktu;
poszukiwanie go doprowadziło do czegoś poważniejszego.

*Prawdziwy wyciek energii, znaleziony przy okazji.* `deltaEnergyPerOrbit`
/`orbitalRadiatedEnergy` — prawdziwa, zmierzona strata energii z
JEDYNEGO mechanicznie całkowanego obiegu na każdym checkpoincie — nigdy
nie trafiała do `elements.specificEnergy` ani do `radiatedEnergyTotal`
dla modelu stochastycznego (dla gałęzi deterministycznej owszem: ten sam
pomiar zasila całą ekstrapolację obwiedni). Zmierzone bezpośrednio na
dwóch trajektoriach: ziarno \(42\) (płytka, \(37\) checkpointów) —
pominięta energia stanowiła \(2\cdot10^{-5}\) energii skredytowanej przez
fotony w całym przebiegu, pomijalne; ziarno \(107\) (ta sama trajektoria
co w punkcie H, spędzająca czas przy wysokim mimośrodzie) —
**\(38{,}5\%\)**, realne, niebagatelne złamanie zachowania energii, nie
błąd zaokrąglenia. Co więcej: to poprzedni fix strażnika (punkt H) czyni
ten wyciek łatwiejszym do osiągnięcia w praktyce — trajektoria, która
wcześniej po prostu padała przy wysokim mimośrodzie, teraz biegnie dalej
i akumuluje ten niekredytowany ubytek.

Naprawione w `crem_collapse.hpp`: `run.finalState.orbitalRadiatedEnergy`
jest teraz doliczane do `elements.specificEnergy` i `radiatedEnergyTotal`
raz na checkpoint dla modelu stochastycznego, DODATKOWO do (nie zamiast)
kredytów fotonowych — te dwa pokrywają rozłączne odcinki (ten jeden
zmierzony obieg kontra `orbitsToSkip` pominiętych po nim), więc to
domknięcie luki księgowej, nie podwójne liczenie, przed którym
wcześniejsze poprawki w tej sekcji musiały się bronić. Zweryfikowane: ta
sama trajektoria z punktu H (ziarno \(107\)) nadal kończy się w tym samym
\(172{,}78\) ps, przy wystarczającym budżecie zegara — uczciwe liczenie
prawdziwej energii kosztuje prawdziwy czas obliczeń (więcej checkpointów
przy wysokim mimośrodzie, bo więcej budżetu energii jest teraz poprawnie
wydawane na każdym z nich). Partia \(30\) trajektorii (ziarno \(99\)) bez
zmian: \(29/30\) dochodzi do granicy, \(0\) awarii numerycznych, mediana
\(100{,}653\) ps — identycznie jak przed tą poprawką. Czysta kompilacja,
`positronium_validation` 33/33.

**J. Skąd bierze się foton, i czy jego produkcja mogłaby być wyzwalana —
dwa pytania koncepcyjne, obie odpowiedzi sprawdzone, nie tylko
wyargumentowane.**

*Czy foton rodzi się w środku masy?* W samej architekturze modelu —
efektywnie tak, i to z konieczności: reprezentacja `crem_collapse.hpp`
(same elementy oskulacyjne — `specificEnergy`, `specificAngularMomentum`)
nie niesie żadnej informacji o pozycji, a mechanizm odrzutu
(`centreOfMassVelocity`, punkt D) explicite traktuje emisję jako
jednorodne kopnięcie całego układu, dokładnie równoważne przypisaniu
źródła do środka masy. Fizycznie to jest standardowe przybliżenie dipola
punktowego, i akurat dla pozytonium (obojętnego elektrycznie,
\(q_1+q_2=0\)) wygodne: moment dipolowy elektryczny takiego układu jest
**niezależny od wyboru początku układu współrzędnych**, więc "w środku
masy" nie jest ani lepszym, ani gorszym wyborem niż jakikolwiek inny, o
ile mówimy tylko o samym momencie dipolowym. Ale prawdziwym źródłem są
przyspieszające ładunki, poruszające się po ruchu WZGLĘDNYM, nie po ruchu
środka masy — jednorodnie poruszający się obojętny układ w ogóle nie
promieniuje E1 (ta sama tożsamość \(\sum m_i\mathbf r_i=0\), która
domknęła poprawkę momentu pędu w punkcie E). Więc ściślej: foton pochodzi
z oscylacji/rotacji separacji ładunków (orbity), a przypisanie go do
środka masy to wygodna konwencja matematyczna przybliżenia dalekiego
pola, nie twierdzenie o realnym, zlokalizowanym punkcie emisji — to ta
sama luka informacyjna (brak anomalii prawdziwej), która już wykluczyła
domknięcie orbitalnej części momentu pędu fotonu w punkcie E.

*Czy generacja jest losowa, czy dałoby się znaleźć wyzwalacz w
konfiguracji pola?* Tak jak zaimplementowano: czysto losowa (proces
Poissona, próg \(\mathrm{Exp}(1)\), zero warunku wyzwalającego związanego
z konfiguracją pola). W ortodoksyjnej (kopenhaskiej) elektrodynamice
kwantowej to nie jest uproszczenie modelu, tylko własność natury: emisja
spontaniczna jest tam nieredukowalnie losowa, bez ukrytej zmiennej, którą
dałoby się z zasady odnaleźć. Istnieje jednak alternatywa —
elektrodynamika stochastyczna (SED) — w której pole punktu zerowego
traktowane jest jako realne, klasyczne (choć nieznane w szczególe) pole
losowe; dla USTALONEJ konkretnej realizacji tego pola odpowiedź cząstki
byłaby w zasadzie w pełni deterministyczna, więc formalnie dałoby się
zdefiniować "wyzwalacz" (konkretna realizacja pola plus faza orbity
przekraczająca próg) — losowość przesuwałaby się tylko do warunku
początkowego samego pola, nie znikałaby. Ten wątek jest w tym repozytorium
już zbadany i zamknięty z wynikiem negatywnym (`--zpf`, patrz niżej
"Wynik eksperymentu z polem punktu zerowego"): sprzężenie klasycznego
pola stochastycznego do orbity nie dało samouzgodnionej równowagi —
pasmo rezonansowe nie robi praktycznie nic, szersze pasma tylko pompują
energię i powodują ucieczkę orbity, nigdy stabilizację. To nie jest
bezpośredni test "czy ZPF wyzwala pojedyncze zdarzenia emisji
deterministycznie", ale pokazuje, że próba samouzgodnionego domknięcia
tej pętli w tej architekturze już się nie powiodła — więc **domyślny** wybór
(losowość Poissona) jest zgodny z główną gałęzią fizyki, nie arbitralnym
uproszczeniem.

*Uzupełnienie: częściowy determinizm już w repozytorium, i czym dokładnie
jest.* Zdanie „obecny wybór to czysta losowość Poissona" opisuje wyłącznie
tryb domyślny. `--emission deterministic` ustawia próg emisji na dokładnie
\(1{,}0\) zamiast losować \(\mathrm{Exp}(1)\), i obejmuje **obie** ścieżki
emisji — orbitalną (`crem_trajectory.hpp`) i analityczny skip
(`crem_collapse.hpp`), a więc całość hazardu, nie jego część. Ponieważ sam
hazard jest deterministycznym funkcjonałem trajektorii, w tym trybie **moment
emisji nie niesie żadnej losowości**: \(n\)-ty foton pada dokładnie wtedy, gdy
scałkowany hazard osiąga \(n\).

Nie jest to jednak wyzwalacz w sensie SED i nie należy tego mylić. Zmienną,
która o emisji decyduje, jest tu **skumulowany hazard samej orbity**, a nie
realizacja pola zewnętrznego; determinizm dotyczy chwili emisji, nie jej
kierunku (ten nadal losuje `sampleRotatingDipolePhotonDirection`) ani warunków
początkowych. Wątek SED — czyli przeniesienie losowości do konkretnej
realizacji pola punktu zerowego — pozostaje zamknięty z wynikiem negatywnym,
niezależnie od tego trybu.

Głębsza reformulacja, w której deterministyczne byłyby także kierunek i widmo,
pozostaje otwartym, trudnym problemem badawczym.

*Jaka konfiguracja pól pozytonium odpowiada polu fotonu — sprawdzone
wprost, nie tylko przez analogię.* Pełne pole Liénarda–Wiecherta
(sekcja 5) rozpada się na człon prędkościowy (\(\propto1/R^2\), bliskie,
quasi-statyczne, księgowane osobno jako energia Schotta) i człon
przyspieszeniowy (\(\propto1/R\), dalekie, poprzeczne, propagujące) — to
ten drugi, w granicy korespondencji klasyczno-kwantowej, staje się
fotonem; to dokładnie ten sam człon, którego strumień Poyntinga daje
`orbitalRadiatedEnergy`. Sprawdzone dalej, konkretniej: czy stan
polaryzacji tego klasycznego pola dalekiego odpowiada rozkładowi
skrętności fotonu już zaimplementowanemu w punkcie E. Policzone
numerycznie (niezależnie od wyprowadzenia w punkcie E, z surowych pól
\(\mathbf E\) wirującego dipola, nie z gotowej formuły): parametr Stokesa
\(V/I\) (stopień polaryzacji kołowej) w funkcji kąta \(\theta\) od osi
wirowania —

| \(\theta\) | \(V/I\) zmierzone | \(2\cos\theta/(1+\cos^2\theta)\) | różnica |
|---|---|---|---|
| \(0°\) | \(1{,}000000\) | \(1{,}000000\) | \(0\) |
| \(30°\) | \(0{,}989743\) | \(0{,}989743\) | \(-1{,}1\cdot10^{-16}\) |
| \(45°\) | \(0{,}942809\) | \(0{,}942809\) | \(-2{,}2\cdot10^{-16}\) |
| \(60°\) | \(0{,}800000\) | \(0{,}800000\) | \(-1{,}1\cdot10^{-16}\) |
| \(90°\) | \(0{,}000000\) | \(0{,}000000\) | \(-4{,}9\cdot10^{-32}\) |

zgodność do precyzji maszynowej, w całym zakresie kąta. To dokładnie ta
sama formuła \(2\cos\theta/(1+\cos^2\theta)\), którą punkt E wyprowadza
niezależnie jako kwantową wartość oczekiwaną skrętności
\(\langle h|\theta\rangle\) dla przejścia dipolowego \(\Delta m=\pm1\) —
nie przypadek, tylko standardowa zasada korespondencji (klasyczny
parametr Stokesa V fali EM = kwantowa wartość oczekiwana spinu fotonu w
granicy klasycznej), tutaj zweryfikowana wprost na tym konkretnym
systemie, nie tylko przywołana z zasady.

**Wnioski, jakie z tego wynikają:**

1. Wzorzec kątowy \((1+\cos^2\theta)\) (punkt C) i rozkład skrętności
   (punkt E) NIE są dwiema osobno zgadniętymi formułami, które akurat
   dobrze ze sobą współpracują — są dwiema twarzami tej samej klasycznej
   konfiguracji pola dalekiego, wzajemnie zweryfikowanymi do precyzji
   maszynowej. To niezależne potwierdzenie poprawności obu mechanizmów
   już zaimplementowanych, nie nowa fizyka.
2. Zarówno "gdzie" (środek masy jako wygodna konwencja, nie realny punkt
   źródłowy), jak i "brak wyzwalacza" (losowość zgodna z ortodoksyjną
   QED, nie uproszczenie) mają solidne uzasadnienie fizyczne — żadne z
   tych dwóch pytań nie ujawniło nowej niespójności modelu, w
   przeciwieństwie do audytu w punkcie I.
3. **Zastrzeżenie, uczciwie:** cała ta weryfikacja polaryzacji jest
   dokładna wyłącznie dla CZYSTO KOŁOWEJ orbity — to samo założenie, na
   którym stoi wzorzec \((1+\cos^2\theta)\) w punkcie C. Dla orbity
   mimośrodowej (a te teraz naprawdę występują, po poprawce spinu fotonu
   w punkcie E, sięgając \(e^2\approx0{,}9\)) prawdziwe pole dalekie
   niesie dodatkową strukturę harmoniczną — moc na wielokrotnościach
   częstości orbitalnej podstawowej, nie tylko na niej samej (klasyczny
   fakt: eliptyczna orbita Keplera ma widmo dipolowe o wielu harmonicznych,
   silniej obsadzonych w miarę wzrostu mimośrodu) — więc polaryzacja
   przestaje być tak czystą funkcją \(2\cos\theta/(1+\cos^2\theta)\), a
   model do niedawna tego nie łapał: każdy foton był losowany tak, jakby
   orbita była kołowa, niezależnie od jej faktycznego mimośrodu w chwili
   emisji. **Zmierzone ilościowo i naprawione w punkcie K poniżej** —
   skala okazała się o rząd wielkości większa, niż sugerowało to wstępne
   zdanie.

**K. Struktura harmoniczna promieniowania orbity mimośrodowej —
zmierzona, a nie tylko przewidziana, i wbudowana jako `CREM_HARMONIC`.**

*Skala problemu, zmierzona.* Rozłożono dipol Keplera na harmoniczne
częstości orbitalnej (numerycznie: równanie Keplera rozwiązywane wprost,
DFT pozycji \(r(t)e^{i\nu(t)}\) względem anomalii średniej, \(N=8192\)
próbek) i policzono, ile mocy promieniowania przypada na harmonikę
podstawową \(n=1\) (tę, na której opierał się dotąd
`photonEnergyReference=\hbar\omega_{orb}`), w funkcji mimośrodu:

| \(e\) | \(e^2\) | moc w \(n=1\) | harmonika szczytowa | mediana mocy przy \(n\le\) | 90% mocy przy \(n\le\) |
|---|---|---|---|---|---|
| \(0{,}10\) | \(0{,}010\) | \(96{,}1\%\) | \(1\) | \(1\) | \(1\) |
| \(0{,}50\) | \(0{,}250\) | \(33{,}1\%\) | \(1\) | \(2\) | \(4\) |
| \(0{,}70\) | \(0{,}490\) | \(8{,}4\%\) | \(3\) | \(5\) | \(11\) |
| \(0{,}80\) | \(0{,}640\) | \(2{,}7\%\) | \(5\) | \(9\) | \(21\) |
| **\(0{,}90\)** | **\(0{,}810\)** | **\(0{,}37\%\)** | **\(16\)** | **\(27\)** | **\(62\)** |
| **\(0{,}945\)** | **\(0{,}893\)** | **\(0{,}07\%\)** | **\(39\)** | **\(66\)** | **\(155\)** |
| \(0{,}970\) | \(0{,}941\) | \(0{,}014\%\) | \(97\) | \(166\) | \(388\) |

Przy mimośrodach, które model teraz rutynowo osiąga po poprawce spinu
(punkt E), poniżej \(0{,}1\%\) mocy jest w ogóle na częstości
podstawowej — rząd wielkości (albo dwa) większy problem niż sugerowała
wstępna wzmianka w punkcie J.

Sprawdzone: `larmorOrbitAveragedPower`/`dipoleEccentricityFactor` liczą
**poprawnie** całkowitą moc dla orbity mimośrodowej (pełny czynnik
\((1+e^2/2)/(1-e^2)^{5/2}\) z bezpośredniego całkowania
\(\langle r^{-4}\rangle\), jawnie odróżniony od czynnika Petersa dla
kwadrupola grawitacyjnego) — budżet energii i tempo hazardu były więc
zawsze poprawne. Błędny był wyłącznie podział tego poprawnego budżetu na
pojedyncze kwanty: przy \(\hbar\omega_{orb}\) dużo mniejszym niż
prawdziwa charakterystyczna częstotliwość, hazard \(=P/\hbar\omega_{ref}\)
wychodził sztucznie zawyżony — model strzelał znacznie za dużo, znacznie
za małych fotonów dla tego samego, poprawnego budżetu energii. Ponieważ
mechanizm spinu (punkt E) zdejmuje dokładnie \(\hbar\) momentu pędu **na
zdarzenie**, niezależnie od jego energii, nadmiar zdarzeń oznaczał
zawyżone tempo utraty momentu pędu właśnie w reżimie, w którym ta
poprawka miała największe znaczenie.

*Dwie osobne poprawki, nie jedna.* Naiwne przemnożenie
`photonEnergyReference` przez harmonikę szczytową \(n_{peak}(e)\) było
tylko prowizorką do pomiaru wrażliwości (patrz niżej). Rygorystyczne
domknięcie wymaga rozróżnienia dwóch rzeczy:

1. **Tempo zdarzeń.** Każda harmonika \(n\) to osobny, niezależny kanał
   Poissona o tempie \(P_n/(\hbar n\omega_{orb})\) — ta sama logika
   korespondencji Bohra, którą model już stosował dla \(n=1\), tylko
   uogólniona. Całkowity hazard \(=(P_{tot}/\hbar\omega_{orb})\cdot S(e)\),
   gdzie \(S(e)=\sum_n w_n/n\) (\(w_n\) — ułamek mocy w harmonice \(n\)).
   Zmierzone: \(S(e)\) spada z \(1\) przy \(e=0\) do \(0{,}061\) przy
   \(e=0{,}9\), \(0{,}022\) przy \(e=0{,}945\), \(0{,}0009\) przy
   \(e=0{,}99\) — znacznie mniej zdarzeń niż szacunek oparty wyłącznie na
   \(n=1\), bo dzielenie przez \(n\) tłumi wysokie harmoniki, które niosą
   większość MOCY, ale liczone po \(n=1\) dawałyby absurdalną liczbę
   zdarzeń.
2. **Która harmonika strzela.** Warunkowo, dany zaszły fakt zdarzenia,
   prawdopodobieństwo przynależności do harmoniki \(n\) jest
   proporcjonalne do \(w_n/n\), NIE do samego \(w_n\) — rzadkie, ogromne
   zdarzenie i tak zdarza się raz. Sprawdzone: nawet ważąc po liczbie
   zdarzeń, harmoniką o największej POJEDYNCZEJ gęstości wciąż jest
   \(n=1\) przy każdym testowanym \(e\) — ale rozkład ma na tyle ciężki
   ogon, że MEDIANA zdarzenia wciąż jest daleko od \(n=1\) (mediana
   \(n=11\) przy \(e=0{,}9\), \(n=27\) przy \(e=0{,}945\), \(n=64\) przy
   \(e=0{,}97\)): tylko \(\sim6\%\) zdarzeń to dosłownie \(n=1\) przy
   \(e=0{,}9\), spadające do \(\sim1{,}4\%\) przy \(e=0{,}97\). Czyli
   "wiele małych fotonów zamiast kilku dużych" było słuszną intuicją
   jakościową, tylko nieprecyzyjną co do dokładnego podziału między
   liczebność a wielkość.

*Co pozostaje bez zmian.* Wzorzec kątowy \((1+\cos^2\theta)\) i podział
skrętności \(P(h=\pm1|\theta)\) (punkty C, E) są używane BEZ ZMIAN na
każdej harmonice, nie uzasadnienie na wyrost: dla płaskiego ruchu
Keplera rozkład zespolony \(z=x+iy\) ma składową "prograde" (zgodną z
kierunkiem orbity) i "retrograde"; czysto obracający się dipol (\(n=1\),
\(e=0\)) jest w \(100\%\) prograde, co czyni \((1+\cos^2\theta)\) i
\(\Delta m=+1\) dokładnymi. Sprawdzone wprost na harmonikach, które
faktycznie niosą moc (szczytowa/medianowa/90-procentowa, \(e=0{,}9\) /
\(0{,}945\) / \(0{,}97\)): moc retrograde wynosi \(0{,}10\)–\(0{,}95\%\)
— te harmoniki są same w sobie w ponad \(99\%\) kołowo spolaryzowane w
TYM SAMYM sensie co orbita, więc użycie niezmienionego wzorca kątowego
przy częstości sprzężonej harmoniki jest przybliżeniem rzędu \(<1\%\),
nie nowym założeniem.

*Implementacja: `CREM_HARMONIC`.* \(S(e)\) jest stablicowane bezpośrednio
per mimośród (bez ekstrapolacji) i wpięte w `crem_collapse.hpp` pod
zmienną środowiskową `CREM_HARMONIC`.

*Pierwsza wersja funkcji kwantylowej dla "która harmonika strzela" była
zbudowana na osobnym założeniu — i to założenie, sprawdzone ponownie na
zadane polecenie, okazało się dziurawe.* Pierwotnie: jedna "uniwersalna"
funkcja kwantylowa zmiennej \(x=n/n_{peak}(e)\)
(\(n_{peak}(e)\approx0{,}504(1-e)^{-1{,}5}\)), skalibrowana na
\(e=0{,}75\)–\(0{,}98\) i zastosowana wszędzie, z twardą bramką
\(e\ge0{,}6\) poniżej progu kalibracji jako zabezpieczeniem. Sprawdzone
dokładniej, drugi raz: bramka na \(e\ge0{,}6\) była jednocześnie **za
ostrożna i niewystarczająca**. Za ostrożna — bezpośredni pomiar mediany i
90. percentyla przy \(e=0{,}2\)–\(0{,}5\) pokazał zgodność z tablicą do
\(0\)/\(+1\) harmoniki, czyli struktura powyżej \(n=1\) była tam już
realna i dobrze uchwytna, a bramka ją całkowicie wyłączała. Niewystarczająca
— nawet przy \(e=0{,}6\)–\(0{,}75\), tuż nad bramką, błąd w ogonie
(\(99\%\), \(99{,}9\%\)) wciąż wynosił \(+1\)/\(+2\) harmoniki. Był to
problem konkretnie ogona rozkładu przy niskim/średnim mimośrodzie, nie
mediany — założenie o jednym uniwersalnym kształcie po prostu nie
trzymało się poza zakresem, w którym zostało sprawdzone.

Naprawione przez usunięcie założenia o uniwersalności, nie przez
przesunięcie bramki: funkcja `eccentricOrbitHarmonicNumber` jest teraz
tablicą DWUWYMIAROWĄ (mimośród × kwantyl, \(18\times16\) wpisów),
każdy wpis policzony NIEZALEŻNIE dla własnego \(e\) (ta sama dekompozycja
DFT, \(N=8192\), do \(n_{max}=6000\)), interpolowaną dwuliniowo — bez
ekstrapolacji z innego zakresu mimośrodu, więc bramka w ogóle nie jest
już potrzebna: każdy wpis jest realną, zmierzoną liczbą, nie
przybliżeniem. Sprawdzone bezpośrednio na tej samej trajektorii, która
ujawniła pierwotny błąd (ziarno \(107\), pierwszy foton przy
\(e^2=0{,}0415\)): przegląd \(37\) zdarzeń fotonowych na \(22\)+ ziarnach
przy niskim mimośrodzie dał harmonikę \(1\) niemal zawsze, z dwoma
wyjątkami — harmonika \(3\) przy \(e^2\approx0{,}04\) i harmonika \(5\)
przy \(e^2\approx0{,}46\) — obie w pełni zgodne z nową tablicą, żadnej z
wcześniejszych absurdalnych wartości (np. harmonika \(9\) przy
\(e\approx0{,}2\), którą dawała pierwsza wersja).

*Zweryfikowane empirycznie, z tablicą dwuwymiarową.* Partia \(15\)
trajektorii (ziarno \(42\), budżet zegara \(60\) s): mediana
Kaplana-Meiera \(0{,}198748\) ns wobec \(0{,}197275\) ns bez poprawki
(\(+0{,}7\%\), w granicach szumu \(N=15\)) — potwierdza, że mediana
wyznacza całkowity budżet energii (poprawny niezależnie od
ziarnistości), nie liczba/wielkość fotonów. RMST \(0{,}304\) ns wobec
\(0{,}284\) ns (\(+7{,}0\%\)) — efekt skupiony w głębokim, mimośrodowym
ogonie, dokładnie tam, gdzie tabela pokazuje, że założenie o częstości
podstawowej najbardziej zawodzi; obie liczby zgodne co do rzędu wielkości
z pierwszą (wadliwą) wersją tablicy, jak należało oczekiwać, skoro obie
tablice zgadzają się dobrze tam, gdzie druga wersja faktycznie coś
zmieniła (ogon przy wysokim \(e\)) — różnią się tam, gdzie pierwsza
wersja była błędna (niski/średni \(e\)), a to wnosi tylko drugorzędną
poprawkę do zagregowanego wyniku. Z `CREM_HARMONIC=0`: wynik **bitowo
identyczny** ze stanem sprzed tej pracy (regresja sprawdzona wprost —
ta sama partia, te same liczby co do ostatniej cyfry, w obu wersjach
tablicy). Czysta kompilacja (zero ostrzeżeń), `positronium_validation`
33/33.

*Rozmiar efektu, przemierzony na większej próbie.* Partia powyżej —
jedno ziarno (42), piętnaście trajektorii — **odtwarza się dokładnie**:
uruchomiona ponownie daje RMST \(306{,}8\) ps wobec \(288{,}9\) ps, a
mediana zgadza się z zapisanym \(0{,}198748\) ns do **piątej cyfry
znaczącej**. Nie jest to więc pomiar błędny. Jest natomiast zbyt mały, żeby
nim rozmiar efektu wyznaczać: obie wartości RMST niosą \(\pm77\) i
\(\pm70\) ps przy różnicy \(18\) ps, czyli pojedyncza partia mierzy efekt
**czterokrotnie mniej precyzyjnie, niż wynosi sam efekt**.

Zmierzone sparowanie po **szesnastu ziarnach po dwadzieścia trajektorii**
(sparowanie jest tu istotne: rozrzut czasu kolapsu wynosi \(\sigma/\mu=0{,}57\)
i pochodzi głównie z losowych warunków początkowych, które przy tym samym
ziarnie są wspólne dla obu konfiguracji, więc w różnicy się kasują):

\[
\Delta_{\rm RMST}=+4{,}82\%,\qquad
\text{95\% CI }[+3{,}64\%;\,+6{,}01\%],\qquad
\text{znak dodatni w }16/16 .
\]

Rozrzut międzyziarnowy wynosi \(2{,}41\%\) (od \(+0{,}45\%\) do
\(+9{,}17\%\)), co stawia oryginalne \(+7{,}0\%\) na \(0{,}90\sigma\) —
zwykłe losowanie, nie rozbieżność. Mediana pozostaje bez zmian, zgodnie z
tym, co ta sekcja mówiła od początku: wyznacza ją budżet energii, którego ta
poprawka nie dotyka.

Osobno: komentarz w `crem_collapse.hpp` cytował tę wartość jako \(+7{,}7\%\),
podczas gdy dwie liczby podane wyżej dają \(0{,}304/0{,}284=+7{,}0\%\).
Poprawione.

*Jak często korekta w ogóle działa.* Zmierzony dobór harmoniki per foton,
\(N=83\) emisji na jedenastu ziarnach: **\(k=1\) w \(96{,}4\%\)
przypadków**, \(k=2\) w \(2{,}4\%\), \(k=3\) w \(1{,}2\%\). Powód jest w
mimośrodzie: w chwili emisji ma on medianę \(0{,}050\) (p90 \(0{,}18\), max
\(0{,}25\)), a mediany harmonik \(11\) i \(27\) z tabeli wyżej odpowiadają
\(e=0{,}9\) i \(0{,}945\), których emisje nie osiągają. Tabela harmonik jest
więc w produkcji niemal bezczynna — co czyniło pytanie, **z którego z dwóch
składników korekty** pochodzi mierzone \(+4{,}82\%\), pytaniem wymagającym
osobnego pomiaru, a nie wnioskowania.

*Rozkład efektu na składniki, zmierzony.* Korekta ma dwa niezależne
składniki i zostały rozdzielone: dobór harmoniki \(k\) z tablicy oraz
tłumienie hazardu \(S(e)\), które mnoży tempo emisji **przy każdym
fotonie, także gdy \(k=1\)**. Cztery konfiguracje (pełna / tylko tablica /
tylko \(S(e)\) / wyłączona), osiem ziaren po dwadzieścia trajektorii,
sparowane po ziarnach, RMST względem konfiguracji wyłączonej:

| składnik | średnia | 95% CI | \(t\) wobec zera |
|---|---|---|---|
| tylko tablica harmonik | \(+0{,}82\%\) | \([-1{,}51;\,+3{,}14]\) | \(0{,}69\) |
| **tylko tłumienie \(S(e)\)** | \(\mathbf{+4{,}44\%}\) | \([+2{,}74;\,+6{,}13]\) | \(\mathbf{5{,}13}\) |
| pełna korekta | \(+5{,}33\%\) | \([+2{,}31;\,+8{,}34]\) | \(3{,}46\) |

Składniki sumują się addytywnie: \(0{,}82+4{,}44=5{,}25\%\) wobec
\(5{,}33\%\) zmierzonych dla pełnej korekty.

**Cały efekt pochodzi z \(S(e)\).** Wkład tablicy harmonik jest zgodny z
zerem — co jest spójne z tym, że \(k=1\) w \(96{,}4\%\) emisji, więc
tablica rzadko ma cokolwiek do zmienienia. Praktyczna konsekwencja: tablica
\(18\times16\) i dekompozycja DFT do \(n_{max}=6000\), które ją wypełniają,
nie wnoszą do zagregowanego wyniku nic mierzalnego przy mimośrodach, jakie
produkcja osiąga; robi to pojedynczy skalar \(S(e)\). Nie jest to argument
za usunięciem tablicy — pozostaje poprawna i zaczęłaby działać przy wyższych
mimośrodach — lecz za tym, żeby nie przypisywać jej efektu, którego nie
wywołuje.

Zastrzeżenie metodyczne, to samo co wyżej: przebiegi pomiarowe zapisują
wykresy do `distributions/`, więc muszą być uruchamiane w osobnym katalogu
roboczym. Partia \(15\)–\(20\) trajektorii podmieniająca figury robione
przy \(N=1000\) jest cichą regresją, którą trzeba było raz cofać
(`c06934b`).

#### `--emission deterministic`: ten sam kwant, bez szumu śrutowego

Bankowany hazard to \(\int (P/E_\gamma)\,dt\), czyli **wypromieniowana
energia wyrażona w kwantach**. Próg emisji rozstrzyga więc tylko o jednym:

- próg \(\mathrm{Exp}(1)\) — emisja spontaniczna, z szumem śrutowym (domyślna),
- próg **stały równy 1** — foton pada dokładnie wtedy, gdy orbita straciła
  ciągle jeden kwant.

Drugi wariant jest tym, co potocznie znaczy „emituj przy przekroczeniu
poziomu": orbita nadal schodzi w sposób ciągły, a samo przekroczenie jest
wyzwalaczem. Różnica w kodzie to jedna stała, nie osobny model.

*Zmierzone,* osiem ziaren po dwadzieścia trajektorii:

| | Poisson | deterministyczny |
|---|---|---|
| \(\sigma/\mu\) | \(0{,}894\) | \(\mathbf{0{,}343}\) |
| rozrzut \(\sigma/\mu\) między ziarnami | \(0{,}778\) | \(0{,}080\) |
| średnia czasu kolapsu | \(161{,}9\) ps | \(167{,}3\) ps |
| rozrzut średnich między ziarnami | \(142{,}2\) ps | \(39{,}7\) ps |
| ziarna z cenzurowaniem | \(2/8\) | \(0/8\) |

Przy założeniu niezależności obu źródeł rozrzutu
(\(\sigma^2=\sigma_{\rm śrut}^2+\sigma_{\rm w.p.}^2\)) daje to rozkład
wariancji: **szum śrutowy odpowiada za \(85\%\)**, a rozrzut warunków
początkowych za \(15\%\). Reszta \(0{,}343\) jest przy tym niemal stała od
ziarna do ziarna (\(0{,}303\)–\(0{,}384\)), podczas gdy pełne \(\sigma/\mu\)
skacze od \(0{,}595\) do \(1{,}373\) — dokładnie tak, jak powinno, jeśli
usunięto składnik losowy, a został deterministyczny rozrzut przygotowania.

**Uściślenie wobec oczekiwania, z jakim ten tryb powstawał.** Spodziewano
się, że zawężenie sprowadzi rozkład „z powrotem do odpowiedzi kontinuum".
Nie sprowadza: średnia pozostaje przy \(167\) ps, czyli na skali trybu
stochastycznego, a nie przy \(36\)–\(40\) ps modeli ciągłych. Powód jest
strukturalny — średnią wyznacza budżet energii, czyli energie kwantów, które
w obu trybach emisji są identyczne; próg zmienia wyłącznie **rozkład czasów**.
Do odpowiedzi kontinuum prowadzi co innego: przywrócenie ciągłej siły
hamującej, a nie zdjęcie losowości z progu.

**Kierunek, któremu to służy: deterministyczne wyznaczanie parametrów
kwantowych.** Energia fotonu jest już wyznaczona przez orbitę — zasadą
korespondencji, a tam gdzie istnieje drabina, odstępem poziomów. Tryb
deterministyczny domyka to samo dla **czasu** emisji, przez co w kwancie nie
zostaje nic losowanego. Oba tryby są modelami produkcyjnymi, nie sondami.

Odstępstwo od opisu konwencjonalnego trzeba przy tym wypowiedzieć wprost, a
nie ukryć: emisja spontaniczna jest zwykle modelowana jako proces Poissona i
z tego powodu `poisson` pozostaje wartością domyślną. To, co daje tryb
deterministyczny, jest natomiast zmierzone: \(85\%\) wariancji czasu kolapsu
okazuje się szumem śrutowym, a nie własnością pary.

**Na wyraźną prośbę: `CREM_HARMONIC` stał się nowym domyślnym modelem
produkcyjnym.** Ta sama zmienna środowiskowa zostaje, ale zmienia się jej
domyślna wartość: poprawka jest teraz WŁĄCZONA, chyba że jawnie ustawi
się `CREM_HARMONIC=0` — to samo "domyślnie włączone, jeden przełącznik
od dawnego zachowania" jak przy promocji samego `stochastic` na domyślny
model reakcji wyżej. Zachowana, nie usunięta, właśnie po to, żeby dawne
zachowanie zostało o jedną flagę od odtworzenia, do porównań
regresyjnych. `CREM_HARMONIC=0` nadal odtwarza dawne zachowanie
dokładnie (dwa fotony w tej samej trajektorii ziarno \(107\) zamiast
jednego, zgodnie z wcześniejszym pomiarem — niezależne od wersji tablicy
powyżej, bo `CREM_HARMONIC=0` w ogóle jej nie dotyka). Czysta kompilacja,
`positronium_validation` 33/33. Pełny
przebieg produkcyjny \(N=1000\) dla wszystkich pięciu eksperymentów pod
tym nowym domyślnym ustawieniem nie został jeszcze przeliczony — mediana
czasu kolapsu nie powinna się zauważalnie zmienić (wyznacza ją budżet
energii, którego ta poprawka nie dotyka), ale rozkład (zwłaszcza ogon,
RMST) tak, więc historyczne liczby rozkładu gdzie indziej w tym
dokumencie pozostają zmierzone pod poprzednim ustawieniem, dopóki ten
przebieg nie zostanie powtórzony.

**L. Para kontra orto-pozytonium: czy któreś zasługuje na inną
harmonikę — zbadane sondami, nie wyargumentowane z fotela.**

*Pytanie wyjściowe.* Skoro para-Ps i orto-Ps różnią się realnie
istniejącym zjawiskiem kwantowym (2γ kontra 3γ), czy ten sam podział
mógłby uzasadniać różne traktowanie harmonik promieniowania orbitalnego
w tym klasycznym modelu?

*Sonda 1 — identyczne warunki startowe.* Dla tego samego ziarna, para i
orto dostają dokładnie te same \(E\), \(L\) na starcie (sprawdzone
bezpośrednio, `checkpoint 0` identyczny co do wyświetlanej precyzji) —
różnią się więc wyłącznie dynamiką, nie próbkowaniem warunków
początkowych.

*Sonda 2 — namierzenie miejsca i skali rozbieżności.* Prześledzone
checkpoint po checkpoincie (ziarno \(42\)): para i orto są **bitowo
identyczne przez pierwsze \(29\) checkpointów**, potem pojawia się
prawdziwa, mikroskopijna różnica (\(E\) różni się o \(\sim6\cdot10^{-6}\)
względnie). Zweryfikowane w kodzie: `runMechanicalTrajectory` używa
pełnego silnika (`ClassicalTrajectoryEngine`), który rzeczywiście
uwzględnia sprzężenie dipol-dipol — jedyną fizyczną różnicę między para
(momenty równoległe) a orto (przeciwrównoległe) w tym modelu.

*Sonda 3 — czy sprzężenie dipol-dipol mogłoby dać strukturę \(3\omega\)
zamiast \(2\omega\).* Rozłożona na harmoniczne (DFT) energia
\(U\propto\boldsymbol\mu_1\!\cdot\!\boldsymbol\mu_2-3(\boldsymbol\mu_1\!\cdot\!\hat{\mathbf r})(\boldsymbol\mu_2\!\cdot\!\hat{\mathbf r})\)
dla wirującego \(\hat{\mathbf r}(t)\):

```
PARA (równoległe):        DC=-0,5000  h1=0  h2=0,7500  h3=0  h4=0
ORTO (przeciwrównoległe):  DC=+0,5000  h1=0  h2=0,7500  h3=0  h4=0
```

Obie konfiguracje dają **dokładnie ten sam** wzorzec — brak pierwszej i
trzeciej harmonicznej, ta sama amplituda drugiej, różny tylko znak
składowej stałej. To nie przypadek: sprzężenie dipol-dipol jest
dwuliniowe w \(\hat{\mathbf r}\) (iloczyn dwóch rzutów na ten sam wirujący
wektor), a taki iloczyn matematycznie daje tylko składową stałą i drugą
harmoniczną — nigdy trzecią, do której potrzebny byłby człon trójliniowy,
którego w tym oddziaływaniu nie ma. **Nie ma więc naturalnego podziału
"2 dla para, 3 dla orto"** — a podobieństwo liczb 2 i 3 do prawdziwej
liczby fotonów anihilacji (2γ/3γ) jest numerologicznym zbiegiem
okoliczności, nie związkiem fizycznym: reguła wyboru 2γ/3γ pochodzi z
zachowania parzystości ładunkowej w zupełnie innym, kwantowym procesie
(anihilacji), którego CREM w ogóle nie modeluje.

*Sonda 4 — czy "zgodny, podwójny" moment magnetyczny należy do orto (jak
podejrzewano), czy do para.* Sprawdzony warunek próbkowania w
`crem_trajectory.hpp`:

```cpp
} while ((sampledScenario == 2 && dot(firstDipole,secondDipole)/(mu1*mu2) < 0.5)    // PARA: wymusza cos>=0.5
      || (sampledScenario == 3 && dot(firstDipole,secondDipole)/(mu1*mu2) >= 0.5)); // ORTO: wymusza cos<0.5
```

To **para** (nie orto) jest wymuszana na kąt \(\le60°\) między momentami
— bliżej zgodnych, "podwójny" moment. Orto dostaje kąt \(>60°\), od
prostopadłych po całkowicie przeciwne. Sens fizyczny: przeciwne ładunki
odwracają naiwną korelację spin–moment, więc antyrównoległe spiny (para,
\(S=0\)) dają zgodne momenty magnetyczne, a równoległe spiny (orto,
\(S=1\)) — przeciwne.

*Sonda 5 — koherentne promieniowanie M1, i dlaczego mimo to nie wraca do
orbity.* Moc dipola magnetycznego jest liczona koherentnie, nie
niezależnie dla każdej cząstki:

```cpp
const Vec3 totalSecondDerivative = first.secondDerivative + second.secondDerivative;
result.power = coefficient * totalSecondDerivative.squaredNorm();
```

Dwa zgodne źródła (para) dają \(4P_1\) (interferencja konstruktywna),
dwa przeciwne (orto) — moc bliską zeru (destrukcyjna) — to realny,
zweryfikowany w walidatorze mechanizm. Ale jego energia jest świadomie
i architektonicznie odcięta od orbity: księgowana do
`dipoleConstraintEnergy`, osobnego rezerwuaru, nigdy nie wchodzi do
budżetu energii orbitalnej ani do doboru harmoniki.

**Ostatnie zdanie jest nieścisłe w obie strony i zostało sprawdzone w
kodzie.** Po pierwsze, `dipoleConstraintEnergy` **jest** składnikiem
`conservativeParticleEnergy`, a więc formalnie wchodzi do budżetu energii
orbitalnej; z orbitą rozłącza je nie architektura, tylko **odejmowanie
tła** — przebieg tła drenuje ten sam rezerwuar, więc po odjęciu zostaje
\(4{,}0\cdot10^{-5}\) energii M1. Wniosek się broni, mechanizm był nazwany
błędnie. Po drugie, w trybie skwantowanym M1 **wchodzi** do doboru fotonu:
`quantizedPower` sumuje E1, M1 i E2 przed podzieleniem przez
\(\hbar\omega\). Rozbiór obu tych ścieżek jest przy omówieniu luki
implementacyjnej wokół kanału M1.

Nieścisła jest też sama dychotomia „para \(4P_1\), orto zero": znoszenie
jest dokładne wyłącznie przy \(\cos=-1\), a uśredniony po orientacjach
udział M1 wynosi \(1{,}88\cdot10^{-3}\) dla para wobec
\(1{,}90\cdot10^{-3}\) dla orto — różnica zgodna z zerem.

*Sonda 6 — czy moment siły reakcji M1 mógłby wracać do orbity inną
drogą, przez ewolucję samych wektorów dipola.* Sprawdzone w kodzie:
`state.firstProperDipole += reaction.firstDipoleTorque * (gyro*dt)` —
moment siły reakcji (koherentny, wrażliwy na wyrównanie) **rzeczywiście**
obraca wektory momentu magnetycznego, które z kolei wchodzą do siły
dipol-dipol napędzającej orbitę. Ścieżka istnieje naprawdę — ale
policzona ilościowo, z prawdziwymi stałymi przy \(a_{Ps}\):

```
Omega_BMT (precesja)                = 5,51e+11 rad/s
Omega_reaction (od reakcji M1)      = 6,26e-33 rad/s   <-- BLAD, patrz nizej
stosunek Omega_reaction/Omega_BMT   = 1,14e-44
```

**Ta liczba jest błędna i została później zmierzona ponownie.** Wartość
\(6{,}26\cdot10^{-33}\) rad/s pochodzi z ewaluacji momentu reakcji na
historii, która nie była jeszcze wypełniona — trzecia pochodna momentu,
\(\dddot{\mathbf m}\), wychodziła wtedy z szablonu numerycznego bliska zeru.
Policzona na wypełnionej historii i sprawdzona na siatce szablonów (\(10\) do
\(640\) kroków wypełnienia, dwa rozstawy różniące się czterokrotnie — wynik
stabilny w granicach czynnika \(1{,}7\)) daje przy \(a_{Ps}\)

```
Omega_reaction (od reakcji M1)      = 8,5e-01 rad/s
stosunek Omega_reaction/Omega_BMT   = 2,1e-12
```

czyli \(12\) rzędów wielkości, nie \(44\). Co więcej stosunek ten **nie
jest stały**: rośnie do \(1{,}2\cdot10^{-3}\) przy promieniu terminalnym.
Wniosek „formalnie istnieje, fizycznie martwa" trzeba więc odrzucić w tej
postaci — kanał nie jest martwy, jest tylko **szybki dopiero tam, gdzie para
prawie nie przebywa** — a ściślej: bo obrót od tego momentu siły jest
**oscylacyjny**, nie sekularny. Ilościowo rozstrzyga to pomiar bezpośredni
\(|\Delta\hat{\boldsymbol\mu}|\) między przebiegiem z momentem siły i bez
niego, ograniczony na poziomie \(\Omega_{\rm reakcji}/\omega_{\rm orb}\);
tabela szybkości i pomiar są przy omówieniu założenia izotropii momentów
wyżej w tym dokumencie. (Scałkowany \(\int\Omega_{\rm reakcji}\,dt\),
podawany tu wcześniej jako \(5{,}2\cdot10^{-3}\) rad, jest złym estymatorem
tej wielkości i został wycofany.)

*Sonda 7 — czy \(6\cdot10^{-6}\) z sondy 2 rzeczywiście pochodzi z
dipol-dipol, i z której jego części.* Tu poprzednia odpowiedź (w tej
samej pracy) była nieprecyzyjna — poprawiona na miejscu: **czysta
oscylacja \(2\omega\) uśrednia się do zera po każdej orbicie, więc nie
może dawać narastającej rozbieżności.** Za rozbieżność odpowiada
składowa STAŁA (DC) sprzężenia — a ta ma przeciwny znak dla para i orto
(\(-0{,}5\) i \(+0{,}5\) w jednostkach znormalizowanych z sondy 3), czyli
działa jak stała, dodatkowa poprawka do potencjału: przyciągająca dla
para, odpychająca dla orto. Policzone ilościowo, z prawdziwymi stałymi:

```
U_coulomb (a_Ps)                      = 13,61 eV
U_dd (skala)                          = 4,53e-05 eV
przewidywana różnica względna (para-orto) = 3,33e-06
zmierzona rozbieżność (checkpoint 29)     = ~6e-06
```

Zgodność co do rzędu wielkości (czynnik \(\sim1{,}8\), w granicach
przybliżenia: \(g\approx2\) zamiast zmierzonego \(2{,}00232\), czynniki
geometryczne rzędu jedności z konkretnej orientacji dipoli, \(a_{Ps}\)
jako przybliżenie promienia oskulacyjnego w chwili checkpointu \(29\)
zamiast dokładnej wartości) — to już nie wniosek przez eliminację, tylko
bezpośrednia, ilościowa zgodność przewidywania z pomiarem.

*Sonda 8 — czy ta mikroskopijna rozbieżność kiedykolwiek zmienia
harmonikę lub wynik, w praktyce.* Partia \(N=30\) dopasowanych
trajektorii para/orto (to samo ziarno): **identyczny rozkład wyników**
(\(26/30\) dotarło do granicy, \(4\) ucięte, w obu przypadkach),
identyczna mediana i RMST Kaplana-Meiera, i **\(44/44\) zdarzeń
fotonowych dało harmonikę \(1\) w OBU przypadkach**, z rozkładami \(e^2\)
zgodnymi do \(5\)–\(6\) cyfr znaczących przy każdym pojedynczym zdarzeniu.

*Sonda 9 — czy składowa DC działa jako precesja apsyd, ściśle
sprawdzone.* Na zadane pytanie, czy sekularna część sprzężenia mogłaby
wzmacniać/osłabiać mimośród (a więc harmonikę) jednego kanału względem
drugiego: zbudowany probierz — orbita Keplera z dodatkowym małym
zaburzeniem centralnym \(\propto1/r^3\) (dokładnie taka postać ma
sekularna część sprzężenia dipol-dipol), scałkowana RK4 przez wiele
orbit. Równanie Bineta dla \(F(r)=-GM/r^2-C/r^3\),
\(u''(\phi)+u(1-C/L^2)=GM/L^2\), jest **dokładnie liniowe niezależnie od
mimośrodu** — rozwiązanie ścisłe \(u(\phi)=\frac{GM}{L^2\Omega^2}
[1+e\cos(\Omega(\phi-\phi_0))]\), \(\Omega=\sqrt{1-C/L^2}\), z precesją
peryhelium **dokładnie** \(2\pi(1/\Omega-1)\) na orbitę — nie
przybliżeniem, tylko ścisłym wynikiem dla tej postaci siły. Pierwszy
pomiar numeryczny (pośredni, przez kierunek wektora Laplace'a-Rungego-Lenza,
który dla zaburzonego problemu nie jest ściśle zachowany) dał zły wynik:
\(-0{,}00528\) rad/orbitę, złego znaku i o czynnik \(\sim2\) różny od
przewidywania. **Błąd znaleziony i naprawiony**: wadliwe rozwijanie kąta
(unwrapping) w kodzie pomiarowym. Naprawione przez bezpośrednie
wykrywanie przejść przez peryhelium (lokalne minima \(r(t)\)):

```
przewidywanie analityczne (ścisłe) = 0,01038258 rad/orbitę
zmierzone (400 przejść przez peryhelium, poprawiona metoda) = 0,01038079 rad/orbitę
zgodność: 99,98%
```

Skoro równanie Bineta jest ściśle liniowe, amplituda oscylacji
\(u(\phi)\) — mimośród — jest ścisłą stałą całkowania: orbita precesuje
ze stałym tempem, jej kształt **nie zmienia się wcale**, do żadnego rzędu
w \(C\). Siła centralna nie wywiera momentu siły (moment pędu zachowany
dokładnie) i pozostaje zachowawcza (energia całkowita, z doliczonym
potencjałem zaburzenia, zachowana dokładnie) — mimośród, będący funkcją
tylko tych dwóch wielkości, nie ma więc żadnej drogi do dryfu.

*Sonda 10 — czy rozbieżność ma charakter skokowy, sprawdzone na surowych
danych.* Rozbita rozbieżność \(E\) checkpoint po checkpoincie (pełna
precyzja, nie \(6\) cyfr wyświetlanych domyślnie): checkpointy \(0\)–\(28\)
identyczne co do wyświetlanej precyzji, w tym przez PIERWSZE zdarzenie
fotonowe (checkpoint \(24\)) — rozbieżność pojawia się dopiero przy
DRUGIM zdarzeniu fotonowym (checkpoint \(29\)), i to w konkretnej
wielkości pochodnej (`cmEnergyKick` różni się o \(11\%\), mimo że
`photonEnergy` jest identyczna) — co wskazywałoby na różny kierunek
próbkowania drugiego fotonu, nie na inną jego wielkość. Wyglądało to na
potwierdzenie: skokowy, nie gładki charakter, dokładnie przy zdarzeniu
progowym procesu Poissona.

*Sonda 11 — sprawdzone jeszcze dokładniej, na żądanie, i obraz okazał
się bardziej subtelny niż sonda 10 sugerowała.* Wyciągnięta
`orbitalRadiatedEnergy` (rzeczywisty wynik całki strumienia pola dla
KAŻDEJ pojedynczej zmierzonej orbity, nie tylko zaokrąglone \(E\)) z
każdego checkpointu:

```
checkpoint 0: para=2,14536e-24J  ortho=2,14550e-24J  różnica względna=6,53e-05
checkpointy 1-24: ta sama różnica względna, konsekwentnie ~6,1-6,5e-05, KAŻDY checkpoint
```

**To obala wniosek sondy 10** — checkpointy \(0\)–\(24\) NIE są bitowo
identyczne, różnią się od samego początku, konsekwentnie, w tym samym
kierunku. Po prostu różnica jest zbyt mała, by być widoczna na \(6\)
cyfrach wyświetlanego `E=`, którego wartość bezwzględna jest dużo
większa. Zsumowana skumulowana różnica (dokładnie tak, jak wchodzi do
`elements.specificEnergy` przez kredytowanie energii zmierzonej orbity,
patrz punkt I) do checkpointu \(28\): \(3{,}45\cdot10^{-27}\) J, czyli
\(\approx7570\) J/kg w jednostkach `specificEnergy` — względnie
\(\approx4\cdot10^{-9}\). **To \(4\)–\(5\) rzędów wielkości za mało**,
żeby wytłumaczyć zmierzoną rozbieżność \(5{,}6\cdot10^{-6}\) przy
checkpoincie \(29\).

*Sonda 12 — głębsza instrumentacja (pełna precyzja podwójna, stan
strumienia PRNG, wektory pośrednie) w poszukiwaniu dokładnego mechanizmu
— i mechanizm znaleziony, ale zupełnie inny niż sondy 9–11 sugerowały.*
Dodana tymczasowa diagnostyka (`CREM_DEBUG_PRECISE`): pełna precyzja
`setprecision(17)` dla \(E\), \(L\), stanu strumienia losowego, oraz
bezpośredni wydruk `angularMomentumDirection` i każdego losowania
(`cardanoQ`, `cosThetaFromAxis`, `photonAzimuth`, `photonDirection`)
tuż przed każdym zdarzeniem fotonowym. Wynik, dla drugiego fotonu
(checkpoint \(29\), ziarno \(42\)):

```
PARA:  Ldir = ( 0,2359,  0,1797,  0,9550)   [PRZED pierwszym fotonem!]
ORTO:  Ldir = (-0,3180, -0,4810,  0,8170)
```

Mimo że `cosThetaFromAxis` i `photonAzimuth` (te same losowania!) są
bit w bit identyczne, `photonDirection` wychodzi kompletnie różny — bo
opiera się na innej osi. **To nie subtelny efekt sekularny — to zupełnie
inny wektor, obecny już przed pierwszym fotonem.**

Źródło znalezione w komentarzu kodu, który sam siebie obalił:
`angularMomentumDirection` był inicjalizowany z
`seedRun.frames.front().noetherAngularMomentum` — **pełnej** wielkości
Noethera, na twierdzeniu "wkład spinu/dipola do niej jest już
udokumentowany gdzie indziej jako \(\sim10^{-5}\) wyrazu orbitalnego".
To twierdzenie nigdy nie zostało bezpośrednio sprawdzone. Sprawdzone
teraz, dokładniej, na osobne żądanie — i wynik jest mocniejszy, niż
wstępnie sądzono: to nie przybliżenie "tego samego rzędu", tylko
precyzyjnie namierzalna pomyłka, z policzoną liczbą i znalezionym
źródłem.

*Dokładna wielkość wkładu spinu — nie oszacowanie, tylko odczytana
tożsamość już istniejąca w kodzie.* Funkcja `noetherAngularMomentum`
liczy `intrinsic = firstDipole/firstGyromagneticRatio +
secondDipole/secondGyromagneticRatio` — dokładnie \(\mu/\gamma_{gyro}\)
na cząstkę. `positronium.cpp`, przy funkcji `gyromagneticRatio`, ma już
istniejący (sprzed tej sesji, z zupełnie innej poprawki — brakującego
czynnika \(g\) w czterech miejscach kodu) komentarz:

```cpp
// The model carries |mu| = (g/2) * magneton, so the correct ratio returns
// S = hbar/2 exactly, as a spin-1/2 particle must.
```

Czyli \(\mu/\gamma_{gyro}=S=\hbar/2\) **dokładnie**, nie w przybliżeniu —
niezależnie ustalone przy zupełnie innej naprawie, teraz tylko odczytane
na nowo. Zestawione z rzeczywistym (nie specyficznym) orbitalnym
momentem pędu tej konkretnej trajektorii
(\(L_{orb}=0{,}00025761469604537469\times\mu_{red}\approx
1{,}1735\cdot10^{-34}\,\text{kg·m}^2/\text{s}\approx1{,}11\,\hbar\)):

```
wkład spinu (2 cząstki, w pełni wyrównane) = 2×(hbar/2) = hbar = 1,0546e-34 kg·m^2/s
L_orbitalny (ta trajektoria)                            ≈ 1,1735e-34 kg·m^2/s
stosunek = hbar / L_orbitalny ≈ 0,90
```

**Nie "tego samego rzędu wielkości" ogólnikowo — wkład spinu sięga aż
\(90\%\) wielkości orbitalnego momentu pędu**, gdy momenty są w pełni
wyrównane (przypadek para).

*Dokładne źródło błędnej liczby "\(\sim10^{-5}\)" — namierzone, nie
domyślane.* `positronium.cpp:1841` ma już istniejące, **poprawne**
zdanie: `"alignment, whose coupling is ~1e-5 of the Coulomb potential"`
— opisujące zupełnie INNĄ wielkość: stosunek energii sprzężenia
dipol-dipol do energii potencjalnej Coulomba (niezależnie zweryfikowany
wcześniej w sondach 3/7/9/11 jako rząd \(10^{-6}\)–\(10^{-5}\)). Pisząc
komentarz przy `angularMomentumDirection` przywołano najwyraźniej tę
samą liczbę "\(\sim10^{-5}\)" z pamięci, zastosowaną do **innego
pytania** (wkład spinu do momentu pędu) bez osobnej weryfikacji — dwie
różne wielkości fizyczne (energia kontra moment pędu), przypadkowo
skojarzone przez tę samą cyfrę zamiast osobno policzone.

Para i orto losują GENUINE różne (nie tylko przeciwnie skierowane)
wektory \(\boldsymbol\mu_1,\boldsymbol\mu_2\) pod warunkiem odrzucania,
który je definiuje — więc `angularMomentumDirection` wychodzi
kompletnie różny, dla KAŻDEJ pojedynczej trajektorii, nie tylko w
porównaniu para/orto.

*Sonda 13 — naprawa i weryfikacja.* Poprawka: zamiast pełnej wielkości
Noethera, `angularMomentumDirection` liczony bezpośrednio z geometrii,
używając dwóch klatek, które `seedOptions` i tak już żąda
(`frameCount=2`, `observationTime=10^{-24}` s — dwie próbki położenia
niemal w tej samej chwili na prawdziwej trajektorii):

```cpp
const Vec3 firstRelativePosition = frames.front().first-frames.front().second;
const Vec3 secondRelativePosition = frames.back().first-frames.back().second;
Vec3 angularMomentumDirection = cross(firstRelativePosition, secondRelativePosition);
```

Uzasadnienie: dla \(r_1=r_0+v\,dt\) przy małym \(dt\),
\(r_0\times r_1=r_0\times r_0+dt\,(r_0\times v)=dt\,(r_0\times v)\) —
dokładnie proporcjonalne do prawdziwego kierunku orbitalnego momentu
pędu, bez żadnego zanieczyszczenia spinem, i bez zależności od
konkretnej konwencji próbkowania (w przeciwieństwie do zakodowania na
sztywno osi \(z\), co akurat działałoby dla sposobu, w jaki TEN plik
próbkuje stany związane, ale cicho zawiodłoby dla każdego innego
wywołania). Sprawdzone bezpośrednio (to samo ziarno): normalizuje się do
\((0,\,2\cdot10^{-13},\,1)\) — oś \(z\) odzyskana do szumu
zmiennoprzecinkowego, nie do różnicy rzędu jedności, jaką dawało
zanieczyszczenie spinem.

**Weryfikacja po naprawie — kompletna zgodność:**

```
Ldir przed fotonem #1:    PARA=(0,-1,2e-13,1)  ORTO=(-0,1,7e-13,1)  [oba ~oś z]
photonDirection foton #1: identyczny do 10+ cyfr znaczących w OBU kanałach
cmEnergyKick foton #1:    1,28403e-23J w OBU (wcześniej: identyczne już i tak)
Ldir przed fotonem #2:    identyczny do wielu cyfr w OBU kanałach
photonDirection foton #2: identyczny do wielu cyfr w OBU kanałach
cmEnergyKick foton #2:    5,46842e-22J w OBU (wcześniej: 5,40e-22 vs 5,99e-22!)
```

Partia \(N=30\) (to samo ziarno co wcześniej): mediana i RMST Kaplana-Meiera
bez zmian względem sprzed naprawy — zgodnie z oczekiwaniem, bo poprawka
dotyczy KIERUNKU pojedynczych fotonów (uśrednia się w statystyce
zbiorczej), nie budżetu energii. `positronium_validation` 33/33.

*Sonda 14 — charakter skokowy zbadany PONOWNIE, po naprawie, żeby
wykluczyć, że sondy 10–11 były same zafałszowane przez błąd sondy 12.*
Zasadne pytanie: skoro `angularMomentumDirection` był kontaminowany
spinem PRZEZ CAŁY CZAS (nie tylko w chwili losowania fotonu), to czy
wniosek sondy 11 („rozbieżność ciągła od checkpointu 0, ale o \(4\)–\(5\)
rzędów wielkości za mała") sam nie był zniekształcony przez ten sam błąd?
Powtórzone od zera na naprawionym kodzie (`CREM_DEBUG_PRECISE`, ziarna
\(42\) i \(7\)), rozbita `measured dE/E` (rzeczywiste tempo strat
zmierzone na pojedynczej orbicie, wielkość napędzająca hazard) checkpoint
po checkpoincie, jako stosunek orto/para:

```
ziarno 42, checkpointy 0-17 (przed fotonem #1): stosunek = 1,000065  (stała wartość na wszystkich 18)
ziarno 42, checkpointy 18-24 (przed fotonem #1): stosunek = 1,000062  (gładki, monotoniczny dryf w dół)
ziarno 42, checkpointy 25-28 (po fotonie #1):    stosunek = 1,000009  (nowe, znów stałe plateau)
ziarno 7,  checkpointy 0-3:                      stosunek ≈ 1,000068  (stałe, ta sama skala)
```

Żadnego skoku — gładki, ciągły, powoli dryfujący ułamkowy offset od samego
początku, dokładnie taki, jakiego oczekuje się po stałym, sekularnym
zaburzeniu \(1/r^3\) (sonda 9), nie po zdarzeniu progowym. Same akty
fotonowe są teraz, po naprawie, niemal identyczne między kanałami:
`cardanoQ`, `cosThetaFromAxis`, `photonAzimuth`, `harmonicNumber`
identyczne co do bitu, `photonDirection` zgodny do \(10\)+ cyfr, a **czas
kolapsu Kaplana-Meiera identyczny w obu kanałach dla tego samego ziarna**
(\(757{,}471\) ps dla ziarna \(42\), \(81{,}7338\) ps dla ziarna \(7\)) —
czego przed naprawą nie było. Wniosek: hipoteza „skoku dokładnie przy
zdarzeniu fotonowym" z sondy 10 była w całości artefaktem błędu z sondy
12 — sam błąd był obecny od \(t=0\) (zły, ale STAŁY kierunek osi),
lecz pozostawał niewidoczny, dopóki nie został wciągnięty do próbkowanej
wielkości (kąt Cardano, azymut) w dyskretnym akcie losowania fotonu; stąd
złudzenie skoku. Prawdziwe sprzężenie dipol-dipol takiego mechanizmu nie
ma — jest, i zawsze było, czysto ciągłe.

*Sonda 15 — sam akt emisji rozbity na część losową i część fizyczną.*
Na żądanie: gdzie dokładnie w akcie emisji (nie tylko w tempie strat
między aktami) jest jakakolwiek różnica para/orto? Rozbite osobno:

Część losowa/dyskretna aktu (który harmonik, pod jakim kątem) —
**identyczna co do bitu** w obu zdarzeniach fotonowych: `threshold`,
`harmonicNumber` (\(1=1\)), `cardanoQ`, `cosThetaFromAxis`,
`photonAzimuth`. Wynika to wprost z architektury: przy tym samym ziarnie
stan strumienia PRNG (`stream=`) jest identyczny w obu kanałach w chwili
każdego losowania, a rozkład kątowy \((1+\cos^2\theta)\), z którego
losowanie korzysta, nie zależy od parametrów orbity ani spinu.

Część ciągła/fizyczna (\(E\), hazard, energia fotonu, `vcm`) różni się —
ale zawsze o ten sam mały, gładko narastający rząd wielkości, policzone
dokładnie (nie na oko: pierwszy ręczny szacunek dla `vcm` dawał
błędnie \(\sim3{,}5\cdot10^{-6}\) — pomyłka o czynnik \(1000\) w
przesunięciu przecinka, wykryta przy przeliczeniu skryptem):

```
                          zdarzenie 1 (chk 24)   zdarzenie 2 (chk 28)
E (różnica względna)          4,16e-09               3,69e-09
hazard (różnica względna)     4,28e-10               7,44e-09
photonEnergy (różn. wzgl.)    3,50e-09               7,91e-09
vcm, każda składowa           -  (=0 przed 1. fotonem)  3,50e-09
Ldir (różnica względna)       ~1e-13 (szum)          ~1e-12 do 1e-16 (szum)
```

`vcm` przy drugim zdarzeniu różni się o dokładnie tyle, ile w tej samej
chwili wynosi różnica \(E\) — bo kopnięcie od fotonu #1 dziedziczy jego
energię, a ta dziedziczy \(E\) orbity w chwili emisji; żadnego
niezależnego, dodatkowego mechanizmu przy `vcm` nie ma. `Ldir` zostaje na
poziomie czystego szumu numerycznego w obu zdarzeniach, potwierdzając, że
naprawa rzeczywiście usunęła zależność kierunku od spinu. Wniosek: sam
**akt** emisji (które zdarzenie, w jakim kierunku względem osi) nie niesie
żadnej niezależnej różnicy para/orto — cała różnica, jaka jest, to
odziedziczona, jedna i ta sama, ciągle rosnąca wartość \(E\), obecna od
\(t=0\) (sonda 11, teraz potwierdzona sondą 14 jako niezafałszowana przez
błąd sondy 12), bez żadnego skoku przypisanego samemu momentowi losowania.

*Sonda 16 — czy ta różnica rośnie w miarę zbliżania się do granicy
kolapsu, sprawdzone na dłuższej trajektorii.* Ziarno \(23\): \(68\)
checkpointów, \(2\) zdarzenia fotonowe, \(|E|\) rosnące (a więc orbita
kurcząca się) niemal \(29\)-krotnie zanim trajektoria dotrze do bariery
Comptona. Skumulowana różnica względna \(E\) między kanałami, cała
trajektoria:

```
checkpoint  0  (start):                    różnica = 0
checkpoint 42  (tuż przed fotonem #1):      różnica = 6,69e-08   |E| ×1,0002
  --- FOTON #1 ---
checkpoint 43  (tuż po fotonie #1):         różnica = 6,06e-08   |E| ×4,05
checkpoint 66  (tuż przed fotonem #2):      różnica = 7,46e-08   |E| ×4,05 (bez zmian)
  --- FOTON #2 ---
checkpoint 67  (bariera Comptona, koniec):  różnica = 1,49e-07   |E| ×28,9
```

Rośnie, ale nie jednorodnie — dwa osobne mechanizmy, rozdzielone
pomiarem: **(a)** wewnątrz jednej fazy (między fotonami) narasta liniowo z
liczbą zmierzonych orbit niemal niezależnie od kurczenia się orbity
(orbita w obrębie fazy ledwo się zmienia — to kumulacja czasu ekspozycji,
nie efekt zbliżania się do bariery); **(b)** sam wewnętrzny offset tempa
strat (`measured dE/E`, oczyszczony z efektu kumulacji czasu) rośnie
naprawdę, ale skromnie: \(|1-\text{stosunek orto/para}|\) idzie z
\(3{,}70\cdot10^{-4}\) (faza \(1\), duża orbita) do \(4{,}67\cdot10^{-4}\)
(faza \(2\), orbita \(4{,}05\times\) mniejsza) — wzrost o \(26\%\) przy
\(4\times\) skurczeniu, dużo mniej, niż sugerowałoby naiwne skalowanie
\(U_{dd}/U_{Coulomb}\sim1/r^2\) (które przewidywałoby \(\sim16\times\)).
**(c)** Same skoki przy fotonach dominują nad powolnym narastaniem
sekularnym: przy fotonie #2 różnica podwoiła się niemal dokładnie
(\(7{,}46\to14{,}9\cdot10^{-8}\)) — nowo wyemitowany foton składa
dotychczasowy skumulowany offset w nową, znacznie większą bazę \(|E|\);
to artefakt dyskretnej architektury kwantyzacji Poissona, nie ciągłego
wzmocnienia fizyki sprzężenia. Wniosek: różnica rośnie w stronę kolapsu
(łącznie o \(\sim2\) rzędy wielkości na całej trajektorii tutaj), ale
wzrost jest zdominowany przez strukturę zdarzeń fotonowych i kumulację
czasu, nie przez gwałtowne wzmocnienie sprzężenia dipol-dipol —
"czysto fizyczny" składnik rośnie realnie, lecz znacznie wolniej, niż
naiwna skala \(1/r^2\) sugerowałaby. Nawet tuż przy barierze Comptona
różnica pozostaje rzędu \(10^{-7}\)–\(10^{-9}\) zależnie od trajektorii —
wciąż o wiele rzędów wielkości za mało, by wpłynąć na wybór harmoniki czy
jakikolwiek obserwowalny wynik (sonda 8).

*Sonda 17 — czy dyskretność emisji (proces Poissona) i narastanie
sprzężenia mogą się ze sobą sprzęgać, dając PRAWDZIWY skok — zbadane
wprost, nie tylko oszacowane.* Zasadne pytanie: hazard, który decyduje,
W KTÓRYM checkpoincie strzeli foton, sam jest zbudowany z \(E\), więc
niesie tę samą maleńką różnicę sprzężenia. Gdyby skumulowany hazard
kiedykolwiek wylądował bliżej progu niż wynosi to zaburzenie, próg
zostałby przekroczony w INNYM checkpoincie w jednym kanale niż w
drugim — prawdziwy, duży, nieciągły skok, jakościowo różny od
wszystkiego dotąd zaobserwowanego, i genetycznie związany ze sprzężeniem
(bo to ono decydowałoby o przechyleniu), nie z samą architekturą.

*Metoda.* Zamiast szacować, przeszukane **\(959\) rzeczywistych zdarzeń
fotonowych** (setki ziaren, jeden kanał): dla każdego margines do progu
(`próg − cumHazard tuż przed odpaleniem`), znormalizowany przez wielkość
przyrostu hazardu na checkpoint (`skipHazard`). Jeśli miejsce
przekroczenia progu wewnątrz przyrostu jest z grubsza losowe, minimum z
\(N\) takich prób powinno skalować się jak \(1/N\) — sprawdzalna
predykcja, nie założenie:

```
N = 959 zdarzeń
zmierzone minimum (margines/skipHazard) = 5,20×10⁻⁴
przewidywanie z rozkładu jednostajnego (1/N) = 1,04×10⁻³
```

Zgodność co do rzędu wielkości (czynnik \(\sim2\)) potwierdza model.
Najciaśniejszy znaleziony przypadek (ziarno \(520\), zdarzenie fotonowe
\(\#2\)): margines \(2{,}75\times10^{-5}\) absolutnie, \(0{,}041\%\)
progu — **przetestowany wprost**, pełne porównanie para/ortho z
`CREM_DEBUG_PRECISE`:

```
para:  foton #2 odpala checkpoint 10->11
ortho: foton #2 odpala checkpoint 10->11   <- TEN SAM checkpoint
photonEnergy, cmEnergyKick: identyczne co do wyświetlanej precyzji w obu kanałach
```

Żadnego przechylenia — ale zmierzone bezpośrednio zaburzenie hazardu w
tym miejscu (\(\sim3\)–\(7\times10^{-9}\) absolutnie) jest wciąż
\(\sim10^4\) razy mniejsze od marginesu, nie \(\sim10^6\), jak sugerował
pierwszy, luźniejszy szacunek z sondy poprzedniej.

*Sprawdzone jeszcze dokładniej, na żądanie — i znaleziony realny błąd
jednostek w pierwszym szacunku.* Pierwsza wersja dzieliła zaburzenie
WZGLĘDNE (unormowane do lokalnej skali hazardu, która waha się od
\(0{,}067\) do \(4{,}16\) między zdarzeniami) przez `skipHazard` tak, jakby
hazard wszędzie miał skalę \(\approx1\) — niespójne jednostki. Poprawka:
liczyć \(P(\text{przechylenie})=|\text{hazard}_{orto}-\text{hazard}_{para}|
/\text{skipHazard}\) bezpośrednio z wielkości bezwzględnych, zdarzenie po
zdarzeniu, na wszystkich \(9\) zmierzonych dotąd zdarzeniach (\(6\)
ziaren):

```
zdarzenie          hazard    próg      skipHazard   zaburzenie(abs)   P(przechylenie)
seed10 ev1        0,687625  0,641708   0,114606      3,04e-10          2,65e-09
seed12 ev1        0,702371  0,653267   0,100341      5,59e-10          5,57e-09
seed42 ev1        2,832412  2,781617   0,113302      1,21e-09          1,07e-08
seed23 ev1        4,252913  4,156797   0,098919      2,77e-09          2,80e-08
seed520 ev1       1,011772  0,998010   0,101180      3,45e-09          3,41e-08
seed42 ev2        0,326824  0,307550   0,069007      2,43e-09          3,52e-08
seed15 ev1        3,061136  2,961908   0,102045      1,40e-08          1,37e-07
seed520 ev2       0,119592  0,066704   0,052916      7,87e-09          1,49e-07
seed23 ev2        1,349880  1,322060   0,052241      1,02e-07          1,95e-06

min = 2,65e-09    mediana = 3,41e-08    max = 1,95e-06    rozrzut = ×738
```

**Odkrycie, którego pierwszy szacunek nie złapał: to nie jest jedna
stała.** Pogrupowane wg tego, czy to pierwszy czy drugi foton danej
trajektorii — "zdarzenie #1" (\(6\) pomiarów): \(2{,}65\times10^{-9}\)
do \(1{,}37\times10^{-7}\), mediana \(\approx1{,}9\times10^{-8}\);
"zdarzenie #2" (\(3\) pomiary): \(3{,}52\times10^{-8}\) do
\(1{,}95\times10^{-6}\), systematycznie wyżej. Drugie zdarzenia fotonowe
(dłużej wystawione na skumulowane sprzężenie) mają wyraźnie wyższe
\(P(\text{przechylenie})\) niż pierwsze — dokładnie zgodnie z sondą 16
(różnica rośnie w stronę kolapsu). Ryzyko bifurkacji więc **rośnie z
każdym kolejnym fotonem tej samej trajektorii**, nie jest płaskie.

*Wniosek, dwukrotnie poprawiony względem pierwszego szacunku.* Poprawny
zakres: \(P(\text{przechylenie})\sim3\times10^{-9}\) (wczesne fotony) do
\(\sim2\times10^{-6}\) (późne fotony, bliżej kolapsu), mediana
\(\sim3\times10^{-8}\) — dolna granica niżej, a rozrzut (\(\times738\),
nie \(\times154\)) szerszy, niż sugerował poprzedni, luźniej policzony
przedział \(10^{-8}\)–\(10^{-6}\). Walidacja order-statystyką (zgodność
minimum marginesu z \(1/N\), czynnik \(\sim2\)) pozostaje słuszna co do
rzędu wielkości, ale to sprawdzian pojedynczej liczby (minimum), nie
pełny test kształtu rozkładu — traktowane jako przybliżona zgodność, nie
dowód ścisły. Mechanizm jest **potwierdzony jako realny i ilościowo
spójny** (bezpośredni test najciaśniejszego z \(959\) przypadków nie dał
przechylenia, dokładnie zgodnie z przewidywaniem \(\sim10^4\)-krotnego
marginesu w tamtym konkretnym przypadku) — ale nieobserwowany w żadnym z
\(959\)+ sprawdzonych zdarzeń ani w partii \(N=30\) (sonda 8). To osobny
od sondy 14 wniosek: sama fizyka sprzężenia jest gładka i ciągła
(sonda 14), ale **architektura kwantyzacji Poissona, w którą jest wpięta,
wprowadza osobny, rzadki, ale niezerowy, wprost zmierzony i systematycznie
rosnący w stronę kolapsu mechanizm bifurkacji** — próg, nie sprzężenie,
jest źródłem potencjalnego skoku, gdyby kiedykolwiek doszło do
wystarczająco ciasnego trafienia.

*Wniosek końcowy — mechanizm w pełni zidentyfikowany i naprawiony, nie
tylko domknięty w statystyce.* Poprzednia hipoteza (precesja apsyd,
sondy 9–11) była matematycznie poprawna jako ZJAWISKO OGÓLNE, ale nie
była właściwym mechanizmem OBSERWOWANEJ rozbieżności — prawdziwą
przyczyną była kontaminacja osi emisji fotonu przez losowo próbkowany
spin, niezwiązana z sekularnym sprzężeniem dipol-dipol wcale. Naprawione,
nie tylko udokumentowane: `angularMomentumDirection` liczony teraz z
czystej geometrii orbitalnej, poprawnie dla KAŻDEJ trajektorii
stochastycznego modelu, nie tylko w porównaniu para/orto. Sprzężenie
dipol-dipol (sondy 3, 9, 11, potwierdzone ponownie sondami 14–16) pozostaje
realne, małe (\(\sim10^{-6}\) energii Coulomba, \(\sim10^{-9}\)–\(10^{-7}\)
w skumulowanej różnicy \(E\), rosnącej w stronę kolapsu głównie przez
strukturę zdarzeń fotonowych, nie przez gwałtowne wzmocnienie samego
sprzężenia — sonda 16), o strukturze DC+\(2\omega\), i — teraz zbadane
wprost na naprawionym kodzie — **samo w sobie bez żadnego charakteru
skokowego, całkowicie gładkie i ciągłe od \(t=0\)** (sonda 14). Jedyny
znaleziony kanał, przez który mogłoby to się zmienić w skok, jest
architektoniczny, nie fizyczny: rezonans ze skończonym progiem procesu
Poissona (sonda 17) — mechanizm realny i ilościowo zmierzony
(\(P\sim3\times10^{-9}\)–\(2\times10^{-6}\) na zdarzenie fotonowe, rosnąco
z każdym kolejnym fotonem tej samej trajektorii), ale w \(959\)
sprawdzonych zdarzeniach i partii \(N=30\) (sonda 8) nigdy nieaktywowany.
Wciąż bez podstaw do systematycznie różnych harmonik między kanałami,
teraz bez domieszki znacznie większego, niepowiązanego artefaktu.
Siedemnaście niezależnych sond, w tym cztery własne pomyłki znalezione i
skorygowane w trakcie (sonda 9: błąd rozwijania kąta; sonda 10→11:
przedwczesny wniosek o skokowym charakterze; `noetherAngularMomentum` z
sondy 12 — błąd sprzed wielu sesji, wykryty dopiero teraz; sonda 15:
pomyłka o czynnik \(1000\) w ręcznym szacunku `vcm`, złapana natychmiast
przeliczeniem) — dokładnie taki proces, jaki ta cała sekcja miała
demonstrować.

**M. Walidacja produkcyjna, pełna skala: wszystkie pięć eksperymentów
statystycznych, \(N=1000\), e⁺e⁻, bez zewnętrznego pola magnetycznego.**
Domknięcie punktów A–L: nie na przykładach czy próbkach po kilkadziesiąt
trajektorii, tylko na pełnej, produkcyjnej partii każdego eksperymentu
(`--mode statistical`, ziarno \(42\) — to samo, którego użył autor kodu
do własnych pomiarów czasu wykonania w komentarzu), po próbie 4 i
naprawie `hFraction`. Łączny czas: \(1\) h \(25\) min (`--phenomenon 1`:
\(35\) min \(47\) s; `2`: \(36\) min \(17\) s; `3`: \(5\) min \(49\) s;
`4`: \(2\) min \(47\) s; `5`: \(4\) min \(17\) s — dłużej niż
udokumentowane w kodzie orientacyjne \(\sim51\) min, prawdopodobnie inna
liczba rdzeni na tej maszynie niż przy oryginalnym pomiarze).

*Eksperymenty 1–2 (CREM, para/orto-Ps) — najsilniejsza dotychczasowa
walidacja tej sekcji.* **Zero awarii numerycznych na \(2000\) łącznie
zmierzonych trajektoriach** (\(1000\) para + \(1000\) orto), obie z
identyczną liczbą ukończeń (\(926/1000\)) i ucięć budżetem czasowym
(\(74/1000\)):
\[
\begin{array}{lcc}
 & \text{para-Ps} & \text{orto-Ps}\\
\text{mediana Kaplana-Meiera} & 119{,}843\text{ ps} & 120{,}843\text{ ps}\\
\text{RMST} & 187{,}428\pm7{,}13\text{ ps} & 187{,}449\pm7{,}14\text{ ps}
\end{array}
\]
Mediany różnią się o \(\approx0{,}83\%\), RMST praktycznie identyczne —
zgodne co do rzędu wielkości z tym, co cała sekcja L już ustaliła
(para/orto nierozróżnialne w granicach tego mechanizmu), teraz
potwierdzone na pełnej próbie produkcyjnej, nie tylko na garści
przykładów. Panele lab-frame (punkt E2): para-Ps dała \(2735\) fotonów
na \(1000\) trajektorii (\(\langle\text{fotony}\rangle=2{,}73\); rozkład
\(2\) z \(1\) fotonem, \(266\) z \(2\), \(727\) z \(3\), \(5\) z \(4\) —
kaskady \(2\)–\(3\)-fotonowe to teraz norma, \(99{,}3\%\) trajektorii),
\(\langle E_{lab}\rangle=738{,}1\) eV, rozkład silnie prawoskośny (ogon
do \(\sim90\) keV) — zgodne z ustaloną w punkcie E4 fizyką kaskad
(kolejne fotony niosą coraz więcej energii).

*Eksperymenty 3–5 (wiązka, poza zakresem tej sekcji, ale ta sama partia
produkcyjna) — dla kompletności.* Eksperyment \(5\) (klasyfikacja
zderzeń, \(N=1000\)): \(216\) kolizji bezpośrednich (\(21{,}6\%\)),
\(691\) rozproszeń (\(69{,}1\%\)), \(93\) stany związane — \(28\)
para-Ps, \(65\) orto-Ps (stosunek \(28{:}65\), oczekiwanie izotropowe
\(1{:}3\) dałoby \(23{:}70\) na tej samej próbie — zgodne w granicach
szumu statystycznego \(N=93\)). Eksperymenty \(3\)–\(4\) (wiązka,
kanał krótkozasięgowy i rozpraszanie sprężyste): własne, niezależne od
CREM diagnostyki rezerwuaru energii (`|E_bound|/|E_rel|`) zakończone bez
awarii, \(0/1000\) przypadków przekraczających energię orbity w
eksperymencie \(4\).

**N. Konwersja czasu kolapsu i mocy promieniowania (CREM, eksp. \(1\)/\(2\))
do układu laboratoryjnego.** Audyt wykazał, że wykresy/statystyki
oznaczone literą `b` (dane wyjściowe) dla CREM — czas kolapsu, moc
promieniowania — były liczone w \(S'\), chwilowym układzie spoczynkowym
pary, a nie w układzie laboratoryjnym, w którym faktycznie następuje ich
pomiar. Dla wewnętrznej fizyki (bilanse energii/pędu/momentu pędu, punkt
L) \(S'\) jest poprawnym i wystarczającym układem — ale sam czas trwania
trajektorii i uśredniona moc promieniowania to wielkości makroskopowe,
którym lab-frame się należy, skoro `centreOfMassVelocity` realnie rośnie z
każdym odrzutem fotonu.

*Implementacja.* `simulatedTimeTotal` (czas własny \(S'\)) ma teraz
odpowiednik `labFrameTimeTotal`, całkowany z tym samym przyrostem na
checkpoint, ale ważony przez \(\gamma(\beta)\) osobno dla każdego
odcinka stałej \(\beta\) między kolejnymi odrzutami fotonu w tym samym
checkpoincie — pozycja `sAtPhoton` każdego fotonu (już liczona dla
kaskady, punkt E4) dzieli czas trwania checkpointu na segmenty w
zamkniętej formie \(T(s)=T_0\,n\,(s/s_{max})(1-s/2)\), a każdy segment
mnożony jest przez \(\gamma\) prędkości CM, jaka obowiązywała PRZED
kolejnym odrzutem. `CremCollapseEstimate` zyskał `lifetimeSecondsLab`,
`meanRadiatedPowerWattsLab` i `calibrationSecondsLab` (ten ostatni — dla
prawoucinanych/nieudanych przebiegów, żeby estymator Kaplana-Meiera nie
mieszał układów między obserwacjami dopełnionymi a cenzurowanymi).
`decayTimes`/`survivalSample`/`calibrationPowers` w `positronium.cpp`
(zasilające `crem_collapse_time`, `collapse_time_distribution`,
`diagnostic_calibration_power`) przełączone na te pola.

*Co świadomie zostało w \(S'\), mimo litery `b`:* `collapse_time_vs_theory`
(panel \(3\)) porównuje zmierzony czas z `classicalInspiralSeconds` —
zamkniętą, beznadziejnie sparametryzowaną (zero parametrów swobodnych)
formułą liczoną wprost z elementów oskulacyjnych, bez pojęcia prędkości
odrzutu do doboostowania; przeliczenie tylko jednej strony porównania
wprowadziłoby sztuczne, niefizyczne przesunięcie. `radiated_power_vs_larmor`
(panel \(4\), `larmorPowerRatio`) to stosunek dwóch wielkości chwilowych z
TEGO SAMEGO checkpointu (zmierzona strata na jednej orbicie vs Larmor dla
tej samej elipsy oskulacyjnej) — nie ma tu żadnego "czasu trwania" do
dylatacji, więc lab-frame nie ma zastosowania.

*Zmierzony efekt.* Zweryfikowane na kilkunastu ziarnach (`CREM_DEBUG=1`,
precyzja `setprecision(12)`): \(t_{lab}/t_{S'}-1\) mieści się w
\(10^{-11}\)–\(10^{-12}\) nawet dla trajektorii z pełną \(3\)-fotonową
kaskadą i szczytową prędkością odrzutu \(\beta\sim0{,}35\%\) (np. ziarno
\(5\): \(t_{S'}=287{,}308925168\) ps, \(t_{lab}=287{,}308925171\) ps).
Przyczyna: podwyższona \(\beta\) utrzymuje się tylko przez znikomy
ułamek całkowitego czasu trajektorii (ostatnie checkpointy tuż przed
kolapsem trwają femtosekundy, podczas gdy cała trajektoria to setki
pikosekund) — ważona czasem poprawka jest więc o rzędy wielkości
mniejsza niż sam szczytowy \(\gamma-1\sim6\times10^{-6}\) sugerowałby.
Efekt jest zatem realny i poprawnie policzony, ale wielokrotnie poniżej
szumu Monte Carlo próby (\(N=1000\)): wykresy/statystyki z punktu M,
wygenerowane przed tą poprawką, pozostają aktualne — ponowny bieg
produkcyjny nie zmieniłby żadnej cyfry w granicach precyzji, w jakiej są
raportowane. `positronium_validation`: \(33/33\) bez regresji.

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
eksperymentu, numer ekranu, literę `a`/`b`, numer padu i opis, rozdzielone
znakami `_`, na przykład:

```text
distributions/1_1_b_1_crem_collapse_time.pdf
distributions/1_2_a_2_dipole_coupling_vs_hyperfine.pdf
distributions/3_1_b_3_energy_loss_cross_section.pdf
distributions/4_1_b_1_differential_cross_section.pdf
```

Litera rozróżnia charakter danych: `a` — rozkład danych WEJŚCIOWYCH
(cechy przygotowanego/wylosowanego stanu, zanim zadziała na nie dynamika:
np. sprzężenie dipol-dipol przygotowanej pary, próbkowana energia
zderzenia, parametr zderzenia, ustawienie dipoli), `b` — rozkład danych
WYJŚCIOWYCH (cokolwiek symulacja faktycznie wyprodukowała lub zmierzyła:
czasy kolapsu, przekroje czynne, residua praw zachowania, widma fotonów).
Klasyfikacja jest jawna w kodzie przy każdym wykresie z osobna
(`root_export::NamedPad::inputOrOutput`), nie domyślna — świadoma decyzja
przy dodawaniu każdego panelu, nie odziedziczony domyślny wybór.
Zrzut ekranu trybu wizualnego (`{N}_1_1_visual_simulation.pdf`) nie
przechodzi przez ten mechanizm i nie ma litery — to migawka całego stanu
symulacji, nie rozkład statystyczny.

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
| `--ground-state-floor` | wyłączone | **Eksperyment, nie część modelu.** Druga próba domknięcia tej samej luki co `--zpf`, przeciwną metodą: zamiast mechanizmu importuje **jeden** fakt kwantowy — drabina Bohra kończy się na \(n=1\) — i odmawia emisji, która związałaby parę ciaśniej niż stan podstawowy. Bramkuje tempo hazardu (nie odrzuca pojedynczych fotonów: bankowanie i ponawianie zakleszcza się, co ten plik odnotowuje osobno), klamruje wszystkie trzy miejsca zapisu `specificEnergy` i przycina ostatni kwant tak, by wylądował dokładnie na podłodze. **Zatrzymuje spiralę na \(a_{Ps}\), ale likwiduje obserwablę czasu kolapsu — patrz niżej.** |
| `--zpf`, `--zpf-band` | `0` (wyłączone) | **Eksperyment, nie część modelu.** Klasyczne pole punktu zerowego elektrodynamiki stochastycznej: losowe fale płaskie o widmie \(\rho(\omega)=\hbar\omega^3/2\pi^2c^3\), 64 mody o równej energii, orientacje i fazy z ziarna `--seed`. `--zpf` skaluje **amplitudę** (1 = poziom fizyczny, moc pochłaniana rośnie jak kwadrat), `--zpf-band lo,hi` ustala pasmo w jednostkach częstości orbitalnej pary (domyślnie `0.3,3`). To jest fluktuacyjna połowa pary fluktuacja–dyssypacja; dyssypacyjną, czyli reakcję promieniowania, model ma od zawsze. Wchodzi w te same trzy miejsca co pole jednorodne, ale próbkowane osobno dla każdej cząstki, bo zależy od położenia i czasu. **Nie odtwarza stanu podstawowego SED — patrz niżej.** |
| `--external-field` | brak (pytanie na starcie) | Jednorodne zewnętrzne pole magnetyczne w mikroteslach; `0` wyłącza. Orientacja jest losowana izotropowo z ziarna `--seed`, więc odtwarza się razem z resztą przebiegu, i jest wypisywana na starcie. Gdy opcji nie podano, a przebieg jest interaktywny, program pyta o to **przed wszystkimi pozostałymi pytaniami** i oferuje 50 µT (skala pola ziemskiego). Przebieg wsadowy z podanym `--mode` i `--phenomenon` nigdy nie pyta i domyślnie nie ma pola. Pole wchodzi w sumę sił chwilowych, w sumę sił retardowanych oraz w pole lokalne widziane przez obie cząstki, przez co obejmuje precesję Thomasa-BMT. Przy 50 µT tempo cyklotronowe \(eB/m\) wynosi 8,8·10⁶ rad/s wobec tempa orbitalnego rzędu 3·10¹⁵ rad/s, więc orbita pozostaje nietknięta, a widocznym kanałem jest precesja dipoli — około 3·10⁻⁴ rad w ciągu 35 ps kolapsu. |
| `--level` | `1` | **Tryb wzbudzony.** Główna liczba kwantowa, na której przygotowywana jest para związana: \(a_n=n^2a_{\rm pary}\), pasmo prędkości stycznej niezmienione względem prędkości kołowej przy tej separacji. Energia fotonu podąża wtedy za odstępem poziomów \(\Delta E(n\to n-1)\), dopóki \(n\ge2\), a poniżej wraca do \(\hbar\omega\). Czas kolapsu rośnie jak \(n^6\), więc dla \(n\ge3\) trzeba podnieść `--crem-wallclock-budget-s`. Wartość `1` jest bit-identyczna z zachowaniem sprzed wprowadzenia opcji. Patrz „Tryb wzbudzony" wyżej. |
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

**Rozstrzygnięte: to niezgodność siły z energią, a nie błąd całkowania.**
Wcześniejsza wersja tego akapitu zostawiała pytanie otwarte, zauważając tylko,
że lustrzana tożsamość wobec rezerwuaru pola nie odróżnia fizycznej wymiany od
błędu schematu. Rozdziela je natomiast prosty eksperyment, którego wtedy nie
wykonano: **przełączyć siłę, a nie tolerancję.**

Silnik ma diagnostyczny przełącznik `useRetardedExternalForces`. Przy `true`
(produkcja) całkowana jest pełna, retardowana siła Liénarda–Wiecherta; przy
`false` — sama siła Coulomba–Darwina, czyli dokładnie ta, której potencjałem
`conservativeParticleEnergy()` **jest**. Zmierzony dryf energii zachowawczej
przez jeden obieg przy \(a_{Ps}\), z wyłączoną reakcją promieniowania i
zerowymi dipolami:

| tolerancja | retardowana (produkcja) | Coulomb+Darwin |
|---|---|---|
| \(10^{-5}\) | \(2{,}314\cdot10^{-6}\) | \(5{,}131\cdot10^{-7}\) |
| \(10^{-6}\) | \(2{,}314\cdot10^{-6}\) | \(5{,}131\cdot10^{-7}\) |
| \(10^{-7}\) | \(1{,}280\cdot10^{-6}\) | \(1{,}283\cdot10^{-7}\) |
| \(10^{-8}\) | \(1{,}827\cdot10^{-6}\) | \(6{,}415\cdot10^{-8}\) |
| \(10^{-9}\) | \(1{,}604\cdot10^{-6}\) | \(3{,}208\cdot10^{-8}\) |

Dopasowana para siła–energia **zbiega monotonicznie**, o czynnik \(2\) na
dekadę tolerancji, bez śladu podłogi. Para niedopasowana staje w miejscu.
Podłogi nie ma więc w integratorze — jest w tym, że mierzy się zachowanie
wielkości, która nie jest całką ruchu całkowanej siły.

*Którego rzędu jest brakujący człon — zmierzone, nie założone.* Podłoga
przemieciona po \(\beta\) (promień od \(a_{Ps}\) do \(a_{Ps}/16\), tolerancja
zacieśniona aż do nasycenia):

| \(r/a_{Ps}\) | \(\beta\) | podłoga | podłoga\(/\beta^3\) | nachylenie lokalne |
|---|---|---|---|---|
| \(1\) | \(7{,}297\cdot10^{-3}\) | \(1{,}604\cdot10^{-6}\) | \(4{,}13\) | — |
| \(1/2\) | \(1{,}032\cdot10^{-2}\) | \(4{,}492\cdot10^{-6}\) | \(4{,}09\) | \(2{,}97\) |
| \(1/4\) | \(1{,}460\cdot10^{-2}\) | \(1{,}190\cdot10^{-5}\) | \(3{,}83\) | \(2{,}81\) |
| \(1/8\) | \(2{,}064\cdot10^{-2}\) | \(3{,}587\cdot10^{-5}\) | \(4{,}08\) | \(3{,}18\) |
| \(1/16\) | \(2{,}919\cdot10^{-2}\) | \(9{,}485\cdot10^{-5}\) | \(3{,}81\) | \(2{,}81\) |

Wykładnik \(2{,}94\pm0{,}15\), a iloraz \(\text{podłoga}/\beta^3\) stały co do
\(7\%\) na czterokrotnym zakresie \(\beta\). Podłoga wynosi więc
\(\approx4\beta^3\) — a \(O(\beta^3)\) to **dokładnie pierwszy rząd, którego
lagranżjan Darwina nie zawiera**: Darwin jest zupełny do \(O(v^2/c^2)\)
włącznie, a następny rząd retardacji nie daje się już zapisać jako potencjał
dwuciałowy, bo to ten rząd, w którym energia realnie przechodzi do pola.

To wyjaśnia komplet objawów naraz. Niezależność od precyzji mantysy — bo to
nie zaokrąglenia. Niezależność od tolerancji — bo to nie błąd kroku.
Deterministyczność — bo to człon fizyczny o ustalonej wielkości. A nawet
**wzrost** przy zacieśnianiu tolerancji z \(10^{-7}\) do \(10^{-9}\): podłoga
nie jest granicą monotoniczną, tylko obwiednią oscylacji zależnej od fazy, więc
inna sekwencja kroków trafia w inne maksimum w tej samej obwiedni.

Praktyczny wniosek stoi bez zmian, ale ma teraz uzasadnienie zamiast samej
obserwacji: próg leży niecałe dwa razy powyżej sygnału
\(2{,}6\cdot10^{-6}\) na obieg i nie usuwa go ani rząd, ani tolerancja, ani
precyzja, więc odejmowanie przebiegu tła jest konieczne — a jest przy tym
**właściwym** narzędziem, bo przebieg tła niesie tę samą obcinkę Darwina w tej
samej fazie, więc odejmowanie ją kasuje, zostawiając stratę fizyczną.

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
wyłącznie obcięciem numerycznym, a wynik od niego silnie zależy. Poprawne
odtworzenie równowagi SED wymaga samouzgodnionej odpowiedzi cząstki na pełne
widmo, a nie dołożenia losowego pola do gotowych równań ruchu; równowaga jest
delikatnym skasowaniem, nie efektem rzędu wiodącego.

> **Sprostowanie dwóch zdań tej sekcji.** Stało tu wcześniej, że „ta
> implementacja nie daje odpowiedzi fizycznej", a kilka akapitów niżej — że
> ekstrapolacja ustalonego trendu czyni bieg przy pełnym odcięciu Comptona
> testem „obarczonym niemal pewnym wynikiem". Oba naraz nie mogą być prawdą.
> Skoro trend z górną krawędzią pasma jest monotoniczny i pewny, a widmo
> fizyczne odpowiada granicy pasm coraz szerszych, to implementacja **daje**
> odpowiedź w tej granicy — brak stanu związanego — tylko nie daje jej przy
> żadnym pojedynczym obcięciu. Werdykt zostaje, uzasadnienie zostało zaniżone.

**Dlaczego to jest trudne — zmierzone, a nie tylko nazwane.** Powyższe „delikatne
skasowanie" było dotąd przypuszczeniem bez liczby, a wszystkie testy tej sekcji
mierzą *skutki* (czas kolapsu, promień końcowy, `trajectory: FAIL`), nigdy samą
wielkość, od której zależy mechanizm. Równowaga SED to warunek
\(\langle P_{\rm pochł}\rangle=\langle P_{\rm wypr}\rangle\), więc zmierzono
wprost pracę pola na obu ładunkach,
\(\int(q_1\mathbf v_1\!\cdot\!\mathbf E_1+q_2\mathbf v_2\!\cdot\!\mathbf E_2)\,dt\),
przeciw \(\int P_{\rm Larmor}\,dt\) — orbita kołowa przy \(a_{Ps}\), amplituda
fizyczna (`scale=1`), osiem orbit:

| pasmo | mody | \(\langle P_{\rm pochł}\rangle\) [W] | \(\langle P_{\rm wypr}\rangle\) [W] | stosunek |
|---|---|---|---|---|
| \([0{,}3;3]\) | \(64\) | \(2{,}45\cdot10^{-6}\) | \(1{,}14\cdot10^{-8}\) | \(+214\) |
| \([0{,}3;3]\) | \(256\) | \(2{,}34\cdot10^{-6}\) | \(1{,}15\cdot10^{-8}\) | \(+204\) |
| \([0{,}5;2]\) | \(64\) | \(-1{,}19\cdot10^{-6}\) | \(1{,}16\cdot10^{-8}\) | \(-102\) |
| \([0{,}5;2]\) | \(256\) | \(-6{,}04\cdot10^{-6}\) | \(1{,}19\cdot10^{-8}\) | \(-506\) |
| \([0{,}9;1{,}1]\) | \(64\) | \(-9{,}02\cdot10^{-7}\) | \(1{,}17\cdot10^{-8}\) | \(-77\) |
| \([0{,}9;1{,}1]\) | \(256\) | \(-2{,}98\cdot10^{-6}\) | \(1{,}19\cdot10^{-8}\) | \(-250\) |
| \([0{,}3;30]\) | \(64\) | \(1{,}28\cdot10^{-5}\) | \(1{,}09\cdot10^{-8}\) | \(+1180\) |
| \([0{,}3;30]\) | \(256\) | \(-6{,}90\cdot10^{-6}\) | \(1{,}19\cdot10^{-8}\) | \(-578\) |

Te liczby są jednak **fluktuacją, nie absorpcją systematyczną**, i to trzeba
powiedzieć od razu, bo pierwsza wersja tego akapitu odczytała je odwrotnie.
Rozstrzyga skalowanie po amplitudzie pola: część fluktuacyjna pracy jest
**liniowa** w polu (prędkość jest w pierwszym przybliżeniu niezaburzona), a
człon systematyczny — ten, który wchodzi do bilansu SED — pochodzi z korelacji
pola z prędkością, którą samo wywołało, więc jest **kwadratowy**. Zmierzone,
pasmo \([0{,}9;1{,}1]\), \(8\) ziaren:

| skala | średnia \(P_{\rm pochł}\) | SD | nachylenie średniej | \(r_{\rm końc}/r_0\) |
|---|---|---|---|---|
| \(1\) | \(-2{,}181\cdot10^{-6}\) | \(1{,}909\cdot10^{-6}\) | — | \(0{,}9951\) |
| \(2\) | \(-4{,}381\cdot10^{-6}\) | \(3{,}900\cdot10^{-6}\) | \(1{,}006\) | \(0{,}9904\) |
| \(4\) | \(-8{,}992\cdot10^{-6}\) | \(8{,}223\cdot10^{-6}\) | \(1{,}037\) | \(0{,}9809\) |
| \(8\) | \(-1{,}890\cdot10^{-5}\) | \(1{,}874\cdot10^{-5}\) | \(1{,}071\) | \(0{,}9606\) |

Nachylenie \(1{,}0\), nie \(2\). To, co wyżej wygląda na „średnią", jest jedną
realizacją przemnożoną przez amplitudę — te same ziarna dają to samo pole przy
każdej skali. Stąd też bierze się rzekoma niestabilność znaku między pasmami:
to własność próbki, nie fizyki.

*Wyodrębnienie członu systematycznego.* Dla ustalonej realizacji część liniowa
jest deterministyczna, więc da się ją usunąć **tożsamościowo**: przesunięcie
faz wszystkich modów o \(\pi\) neguje pole dokładnie, a w sumie
\(P_{\rm pochł}(+\mathbf F)+P_{\rm pochł}(-\mathbf F)=2b\,s^2\) znikają
wszystkie rzędy nieparzyste. Kasowanie działa (redukcja o czynnik \(\approx7{,}5\)
przy \(s=2\)). Wynik na \(16\) ziarnach, w jednostkach mocy promieniowanej:

```
średnia   -2,98   SEM 1,61
mediana   -0,36
średnia ucinana 25%   -1,29
zakres    -22,54 .. +2,34      ujemnych 9 z 16
```

Niezależne dopasowanie \(P=as+bs^2\) na pięciu amplitudach i \(10\) ziarnach
daje zgodne \(-3{,}61\pm1{,}84\).

**Co z tego wynika, a co nie.** Wynika **rząd wielkości**: człon systematyczny
jest tego samego rzędu co moc promieniowana, \(|b|\sim P_{\rm wypr}\). To jest
właśnie reżim „delikatnego skasowania" — teraz z liczbą, a nie jako
przypuszczenie. **Nie wynika** ani znak, ani wartość: obie niezależne oceny
leżą \(\approx1{,}9\sigma\) od zera, a mediana i średnia różnią się
ośmiokrotnie, więc rozkład jest ciężkoogonowy i zdominowany przez dwa
odstające ziarna. Doprowadzenie SEM poniżej \(0{,}5\,P_{\rm wypr}\)
wymagałoby \(\approx165\) ziaren przy założeniu gaussowskim, którego ten ogon
nie spełnia.

> **Sprostowanie mojego własnego wniosku.** Zapisałem tu wcześniej, że „wąskim
> gardłem jest okno uśredniania, **nie** górna krawędź pasma". To skleja dwa
> pytania o różnych odpowiedziach. Krawędź pasma decyduje o **wyniku** — trend
> \(43{,}0\to152{,}5\) ps jest systematyczny i monotoniczny. Okno uśredniania
> decyduje o **mierzalności** warunku równowagi. Drugie zdanie się broni,
> przeczenie w pierwszym było błędne.

Uczciwy stan rzeczy jest więc taki: **ta sekcja nie zmierzyła warunku
równowagi SED i przy tym nakładzie go nie zmierzy** — ustaliła natomiast, że
absorpcja systematyczna i strata promienista są w rezonansie tego samego
rzędu. To nie jest to samo co wykazanie, że równowagi nie ma.

Werdykt sekcji zostaje bez zmian — ucieczka orbity przy szerokich pasmach jest
zmierzona i realna, a pięć prób kończących się `trajectory: FAIL` na trzech
ziarnach nadal obowiązuje. Zmienia się diagnoza przyczyny.

*Umiejscowienie w literaturze — i granica tego, co tu twierdzę.* Sekcja
opierała się na pojęciu elektrodynamiki stochastycznej, a bibliografia nie
miała ani jednej pozycji z tej dziedziny; czytelnik nie miał więc jak
sprawdzić ram ani znaleźć miejsca, gdzie to samo pytanie jest dyskutowane.
Do katalogu dopisano sześć pozycji, wszystkie z metadanymi zweryfikowanymi
w Crossref przez `make references-check` (\(11\) DOI sprawdzonych, \(0\)
niezgodności):

| klucz | pozycja |
|---|---|
| `boyer_1975` | T. H. Boyer, *Random electrodynamics: The theory of classical electrodynamics with classical electromagnetic zero-point radiation*, Phys. Rev. D **11**, 790–808 (1975) |
| `cole_zou_2003` | D. C. Cole, Y. Zou, *Quantum mechanical ground state of hydrogen obtained from classical electrodynamics*, Phys. Lett. A **317**, 14–20 (2003) |
| `boyer_2003_comment` | T. H. Boyer, *Comments on Cole and Zou's Calculation of the Hydrogen Ground State in Classical Physics*, Found. Phys. Lett. **16**, 613–617 (2003) |
| `milonni_2003_comment` | P. W. Milonni, *Comment on Cole and Zou's Classical Computations of the Hydrogen Ground State*, Found. Phys. Lett. **16**, 619–621 (2003) |
| `nieuwenhuizen_liska_2015` | T. M. Nieuwenhuizen, M. T. P. Liska, *Simulation of the hydrogen ground state in stochastic electrodynamics*, Phys. Scr. **T165**, 014006 (2015) |
| `nieuwenhuizen_liska_2015b` | T. M. Nieuwenhuizen, M. T. P. Liska, *…-2: Inclusion of Relativistic Corrections*, Found. Phys. **45**, 1190–1202 (2015) |

Weryfikacja od razu złapała błąd, który powstałby przy cytowaniu z pamięci:
tytuł pracy Cole'a i Zou kończy się na „from classical **electrodynamics**",
nie „from classical physics".

**Czego te pozycje tutaj nie uzasadniają.** Crossref nie zwraca dla nich
abstraktów, a prace nie zostały tu przeczytane, więc **nie parafrazuję ich
wniosków** i żaden wynik tej sekcji się na nich nie opiera. Twierdzę
wyłącznie to, co wynika z samych tytułów: że Cole i Zou ogłosili wynik
pozytywny dla stanu podstawowego wodoru, że ukazały się dwa opublikowane
**komentarze** do tego rachunku, oraz że problem był potem symulowany
ponownie, w drugiej pracy z poprawkami relatywistycznymi. Innymi słowy —
pytanie, które ta sekcja bada, jest w literaturze **sporne**, a nie
zaniedbane; kto ma rację, nie jest tu rozstrzygane ani zakładane.

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

### Alternatywa dla `--zpf`: podłoga emisyjna w stanie podstawowym

Skoro `--zpf` nie dostarczył mechanizmu zatrzymującego spiralę, spróbowano
drogi przeciwnej. `--zpf` szukał **mechanizmu** — równowagi fluktuacyjno-
dyssypacyjnej z realnym polem. `--ground-state-floor` żadnego mechanizmu nie
proponuje: importuje **jeden** fakt kwantowy, że pod \(n=1\) nie ma stanu, do
którego foton mógłby parę zostawić, i odmawia takiej emisji.

**To jest domknięcie, nie wyprowadzenie**, i tak trzeba czytać każdy wynik z
tą flagą. Nie tłumaczy, dlaczego niższego stanu nie ma — stwierdza to.
W szczególności para kończy przy \(6{,}802847\) eV **dlatego, że tę liczbę
bramce podano**: \(a_{Ps}\) jest zdefiniowane jako \(\hbar^2/(\mu k)\), więc
\(k/(2a_{Ps})\) jest algebraicznie \(\mu k^2/2\hbar^2\), czyli sam wzór
Bohra-Rydberga; obie drogi zgadzają się do \(1{,}8\cdot10^{-16}\), bo to ten
sam wzór. Żadna dynamika CREM w to nie wchodzi.

*Co daje.* Spirala **zatrzymuje się** — \(0\) kolapsów na \(200\) wobec
\(186/200\) bez bramki, \(S(t)=1\) przez \(4398\) ps. Przy budżecie
zegarowym pozwalającym dojść do końca, \(7\) z \(8\) trajektorii siada na

```
E/E_gs = 1,000001     r/a_Ps = 0,999999
```

czyli na stanie podstawowym co do sześciu cyfr, z rozmaitych warunków
początkowych. Ósma po \(16{,}7\) ns wciąż schodzi. Tego `--zpf` nie osiągnął
ani razu: tam szersze pasma zawsze pompowały orbitę ku ucieczce.

*Dwie usterki wykryte pomiarem, obie niewidoczne w teście „czy się
zatrzymuje".* Warto je zapisać, bo pokazują, że samo zatrzymanie niczego nie
dowodzi:

1. **Przeciek.** Bramkowanie hazardu fotonowego nie wystarcza — ścieżka
   stochastyczna ma drugi kanał, `elements.specificEnergy += deltaEnergyPerOrbit`,
   omijający maszynerię fotonową. Trajektoria zeszła \(19\%\) za głęboko, na
   \(0{,}8407\,a_{Ps}\). Zamknięte klamrą przy każdym zapisie.
2. **Utknięcie powyżej podłogi.** Odrzucanie zbyt dużego fotonu zostawiało
   parę na \(1{,}2478\,a_{Ps}\): kwant korespondencyjny \(\hbar\omega_{\rm orb}\)
   wynosi tam \(1{,}43\,E_{gs}\), a miejsca zostało \(0{,}20\,E_{gs}\), więc
   żaden foton już się nie mieścił. Właściwą regułą nie jest „odmów", tylko
   „ostatnie przejście jest skokiem poziomowym": para o krok nad stanem
   podstawowym emituje \(E(n)-E(1)\), nie \(\hbar\omega_{\rm orb}\). Kwant jest
   więc przycinany do odległości od podłogi — ta sama reguła, którą plik
   stosuje już dla \(n\ge2\).

#### Przewiązanie anihilacji i kwantowanie momentu pędu

Sama podłoga energetyczna zatrzymuje spiralę, ale zostawia parę bez punktu
końcowego: nie docierając do bariery Comptona, nie anihiluje, więc wszystkie
trajektorie są cenzurowane i program sygnalizuje to kodem \(2\). Anihilacja
została więc przewiązana z bariery do stanu podstawowego.

**Przewiązanie nie przemyca tempa anihilacji, bo model żadnego nie ma — ani
przed, ani po.** W obu wariantach para anihiluje **deterministycznie** w
chwili dotarcia; zmienia się dokąd dociera, nie czym jest to zdarzenie.
Sprawdzone zostały trzy możliwe źródła tempa i żadne nie działa: kanał
kontaktowy nie istnieje (klasyczna orbita Keplera z \(L\ne0\) nigdy nie
osiąga \(r=0\), a z podłogą orbita jest dodatkowo zamrożona), hazard E1 w
stanie podstawowym opisuje kwanty \(13{,}6\) eV a nie dwa po \(511\) keV,
a \(\Gamma_{\rm para}\) z QED byłoby importem wyniku. Raportowany czas pod
tą flagą jest więc **czasem kaskady**, a nie czasem życia anihilacyjnego.

*Zysk, który zostaje niezależnie od reszty: linia anihilacyjna.* Bez
przewiązania para anihilowała ze stanu \(375\times\) za głęboko związanego i
linia wychodziła \(510{,}4\) keV — przesunięta o \(0{,}6\) keV. Po
przewiązaniu:

```
annihilation W    1021,99 keV wobec 1022 keV      (-0,00067%)
W/2 (linia 2g)     510,996 keV wobec 510,999 keV
```

czyli przesunięcie \(\approx3{,}4\) eV, dokładnie skala energii wiązania
pozytonium. Poprzednia wartość była **\(160\times\) za duża** — nie dlatego,
że kinematyka była zła, tylko dlatego, że punkt końcowy był zły.

*Kwantowanie momentu pędu — dlaczego okazało się konieczne.* Z samą podłogą
energetyczną zmierzone średnie wiązanie terminalne wyszło \(7{,}511\) eV
wobec podłogi \(6{,}803\) eV, czyli \(10\%\) za głęboko. Powód: podłoga na
ENERGII ustala półoś wielką i nie mówi nic o peryapsis, a orbita mimośrodowa
przy \(a=a_{Ps}\) ma peryapsis \(a(1-e)\), które dla dużego \(e\) sięga
bariery — i część trajektorii kończyła tam.

Zamyka to Bohr-Sommerfeld: azymutalna liczba kwantowa biegnie \(k=1..n\), więc
\(n=1\) dopuszcza wyłącznie \(k=1\), czyli \(L=\hbar\) i mimośród **zero**.
Podłoga \(L\ge\hbar\) czyni ze stanu podstawowego orbitę kołową o peryapsis
równym \(a_{Ps}\), dla której bariera jest nieosiągalna. Zmierzone:

| | podłoga na E | + podłoga na L |
|---|---|---|
| promień terminalny | \(95\,859\) fm \((0{,}906\,a_{Ps})\) | \(105\,835\) fm \(=a_{Ps}\) |
| wiązanie terminalne | \(7{,}511\) eV | \(6{,}80285\) eV \(=E_{gs}\) |

Oba trafione co do cyfry, kanał barierowy zamknięty.

*Czas kaskady.* Przy domyślnym seedingu mediana wynosi \(0\) ps i **jest to
poprawne**: pasmo startowe jest wyśrodkowane na \(a_{Ps}\), więc połowa par
rodzi się w stanie podstawowym i nie ma kaskady do przejścia. (Seeding
przepuszczono przez obie klamry, żeby para wylosowana *poniżej* podłogi
startowała na niej, a nie wewnątrz niej; startowania *na* podłodze to nie
zmienia i zmienić nie może.) Kaskadę mierzy się dopiero od \(n\ge2\), i tam
\(60\) trajektorii z \(n=2\) nie ukończyło jej w \(50\,832\) ps czasu
symulacji — mamy więc **dolne ograniczenie \(>50\) ns**, nie wartość.

Powolność jest zrozumiała i wynika z samej podłogi na \(L\): orbity
mimośrodowe promieniują nieporównanie mocniej (czynnik
\((1+e^2/2)/(1-e^2)^{5/2}\) rozbiega się przy \(e\to1\)), a kwantowanie
\(L\) trzyma je okrągłymi. Kaskada zwalnia przez to o rząd wielkości wobec
naiwnego oszacowania z prawa \(a^3\), które dawałoby \(\approx7{,}6\) ns.

#### Kwantowanie spinu, i czy 2γ/3γ daje wskazania do modelu anihilacji

Liczba fotonów anihilacji nie ma dziś w modelu **żadnej** treści dynamicznej:
kanał wybiera flaga `--phenomenon`, a `annihilationPhotonEnergiesFor` wyprowadza
z tego wyboru \(2\gamma\) albo \(3\gamma\) zakodowaną regułą. Ustawienie
momentów służy do przygotowania pary zgodnie z wybranym kanałem, nie do jego
wyznaczenia.

Model zawiera natomiast **jeden dokładny mechanizm** różnicujący kanały, i to
o właściwej strukturze. Moc M1 liczona jest koherentnie z sumy
\(\mathbf m=\boldsymbol\mu_1+\boldsymbol\mu_2\), więc przy \(\cos=-1\)
kasuje się **ściśle do zera** (zmierzone: `0,0000e+00`, nie „małe"), a przy
\(\cos=+1\) interferuje konstruktywnie. To ten sam kształt, co reguła wyboru
w QED, gdzie \(C|n\gamma\rangle=(-1)^n\) wobec \(C|Ps\rangle=(-1)^{L+S}\)
zamyka wiodący kanał dla jednego stanu spinowego i spycha go o rząd
\(\alpha\) wyżej.

Przy losowaniu kąta z przedziału mechanizm się jednak rozmywał: udział M1
wychodził \(1{,}875\cdot10^{-3}\) dla para wobec \(1{,}903\cdot10^{-3}\) dla
orto, różnica zgodna z zerem. Trzeci floor usuwa przedział — \(S=0\) i
\(S=1\) to stany dokładne, więc \(\cos=\pm1\) ściśle.

*Zmierzony skutek, po sprostowaniu.* Pierwsza wersja tego akapitu głosiła, że
człon dipol-dipol „przestaje uśredniać się do zera", na podstawie przebiegu
dającego \(-6{,}053\cdot10^{-9}\) keV dla para i \(+6{,}411\cdot10^{-9}\) keV
dla orto. **To była nadinterpretacja.** Oba przebiegi szły z tego samego
ziarna, więc miały te same orientacje względem \(\hat{\mathbf L}\), a różnił
je wyłącznie znak \(\cos\). Niemal dokładnie przeciwne wartości to **jedna
próbka wielkości zerośredniej z odwróconym znakiem**, a nie dwa przesunięcia.

Sprawdzone wprost, \(20\,000\) realizacji ze skwantowanym spinem
(\(\cos=\pm1\) ściśle, losowa tylko orientacja względem \(\hat{\mathbf L}\)):

| \(r\) | średnia \(U_{dd}\) [eV] | SD [eV] | odchylenie od zera |
|---|---|---|---|
| \(a_{Ps}\) | \(-1{,}15\cdot10^{-8}\) | \(2{,}04\cdot10^{-5}\) | \(0{,}08\sigma\) |
| \(a_{Ps}/100\) | \(+0{,}255\) | \(20{,}4\) | \(1{,}77\sigma\) |
| \(1000\) fm | \(-0{,}182\) | \(24{,}0\) | \(1{,}07\sigma\) |
| \(282\) fm | \(+3{,}46\) | \(1071\) | \(0{,}46\sigma\) |
| bariera | \(-81{,}3\) | \(3316\) | \(3{,}5\sigma\), SD \(40\times\) większe |

*Zmierzone osobno dla obu kanałów, bo poprzednia wersja mierzyła tylko para i
zakładała, że orto jest jej odbiciem.* Odbiciem jest, i to **dokładnym**:
\(\max|U_{\rm para}+U_{\rm orto}|/|U_{\rm para}|=0{,}000\) bitowo, przy każdym
promieniu. Ale wartość oczekiwana zachowuje się różnie w różnych reżimach
(\(N=200\,000\) orientacji na wiersz):

| \(r\) | para: średnia [eV] | odchylenie od zera |
|---|---|---|
| \(a_{Ps}\) | \(+3{,}93\cdot10^{-9}\) | \(0{,}09\sigma\) |
| \(a_{Ps}/10\) | \(-5{,}90\cdot10^{-5}\) | \(1{,}30\sigma\) |
| \(a_{Ps}/100\) | \(-2{,}03\cdot10^{-2}\) | \(0{,}45\sigma\) |
| \(1000\) fm | \(+8{,}01\cdot10^{-3}\) | \(0{,}15\sigma\) |
| \(282\) fm | \(-2{,}38\) | \(0{,}99\sigma\) |
| **bariera** | \(\mathbf{-59{,}1}\) | \(\mathbf{7{,}98\sigma}\) |

**Przy barierze Comptona wartość oczekiwana NIE znika**, i to nie jest
fluktuacja: pięć niezależnych ziaren daje \(-52{,}6\), \(-59{,}7\),
\(-63{,}3\), \(-51{,}3\), \(-67{,}6\) eV, każde \(6{,}9\)–\(9{,}2\sigma\) od
zera, średnio \(-58{,}9\) eV \(=-0{,}79\%\) potencjału kulombowskiego
(\(7449\) eV). Różnica kanałów wynosi więc \(118\) eV, czyli \(1{,}6\%\)
Coulomba.

Powód jest w regularyzacji. Dla czystego \(1/r^3\) kasowanie jest
tożsamościowe — przy \(\langle\cos^2\theta\rangle=1/3\) człon poprzeczny i
radialny znoszą się dokładnie — i to tłumaczy zera przy dużych promieniach.
Regularyzacja łamie to kasowanie i przy barierze zostawia realną resztkę.

*Odpowiedź, zawężona.* Różnica kanałów **nie** wchodzi przez budżet
energetyczny **przy promieniach rzędu \(a_{Ps}\)** — tam wartość oczekiwana
znika (\(0{,}09\sigma\)), więc człon daje wariancję, nie przesunięcie. Ale
**przy barierze wchodzi**, systematycznie, na poziomie \(1{,}6\%\) Coulomba.
Poprzednia wersja tego akapitu twierdziła, że znika wszędzie; to było za
mocne.

*Dlaczego mimo to nie widać tego w żadnym przebiegu.* Rozrzut wynosi
\(3316\) eV wobec średniej \(58{,}9\) eV, czyli \(\mathrm{SD}/|\text{średnia}|
=56{,}3\). Wykrycie efektu wymaga \(\approx3170\) trajektorii na \(1\sigma\) i
\(\approx28\,500\) na \(3\sigma\). Przebiegi produkcyjne po \(60\)–\(1000\)
zdarzeń są o rząd do dwóch poniżej progu, więc **każde porównanie kanałów w
tej dokumentacji jest zdominowane szumem** — co jest niezależnym wyjaśnieniem,
dlaczego mediany obu kanałów zawsze wychodziły takie same.

*Konsekwencja dla trzech podłóg — i skąd naprawdę bierze się ta resztka.*
Reżimy różnią się jakościowo, ale zbadane głębiej okazują się różnić inaczej,
niż wyglądało. Resztka przy barierze **nie jest efektem fizycznym, który
podłogi tłumią** — jest artefaktem regularyzacji, poza który podłogi model
wyprowadzają.

Regularyzacja pola dipolowego ma promień \(r_{\rm reg}=68{,}47\) fm i wykładnik
\(6\), a bariera Comptona leży przy \(r/r_{\rm reg}=2{,}82\), gdzie waga
\(w=0{,}998\). Resztka śledzi tę wagę:

| \(r/r_{\rm reg}\) | \(w(r)\) | średnia [eV] | średnia\(/U_{\rm coul}\) |
|---|---|---|---|
| \(8\) | \(0{,}999996\) | \(+0{,}29\) | \(1{,}1\cdot10^{-4}\) |
| \(4\) | \(0{,}99976\) | \(-0{,}20\) | \(-3{,}7\cdot10^{-5}\) |
| \(2\) | \(0{,}98462\) | \(-1251\) | \(-0{,}119\) |
| \(1{,}5\) | \(0{,}91929\) | \(-14703\) | \(-1{,}05\) |
| \(1\) | \(0{,}5\) | \(-167604\) | \(-7{,}97\) |

*Test rozstrzygający.* Gdyby resztka była fizyczna, nie zależałaby od
arbitralnego \(r_{\rm reg}\). Zależy — mierzone przy **stałym** promieniu
(bariera), \(\mathrm{SEM}\approx7{,}4\) eV:

| \(r_{\rm reg}\) [fm] | średnia przy barierze [eV] |
|---|---|
| \(34{,}2\) | \(+5{,}80\) (zgodne z zerem) |
| \(47{,}9\) | \(-0{,}21\) (zgodne z zerem) |
| \(\mathbf{68{,}5}\) (rzeczywiste) | \(\mathbf{-51{,}9}\) (\(7\sigma\)) |
| \(95{,}9\) | \(-424\) |
| \(136{,}9\) | \(-2964\) |

Zmiana \(r_{\rm reg}\) o czynnik \(2\) zmienia resztkę o czynnik \(57\) i
odwraca jej znak, a przy odsunięciu regularyzacji resztka **znika** — zgodnie
z tożsamościowym kasowaniem dla czystego \(1/r^3\). To nie jest wielkość
fizyczna.

Zwraca uwagę wzmocnienie: odchylenie wagi od jedynki o \(0{,}2\%\) przy
barierze daje efekt \(0{,}7\%\) potencjału kulombowskiego. Odpowiada za to
pochodna profilu — przy wykładniku \(6\) logarytmiczna pochodna wagi mnoży
odchylenie kilkukrotnie.

*Stąd właściwy podział.* Z trzema podłogami para siedzi przy \(a_{Ps}\), gdzie
\(r/r_{\rm reg}=1546\) i \(w=1\) do granic maszyny: kasowanie jest ścisłe, a
znikanie wartości oczekiwanej **fizyczne**. Bez podłóg para kończy przy
barierze, gdzie regularyzacja pracuje i nieznikanie jest **artefaktem**.
Podłogi nie ukrywają więc efektu kanałowego — trzymają model poza strefą, w
której jego własna regularyzacja taki efekt fabrykuje.

*Konsekwencja praktyczna — i naprawa.* Raportowany `terminal dipole-dipole`
niósł systematyczne \(-59\) eV pochodzące z doboru \(r_{\rm reg}\), na tle
\(3316\) eV rozrzutu; każde porównanie kanałów oparte na tej wielkości mierzyło
regularyzację, nie oddziaływanie dipol-dipol.

#### Czego warstwa sekularna nie widzi: człon dipolowy w punkcie zwrotnym

Rozrzut \(3332\) eV przy barierze, wobec \(7449\) eV kulombowskich, każe
zapytać, gdzie ta wielkość wchodzi, a gdzie jest pomijana. Prześledzone:

| miejsce | czy niesie człon dipolowy |
|---|---|
| siła w integratorze (`mutualForces`) | **tak** |
| `conservativeParticleEnergy` | **tak** |
| niezmiennik anihilacyjny \(W\) | **tak**, dodawany na końcu |
| `elements.specificEnergy` (element Keplera) | **nie**, świadomie |
| energia startowa `relativeEnergy` | **nie** |
| `regularizedPotentialEnergy` (punkty zwrotne) | **nie** |
| `osculatingPeriapsis` — **reguła zatrzymania** | **nie** |

Ostatni wiersz jest istotny. Trajektoria kończy się, gdy

```cpp
const double periapsis=osculatingPeriapsis(elements,attractionParameter);
if(periapsis<=comptonBarrierRadius || ...)
```

czyli decyduje **czysto keplerowskie** peryapsis, liczone z potencjału
zawierającego wyłącznie Coulomba — a w tym właśnie miejscu pominięty człon
dipolowy stanowi \(44{,}7\%\) potencjału kulombowskiego i jest **większy niż
raportowane wiązanie terminalne** (\(3332\) eV wobec \(2561\) eV).

*Test wrażliwości.* Dla orbity terminalnej z przebiegu produkcyjnego
(\(a=281{,}078\) fm, \(e=0{,}3123\)) rachunek bez dipola odtwarza próg
\(193{,}3035\) fm co do cyfry. Dodanie zmierzonego rozkładu \(U_{dd}\) przy
ustalonych \((E,L)\), \(20\,000\) orientacji:

| punkt zwrotny | udział |
|---|---|
| głębszy o \(>5\%\) | \(57{,}0\%\) |
| brak dozwolonej orbity przy tych \((E,L)\) | \(36{,}4\%\) |
| odsunięty o \(>5\%\) | \(5{,}1\%\) |
| **w granicach \(\pm5\%\)** | **\(1{,}5\%\)** |

Czyli dla \(98{,}5\%\) orientacji przewidywany punkt zwrotny leży poza
\(\pm5\%\) od tego, którym model się posługuje.

*Co to znaczy, a czego nie znaczy.* **Nie** znaczy, że trajektorie są źle
całkowane — integrator niesie pełną siłę i ruch jest poprawny. Znaczy, że
**warstwa sekularna** opisywała układ potencjałem, w którym w chwili
podejmowania decyzji brakowało członu porównywalnego z wiodącym.

To jest problem **wyłącznie terminalny**: człon dipolowy skaluje się jak
\(1/r^3\) wobec \(1/r\) Coulomba, więc przy \(a_{Ps}\) stanowi \(10^{-6}\)
potencjału. Gryzie dopiero tam, gdzie trajektoria się kończy — czyli dokładnie
tam, gdzie zapada decyzja.

*Naprawione: `dipoleAwarePeriapsis`.* Reguła zatrzymania rozwiązuje teraz
numerycznie

\[E=U_{\rm coul}(r)+U_{dd}(r)+\frac{L^2}{2r^2}\]

zamiast odwracać relację Keplera z potencjału bez dipola.

**Element nie wymagał zmiany**, i to jest sedno rozwiązania. Wcześniejsza
wersja tego dokumentu odkładała naprawę, twierdząc, że wymaga nadania elementom
niekeplerowskiego potencjału, bo osiem miejsc czyta je jako \(a=-K/2\varepsilon\),
\(e^2\) i okres. Okazało się to rozwiązywalne bez ich dotykania:
`elements.specificEnergy` jest inicjowane jako \(KE+U_{\rm coul}\) przy
\(a_{Ps}\), gdzie dipol to \(10^{-6}\) potencjału, a potem maleje wyłącznie o
wypromieniowaną energię — a zachowana jest suma \(KE+U_{\rm coul}+U_{dd}\).
Zatem

\[E_{\rm elem}=\big(KE+U_{\rm coul}\big)\Big|_{\rm teraz}+U_{dd}\Big|_{\rm teraz},\]

czyli **element już jest pełną energią całkowitą**; wadliwe było wyłącznie
*odwracanie* go z powrotem na promień. Dipol wchodzi więc tam, gdzie promień
jest faktycznie rozwiązywany, a znaczenie keplerowskie elementu zostaje
nietknięte.

Nowa funkcja wyraża też stan, którego stara reguła nie umiała: gdy człon
dipolowy jest odpychający na tyle, że żaden promień poniżej bieżącej orbity nie
jest dostępny, para **zatrzymuje się** i zwracany jest promień zatrzymania.
Wcześniej takie trajektorie ogłaszano jako „dotarły do bariery".

#### Wymuszony check bilansu kanału skwantowanego

Dwie gałęzie promieniste całkują **tę samą** kopertę Larmora po orbitach
pominiętych w checkpoincie: deterministyczna zdejmuje ją wprost, stochastyczna
ma dostarczyć ją w dyskretnych fotonach, przy hazardzie skalibrowanym tak, że
liczba fotonów skaluje się jak \(1/\hbar\omega\), a każdy niesie
\(\hbar\omega\). Jeśli to zachodzi, usunięta energia jest **niezmiennicza**
względem wyboru kwantu — czyli dokładnie teza, którą podważało wycofane 45%.

*Oczywisty check jest bezwartościowy.* Porównanie tego, co faktycznie
wyemitowano, z kopertą zmierzono: **1,87 / 2,49 / 2,92** na trzech ziarnach.
To nie jest wyciek. Trajektoria emituje **\(\approx2{,}6\) fotonu w całości**,
energia fotonu rośnie jak \(u^{3/2}\) wzdłuż spirali, w której \(u\) zmienia
się o dekady, a przebieg **kończy się na fotonie**. Suma jest zdominowana przez
własny ostatni wyraz i wyselekcjonowana na to, że jest duży — estymator o
wariancji \(O(1)\) z obciążeniem stopu, nie \(1/\sqrt N\). Odczytanie tego
rozrzutu jako fizyki byłoby błędem.

*Co da się wymusić.* Dwie tożsamości o **zerowej wariancji**, każda
sprawdzona osobno tym, że zawodzi po zepsuciu:

**1. Bilans koperty, na każdym checkpoincie, w produkcji.** Strona hazardowa
złożona z własnych zmiennych ścieżki emisji — `skipHazard` razy
`hazardReference` razy średni czynnik narastania \((1-s)^{-1}\) w skoku —
wobec koperty \(u\,[(1-J)^{-2/3}-1]\,\mu\). Algebraicznie równe, więc każda
różnica jest usterką kodu, nigdy fluktuacją. Zmierzone: **1 z dokładnością do
piętnastu cyfr**. Egzekwowane przy \(10^{-9}\), koszt kilku flopów, zero
ryzyka fałszywego alarmu. Checkpointy z przyciętym `jumpParameter` są
wykluczone z obu stron, bo tam obie opisują naprawdę różne energie.

**2. Tożsamość tabel harmonicznych, raz na proces.** Check (1) jest na nią
strukturalnie **ślepy**: hazard dzieli przez \(E_{\rm ref}/S(e)\), a liczba
mnoży przez \(S(e)\), więc \(S\) się skraca. Tymczasem \(S(e)\) jest
poprawne tylko wtedy, gdy równa się \(\langle1/n\rangle\) po widmie mocy —
bo wtedy i tylko wtedy harmonika losowana z rozkładu **zliczeń**
\((P_n/n)\) spełnia \(\langle n\rangle=1/S(e)\), co jest dokładnie tym, co
sprawia, że \(N\langle n\rangle E_{\rm ref}\) odtwarza \(P\) niezależnie od
rozłożenia mocy. Obie tabele powstały z **osobnych** dekompozycji
numerycznych, więc nic poza tą tożsamością ich nie wiąże. Zmierzone:
\(S(e)\langle n\rangle=1{,}00\pm0{,}03\), najgorzej \(1{,}065\) przy
\(e=0{,}95\); pas \(0{,}10\), zakres do \(e=0{,}97\).

Check (2) stoi przy tabelach, a nie w zestawie walidacyjnym, bo
`crem_collapse.hpp` siedzi za produkcyjnym `#ifndef` i **nie jest w ogóle
kompilowany** do walidatora.

*Czego szukałem i nie znalazłem.* Pierwszy pomiar dał deficyt \(16{,}5\%\)
— i był to **błąd w mojej sondzie**, nie w modelu: pominąłem czynnik
narastania energii fotonu w trakcie skoku. Rozpoznane po tym, że wyszedł
identyczny na trzech ziarnach (0,835318 / 0,835354 / 0,835283), co dla wielkości
statystycznej jest niemożliwe. Analitycznie
\(2(1-(1-J)^{1/3})/((1-J)^{-2/3}-1)=0{,}8353\) daje \(J=0{,}2998\), czyli
`maximumJumpParameter` = 0,30 — sonda sama wskazała, czego jej brakuje.

*Trop w sprawie 45%.* Bilans jest domknięty, więc te 45% **nie są wyciekiem
energii**. Przy \(\approx2{,}6\) fotonu na trajektorię każdy zdejmuje
\(\sim2\times\) bieżące wiązanie (\(\hbar\omega_{\rm orb}/|E|=2/n\)), a
przebieg kończy się, gdy periapsis dosięgnie bariery — więc o stanie końcowym
decyduje **przeskok** ostatniego fotonu za barierę. Większy kwant to skok z
dalszej pozycji. Test bezpośredni: z poprawką mediana spada
\(158{,}8\to26{,}6\) i \(87{,}1\to39{,}5\) ps, a fotonów **przybywa**
(\(21\to29\), \(23\to27\)) — i wiązanie terminalne idzie **głębiej**
(\(2{,}667\to3{,}883\) keV). Głębsze wiązanie przy mniejszym kwancie to
sygnatura przeskoku, nie wycieku. To ziarnistość reguły stopu, nie budżet
energii. Nie zamyka sprawy, ale przenosi ją z bilansu na regułę stopu.

#### Reguła stopu: co naprawdę kończy trajektorię

Bilans energii jest domknięty (poprzednia sekcja), więc czułość na kwant musi
siedzieć w regule stopu. Siedzi, i obraz jest ostrzejszy, niż się spodziewałem.

**Bariera Comptona kończy mniejszość przebiegów.** Zmierzone na 39 kolapsach:

```
stopped by             9 Compton barrier (23.1%), 30 retardation limit (76.9%)
landed at periapsis    1.85 r* (median), range 0.313-16.5 r*
period/light-crossing  66.1 at the stop (median), lowest 36.1 against 150
```

Trzy czwarte trajektorii kończy na warunku \(T/t_{\rm light}\le150\) — a ten
próg jest, wedle własnego uzasadnienia w kodzie, **marginesem bezpieczeństwa**
(czynnik \(\approx2\) nad pasmem 37,8–71,6, w którym zaobserwowano prawdziwe
awarie numeryczne), nie skalą fizyczną. Główna obserwabla modelu jest więc dla
większości przebiegów wyznaczona przez zapas numeryczny, a nie przez barierę,
którą model deklaruje jako swoją granicę.

**Model wychodzi z wnętrza własnego obszaru niedozwolonego.** Mediana
\(T/t_{\rm light}\) w chwili stopu to 66 przy progu 150, minimum 36 — czyli
przekroczenie o czynnik 2–4, dokładnie w paśmie, w którym udokumentowano
23/23 awarie. Powód jest jeden: przy \(n\lesssim1\) jeden foton niesie około
podwójnego bieżącego wiązania, więc pojedyncza emisja przesuwa periapsis o
czynnik 10–60. Zmierzone wprost wewnątrz jednej kaskady:

```
CASCADE n=2  r_p/r* = 21.71
CASCADE n=3  r_p/r* = 0.347        <- czynnik 62 w JEDNYM fotonie
```

Reguła sprawdzana między checkpointami nie ma jak rozdzielczości poniżej tego
skoku. Sprawdzone osobno i **wykluczone**: żaden foton nie jest emitowany ze
stanu już nieważnego — naruszenie zawsze przypada na ostatni foton kaskady,
łapany na szczycie następnego checkpointu. Przerwanie kaskady nic by nie dało.

**Cały kolaps to 2 albo 3 fotony, i to rozstrzyga wszystko.** Populacje są
rozdzielone bez zakładek:

| | fotony = 2 | fotony = 3 |
|---|---|---|
| wiązanie terminalne | 0,23–0,35 keV | 1,24–3,40 keV |
| \(L/\hbar\) | 0,140–0,173 | 0,045–0,074 |
| \(T/t_{\rm light}\) | 122–150 | **38–70** |
| \(r_p/r^*\) | 9,8–16,5 | 0,31–3,0 |
| przyczyna stopu | retardacja (100%) | mieszane |

Mimośród na końcu to \(0{,}001\)–\(0{,}019\), więc orbita jest kołowa i
\(r_p\propto L^2\); para traci 83–96% momentu pędu. O tym, czy trajektoria ma
dwa fotony czy trzy, decyduje **ostrze noża**: czy po drugim fotonie
\(T/t_{\rm light}\) wypadnie powyżej, czy poniżej 150. Przesuń kwant, a
przesuniesz stan po drugim fotonie i przerzucisz część próbki przez ostrze —
stąd 45%, i stąd głębsze wiązanie terminalne przy mniejszym kwancie.

*Czego NIE zrobiłem i dlaczego.* Napraszało się przycięcie ostatniego fotonu
tak, żeby lądował dokładnie na barierze — mechanizm już w pliku jest
(`roomToFloor` pod `--ground-state-floor`). Odrzucone: model deklaruje, że
poniżej \(r^*\) elektrodynamika punktowa przestaje obowiązywać, więc
dopisanie przejścia kończącego się dokładnie na \(r^*\) **wymyśliłoby**
przejście, którego model nie ma, i ubrałoby granicę ważności w fizykę. Z tego
samego powodu nie przycinam do progu retardacyjnego, który jest wprost
liczbą numeryczną.

*Co zrobiłem.* `stopCause`, `terminalPeriapsisOverBarrier` i
`terminalPeriodToLightCrossing` są teraz zapisywane i raportowane. To zamienia
ukryty artefakt w podaną liczbę: czytelnik widzi, że `terminal binding` nie
jest wiązaniem NA granicy, tylko własnością ziarnistości ostatniego skoku, i
że w 77% przypadków granicą był zapas numeryczny.

*Konsekwencja dla drabiny \(\alpha\).* Sekcja niżej podaje wiązanie na
barierze jako \((2/g)\alpha m_ec^2=3{,}7246\) keV. To jest wartość **przy**
\(r^*\) i pozostaje poprawna jako własność bariery — ale nie jest tym, co
przebiegi raportują jako `terminal binding` (\(\approx1{,}6\)–\(2{,}7\) keV,
zależnie od cenzury), bo większość z nich zatrzymuje się powyżej \(r^*\), na
progu retardacyjnym. Te dwie liczby opisują różne rzeczy i nie należy ich
mylić.

#### Dlaczego żadna zasada zachowania nie wykryła brakującego oddziaływania

Pytanie zasadne: człon siedział w całkowaniu bez odpowiednika w energii, a
wszystko przechodziło 39/39. Przyczyny są trzy i każda jest strukturalna.

**1. Pęd był ślepy z konstrukcji, nie przez niedopatrzenie.**
`chargeDipolePairForces` zwraca `{onCharge, -onCharge}` — reakcja jest
**przypisana**, a nie wyprowadzona z niczego. Zmierzone wcześniej
\(|\mathbf F_1+\mathbf F_2|=0\) co do zera maszynowego na każdym promieniu.
Zasada zachowania, która jest **narzucona**, a nie emergentna, ma zerową moc
diagnostyczną: nie może zawieść, więc nie może o niczym donieść. To był jedyny
kandydat, który wyglądał na kontrolę, a kontrolą nie był.

**2. Energia to widziała, ale próg był 43× za luźny.** Check istnieje
(`rawEnergyResidual < 1.0e-4`) i jest czuły — po dodaniu członu reszta spadła z
\(3{,}97\cdot10^{-6}\) do \(1{,}64\cdot10^{-6}\), czyli brakujący człon
odpowiadał za \(59\%\) jej wartości. Ale \(2{,}33\cdot10^{-6}\) wobec progu
\(10^{-4}\) to **43× poniżej punktu zadziałania**. Sygnał był, tolerancja go
przykryła.

**3. Test, który by to rozstrzygnął, istniał — zastosowany do innej
wielkości.** Przebieg z połową kroku jest liczony od dawna, ale używany
wyłącznie do `stepConvergence`, czyli zbieżności **położeń i prędkości**.
Reszta energii nigdy nie była rafinowana, więc stałe \(4\cdot10^{-6}\)
czytało się jako błąd całkowania. Tymczasem projekt sam formułuje kryterium
przy bilansie długozasięgowym: *fizyczne niezbilansowanie stoi w miejscu, błąd
dyskretyzacji maleje*. Zastosowane tutaj daje odpowiedź natychmiast:

| | reszta (pełny krok) | reszta (pół kroku) | stosunek |
|---|---|---|---|
| bez członu energii | \(3{,}97467\cdot10^{-6}\) | \(3{,}97572\cdot10^{-6}\) | **0,9997** |
| z członem energii | \(1{,}6429\cdot10^{-6}\) | \(1{,}64396\cdot10^{-6}\) | **0,9994** |

Reszta **nie maleje** przy połowieniu kroku. Bez członu było to więc od razu
rozpoznawalne jako **fizyka, nie arytmetyka** — wystarczyło policzyć iloraz z
danych, które już były liczone.

**I to samo kryterium mówi coś nowego** — ale wersja, którą tu najpierw
zapisałem, wymagała trzech poprawek.

*Reszta jest rozłożona na dwa człony, których sumowanie zaciera oba.*
`rawEnergyResidual` dodaje energię **wypromieniowaną** do zmiany energii
mechanicznej. Przy modelu `disabled` żadna siła tego strumienia z mechaniki nie
zabiera, więc dodawanie go nie testuje żadnego bilansu — częściowo znosi dryf
mechaniczny i sprawia, że reszta wygląda na **mniejszą**, niż dryf naprawdę
jest:

```
dE_mech/|E| = -2.4840e-06        <- to jest wielkosc fizyczna
E_rad/|E|   = +8.4112e-07        <- strumien, nieodejmowany
suma        = -1.6429e-06        <- to bylo raportowane
```

Poprawna liczba to więc \(\mathbf{2{,}48\cdot10^{-6}}\), nie
\(1{,}64\cdot10^{-6}\).

*Ale wniosek się broni, i to na czystszej wielkości.* Rafinacja zastosowana do
samego \(dE_{\rm mech}\): pełny krok \(-2{,}48401\cdot10^{-6}\), pół kroku
\(-2{,}48508\cdot10^{-6}\), stosunek \(0{,}9996\) — **płaskie**. Energia
mechaniczna spada przy **wyłączonej** reakcji, gdzie powinna być zachowana
dokładnie, i nie znika to przy zagęszczeniu siatki. Niezbilansowanie jest
realne.

*Wykluczenie Darwina było nieprawidłowe logicznie.* Napisałem, że to „nie jest
próg Darwin-vs-retardacja, bo ten leży dwa rzędy wyżej". Jeśli jakiś próg leży
**wyżej**, to bycie pod nim go nie wyklucza — znaczy tylko, że nie działa. A
sonda **całkuje siły retardowane** (Liénard-Wiechert) przeciwko księdze energii
rzędu **Darwina**, więc ta niezgodność jest tu wręcz kandydatem wiodącym, nie
wykluczonym.

*Próba potwierdzenia jej skalowaniem po \(\beta\) jest zafałszowana przez
konstrukcję sondy.* Zmiana prędkości początkowej zmienia całą trajektorię, a
nie samo \(\beta\); zmierzone \(dE_{\rm mech}\) rośnie, a potem maleje
(\(3{,}4\cdot10^{-7}\to2{,}5\cdot10^{-6}\to2{,}3\cdot10^{-6}\to
1{,}7\cdot10^{-6}\) dla \(\beta_{\rm rel}=0{,}025\)–\(0{,}2\)), więc żadnego
wykładnika z tego nie odczytam.

*I to też było błędne.* Napisałem, że kandydat wiodący jest „niesprawdzony" i
że sonda nie umie go sprawdzić. **Był sprawdzony** — wcześniej w tej samej
sesji, a wynik stoi zapisany w komentarzu przy `conservativeParticleEnergy`, w
pliku, który wtedy edytowałem:

> jedna orbita przy \(a_{Ps}\), reakcja wyłączona; przy dopasowanej sile
> Coulomb+Darwin dryf spada \(5{,}13\cdot10^{-7}\to1{,}28\cdot10^{-7}\to
> \(6{,}42\cdot10^{-8}\to3{,}21\cdot10^{-8}\) przy zaostrzaniu tolerancji
> \(10^{-6}\to10^{-9}\), zbiegając czysto; przy sile retardowanej **stoi na
> \(1{,}3\)–\(2{,}3\cdot10^{-6}\) i nie zbiega wcale**. Zamiecione po
> \(\beta\): wykładnik \(2{,}94\pm0{,}15\), \(\text{floor}/\beta^3\)
> stałe do \(7\%\) na \(3{,}8\)–\(4{,}1\).

Czyli mechanizm jest **zidentyfikowany i skwantyfikowany**: to człon
retardacyjny \(O(\beta^3)\), pierwszy rząd, który Darwin pomija, i rząd, w
którym oddziaływanie dwuciałowe przestaje w ogóle wynikać z potencjału. Nie
brakuje niczego nowego — brakowało mi połączenia reszty z własnym pomiarem
sprzed kilku godzin.

*Przełącznik też istnieje* (`useRetardedExternalForces`), więc „sonda nie umie"
było nieprawdą co do litery. Ale wniosek co do meritum się broni z **innego**
powodu, niż podałem: przełączenie siły zmienia **trajektorię**, więc dwa
przebiegi nie są tą samą orbitą. Zmierzone na tej sondzie, Darwin wypada
**stale o \(\approx4{,}05\cdot10^{-6}\) gorzej**, z członem ładunek-dipol i
bez — czyli odwrotnie, niż dawałoby dopasowanie sił, i jest to sygnatura dwóch
różnych ścieżek, a nie dwóch ksiąg jednej ścieżki. Czysta wersja porównania to
ta zapisana wyżej: jedna orbita przy \(a_{Ps}\), sterowana tolerancją.

*Efekt uboczny, wart zapisania.* Cztery przebiegi dały też niezależne
potwierdzenie samego członu ładunek-dipol:

| | z członem | bez | poprawa |
|---|---|---|---|
| siły retardowane | \(-2{,}484\cdot10^{-6}\) | \(-4{,}809\cdot10^{-6}\) | \(2{,}325\cdot10^{-6}\) |
| siły Darwinowskie | \(-6{,}564\cdot10^{-6}\) | \(-8{,}838\cdot10^{-6}\) | \(2{,}274\cdot10^{-6}\) |

Poprawia **obie** ścieżki o tyle samo, zgodnie do \(2\%\). Jest więc
dopasowany do fizyki, a nie do jednej implementacji — co osobno rozwiewa obawę
zgłoszoną przy jego wdrażaniu, że przypisana z definicji reakcja
(`onDipole = -onCharge`) może go czynić zależnym od ścieżki.

*Wdrożone:* oba człony są teraz raportowane osobno obok reszty, więc jej skład
nie jest już ukryty za sumą.

*Wdrożone:* iloraz rafinacji jest teraz **stale raportowany** obok reszty, z
opisem, jak go czytać. Nie jako bramka — bramka dziś by zawiodła, a wypuszczenie
czerwonego builda po to, żeby zaznaczyć otwarte pytanie, jest gorsze niż
wypisanie liczby, która to pytanie stawia.

#### Dwa udokumentowane braki, sprawdzone: jeden trzyma, w drugim liczba była błędna

Tabela audytu niżej ma dwie pozycje oparte na **uzasadnieniu**, a nie na
pomiarze: M2/E3 i dynamika anihilacji. Obie sprawdzone osobno.

**M2/E3 — struktura trzyma dokładnie, magnitude nie.**

Struktura da się zapisać w postaci zamkniętej. Dla równych mas moment
multipolowy rzędu \(n\) idzie jak \((r/2)^n[q_1+q_2(-1)^n]\), co dla
\(q_2=-q_1\) **znika dla parzystych \(n\)** i przeżywa dla nieparzystych;
szereg magnetyczny dostaje jeszcze jedno odwrócenie znaku z \(\mathbf
v_2=-\mathbf v_1\), więc ma **odwrotną parzystość**:

| | E1 | E2 | E3 | M1 | M2 | M3 |
|---|---|---|---|---|---|---|
| | przeżywa | **znika** | przeżywa | **znika** | przeżywa | znika |

Następny niezerowy rząd po E1 to więc \(\{E3, M2\}\), oba \(\beta^4\) w
mocy — dokładnie jak twierdzi kod.

Ale komentarz dodawał, że \(\beta^4\) „rośnie do rzędu jedności przy samej
barierze Comptona". **Nie rośnie.** Orbita kołowa przy \(r^*\) ma
\(|E|=K/2r^*=(2/g)\alpha m c^2=\tfrac12\mu v^2\), skąd \(\beta^2=4\alpha\):

| miejsce | \(\beta\) | \(\beta^4\) | komentarz |
|---|---|---|---|
| \(a_{Ps}\) | 0,007297 | \(2{,}84\cdot10^{-9}\) | „~3e-9" ✔ |
| próg retardacyjny (150) | 0,041888 | \(3{,}08\cdot10^{-6}\) | — |
| bariera \(r^*\) | 0,170849 | \(8{,}52\cdot10^{-4}\) | „rzędu jedności" ✘ **1174×** |

Błąd jest w stronę **bezpieczną**: zaniedbany rząd jest trzy rzędy mniejszy,
niż głosiło jego własne uzasadnienie, więc decyzja o niewdrażaniu stoi mocniej,
niż sama twierdziła. I znaczy dziś mniej: produkcja kończy na stanie
podstawowym, a nie na barierze, więc punkt końcowy siedzi przy
\(\beta^4=2{,}8\cdot10^{-9}\), nie \(8{,}5\cdot10^{-4}\).

**Dynamika anihilacji — brak trzyma, i da się go podać liczbą.**

Sprawdzone: w kodzie nie ma żadnego rate'u ani czasu życia anihilacji —
jedyne trafienia na „Gamma" to gamma Lorentza i tytuł w bibliografii. Moduł
`qed_reference` niesie **kinematykę** anihilacji (kanał \(s\), Bhabha), nie
jej tempo.

Brak nie jest brakującym współczynnikiem, tylko brakującym pojęciem: tempo
anihilacji jest w QED proporcjonalne do \(|\psi(0)|^2\), a klasyczna orbita o
\(L\neq0\) ma w zerze gęstość dokładnie zero. Konkretnie, jak daleko jest
model od skali, na której anihilacja w ogóle zachodzi:

| skala | wartość | \(a_{Ps}/\text{skala}\) |
|---|---|---|
| \(a_{Ps}\) (postój pod podłogą) | \(105{,}835\) pm | — |
| bariera \(r^*\) | \(193{,}303\) fm | 547 |
| zredukowana dł. Comptona \(\bar\lambda_C\) | \(386{,}159\) fm | **274** |
| klasyczny promień elektronu | \(2{,}818\) fm | 37 558 |

Pod domyślną konfiguracją para parkuje \(274\) zredukowane długości Comptona
od kontaktu.

*Dalszy ciąg tego akapitu był błędny w trzech miejscach* i wymagał osobnego
sprawdzenia; poprawki niżej.

**(a) „W trybie barierowym schodzi do \(r^*=0{,}501\,\bar\lambda_C\)" —
nie.** Zmierzone, periapsis w chwili stopu ma medianę \(1{,}347\,r^*\) i
zakres \(0{,}313\)–\(2{,}275\,r^*\), więc najbliższe podejście to
\(0{,}313\,r^*=\mathbf{0{,}157\,\bar\lambda_C}\) — **wewnątrz** skali
Comptona, nie na niej. Do tego tylko \(42\%\) trajektorii w ogóle kończy na
barierze; \(58\%\) staje wcześniej na progu retardacyjnym. Zacytowałem
nominalną granicę zamiast tego, gdzie trajektorie faktycznie lądują — wbrew
własnemu pomiarowi z sekcji o regule stopu.

**(b) „o \(274\times\) dalszy" — zła arytmetyka.** \(274\) to
\(a_{Ps}/\bar\lambda_C\), a nie stosunek między trybami. Ten wynosi
\(a_{Ps}/r^*=\mathbf{547{,}5}\) wobec nominalnej bariery, a wobec
rzeczywistego najbliższego podejścia (\(0{,}313\,r^*\)) —
\(\mathbf{1749\times}\).

**(c) Liczba, na której to oparłem, była zepsuta przez błąd, który przy tym
sprawdzeniu wyszedł.** `dipoleAwarePeriapsis` zwracał **\(2a\) zamiast
\(a\)** dla orbit kołowych. Dla orbity kołowej apoapsis \(=\) periapsis, więc
ogranicznik `outer` spada na arbitralne \(2\times\) periapsis, a wyjście
„nigdzie nie ma dostępnego obszaru" **zwraca** ten ogranicznik. Zmierzone na
domyślnej konfiguracji, gdzie podłoga parkuje parę w dokładnie kołowym stanie
podstawowym (\(a/a_{Ps}=1\), \(L/L_{\rm circ}=1\), dyskryminant
\(-10^{-14}\)): \(r_p/a\) wychodziło \(2\), gdy człon dipolowy był w
periapsis odpychający, i \(0{,}998\), gdy przyciągający — czyli raportowany
periapsis **przeskakiwał między \(a\) a \(2a\) na znaku członu dipolowego**.
Mediana partii wynosiła \(1095\,r^*\), dokładnie \(2a_{Ps}\), dla pary
siedzącej w \(a_{Ps}\).

Waży to dziś więcej niż kiedyś: ostre przygotowanie i podłoga stanu
podstawowego robią orbity dokładnie kołowe **z konstrukcji**.

*Naprawione* w dwóch miejscach: sprawdzenie zdegenerowanego podwójnego
pierwiastka, oraz fallback zwracający periapsis keplerowski zamiast
arbitralnego ogranicznika, gdy `outer` nie jest fizyczny. Po poprawce tryb
domyślny raportuje \(547{,}509\,r^*\) (zakres \(546{,}5\)–\(547{,}5\)),
czyli \(a_{Ps}\), a górna krawędź zakresu barierowego spada z \(4{,}551\) na
\(2{,}275\,r^*\) — też była artefaktem \(2\times\). Mediana kolapsu, średnia i
podział przyczyn stopu **bez zmian w obu trybach**: poprawka dotyczy raportu,
nie dynamiki.

*Co z tego zostaje w mocy:* brak anihilacji jest strukturalny w obu trybach.
Pod domyślną konfiguracją para stoi \(274\,\bar\lambda_C\) od kontaktu; w
trybie barierowym dochodzi do \(0{,}157\,\bar\lambda_C\) — czyli **wchodzi**
w skalę Comptona, ale robi to w obszarze, który model deklaruje jako
nieważny, więc nic z tego nie wynika poza tym, że tam się zatrzymuje.

(Przy okazji test spójności: \(r^*=g\hbar/4m_ec=(g/4)\bar\lambda_C\), a
zmierzone \(r^*/\bar\lambda_C=0{,}501\) — bariera to dokładnie pół
zredukowanej długości Comptona, z dokładnością do \(g/2\).)

#### Audyt kompletności: który człon jest w której warstwie

Pytanie: czy wszystkie policzone efekty — spin-orbita, orbita-spin i tak dalej
— mają odzwierciedlenie w CREM. Model ma **trzy warstwy**, w których człon może
być obecny lub nie, i to nie to samo: siła w całkowaniu mechanicznym, energia w
księdze zachowawczej `conservativeParticleEnergy`, oraz reprezentacja w
warstwie sekularnej (elementy \((E,L)\), `dipoleAwarePeriapsis`, hazard
emisji).

| człon | siła | energia | warstwa sekularna |
|---|---|---|---|
| ładunek–ładunek (kulombowski) | ✔ | ✔ | ✔ |
| ładunek–ładunek \(v^2/c^2\) (Darwin) | ✔ | ✔ | — (przez zmierzoną orbitę) |
| dipol–dipol (spin–spin) | ✔ | ✔ | ✔ `azimuthAveragedDipoleEnergy` |
| ładunek–dipol (spin–orbita, oba parowania) | ✔ | **BRAK** | ✔ tylko jako moment siły |
| precesja BMT + Thomas | ✔ (moment) | — (nie wykonuje pracy) | — |
| moment reakcji promienistej na dipole | ✔ | ✔ | — |
| promieniowanie E1, retardowane | ✔ | ✔ | ✔ |
| promieniowanie M1 (koherentne) | ✔ | ✔ | ✔ |
| M2/E3 | — | — | — (udokumentowane jako \(\beta^4\)) |
| dynamika anihilacji | — | — | — (udokumentowany brak) |

Odpowiedź brzmi więc: **prawie wszystkie, z jednym wyjątkiem, i ten wyjątek
nie jest mały.**

*Znaleziona luka: energia oddziaływania ładunek–dipol nie istnieje nigdzie w
modelu.* Siła jest (`chargeDipoleForces`, oba parowania, domknięta para
akcja–reakcja), ale `conservativeParticleEnergy` niesie tylko kinetyczną,
kulombowską, dipol–dipol, Darwina i więzy. Nie ma funkcji
`chargeDipoleInteractionEnergy` — sprawdzone przez wszystkie moduły.

*Rozmiar.* Z \(B=\frac{\mu_0}{4\pi}\frac{qv}{r^2}\) i
\(\mu=\frac{e\hbar}{2m}\) wychodzi zamknięcie
\(U_{\rm so}/U_{\rm dd}=4L/\hbar\), potwierdzone numerycznie:

| \(r/r^*\) | \(L/\hbar\) | \(|U_{dd}|/|U_C|\) | \(|U_{\rm so}|/|U_C|\) |
|---|---|---|---|
| 1 (bariera) | 0,043 | 2,0 | **0,342** |
| 1,5 | 0,052 | 0,889 | 0,186 |
| 3 | 0,074 | 0,222 | 0,066 |
| \(a_{Ps}\) | 1,0 | \(6{,}7\cdot10^{-6}\) | \(2{,}7\cdot10^{-5}\) |

(orientacje maksymalizujące, więc górne ograniczenia). Człon nieobecny w
księdze jest wart do \(34\%\) członu kulombowskiego przy barierze — a człon
dipol–dipol, który **jest** w księdze, w tym samym miejscu.

*Struktura — z poprawką, bo pierwsze wyprowadzenie było błędne.* Napisałem tu,
że „ładunki i prędkości obu cząstek są przeciwne, więc oba pola dodają się",
i wyszło mi \(U_{\rm so}=-(\boldsymbol\mu_1+\boldsymbol\mu_2)\cdot\mathbf
B\), czyli struktura M1: maksimum dla para, zero dla ortho. **To jest
odwrócone.** Policzyłem trzy odwrócenia znaku jako dwa: odwraca się ładunek,
prędkość **oraz** kierunek do punktu pola. Rachunek do końca:

\[U_1=-\frac{\mu_0}{4\pi}\frac{q_2\,\boldsymbol\mu_1\cdot(\mathbf
v_2\times\hat r_{12})}{r^2},\qquad
U_2=+\frac{\mu_0}{4\pi}\frac{q_1\,\boldsymbol\mu_2\cdot(\mathbf
v_1\times\hat r_{12})}{r^2},\]

a po podstawieniu \(q_2=-q_1\), \(\mathbf v_2=-\mathbf v_1\):

\[U_1+U_2=\frac{\mu_0}{4\pi}\frac{q_1\,(\boldsymbol\mu_2-\boldsymbol
\mu_1)\cdot(\mathbf v_1\times\hat r_{12})}{r^2}.\]

**Różnica momentów, nie suma** — czyli kombinacja **antykoherentna**: zerowa
dla para (równoległe) i maksymalna dla ortho (antyrównoległe). Dokładne
dopełnienie M1, a nie jego powtórzenie.

*Zmierzone, i to właśnie ten pomiar wykrył błąd.* Reszta bilansu energii
(`reaction raw off`, czyli zachowanie energii przy wyłączonej reakcji
promienistej), z członem i bez:

| konfiguracja | z członem | bez członu | różnica |
|---|---|---|---|
| ortho (antyrównoległe) | \(3{,}9746716\cdot10^{-6}\) | \(3{,}9649546\cdot10^{-6}\) | \(+9{,}7\cdot10^{-9}\) |
| para (równoległe) | \(4{,}0128854\cdot10^{-6}\) | \(4{,}0132439\cdot10^{-6}\) | \(-3{,}6\cdot10^{-10}\) |

Efekt jest **27× większy dla ortho** — dokładnie ta struktura, którą daje
\((\boldsymbol\mu_2-\boldsymbol\mu_1)\), i przeciwna do tej, którą
zapowiadałem.

*Dwie rzeczy, które czyniły ten pomiar pustym, dopóki ich nie usunąłem.*
Sonda walidacyjna ma momenty ustawione **antyrównolegle**, więc dla mojej
(błędnej) hipotezy była ślepa z konstrukcji; i sprzężenie ładunek–dipol wchodzi
tam **ścieżką retardowaną**, jako pole dipolowe partnera wewnątrz siły
Lorentza, a nie przez `chargeDipoleForces`. Zabramkowanie samej
`chargeDipoleForces` dawało wynik identyczny co do cyfry — trzeba było
zabramkować `retardedMagneticDipoleField`.

*Skala.* Człon zaburza resztę bilansu o \(\approx0{,}24\%\) jej wartości
przy promieniu sondy (\(\approx6000\,r^*\)), gdzie sam jest wart
\(\approx7\cdot10^{-7}\) członu kulombowskiego. Czyli **około 1% jego
wielkości ujawnia się jako niezbilansowanie** — reszta pracy nie wykonuje, co
jest oczekiwane dla siły magnetycznej. Ekstrapolacja do bariery, gdzie człon
sięga \(0{,}34\) kulombowskiego, dawałaby \(\sim3\cdot10^{-3}\) — ale to
jest ekstrapolacja, nie pomiar: sonda nie sięga tego promienia.

Model ma więc dwie dyskryminanty kanału w dwóch miejscach:
\(\boldsymbol\mu_1\cdot\boldsymbol\mu_2\) w sekularnym momencie siły i
\((\boldsymbol\mu_2-\boldsymbol\mu_1)\) w energii spin-orbita.

*Wdrożone — i zweryfikowane pomiarem, który mógł je obalić.*
`chargeDipoleInteractionEnergy` dokłada ten człon do
`conservativeParticleEnergy`, z tą samą wagą regularizacyjną, której używa pole
dipolowe, więc przy małych \(r\) tłumi go ten sam regulator, a nie drugi,
niezależny.

Ryzyko było realne: siła jest zaimplementowana jako \(q(\mathbf
v\times\mathbf B)\) z reakcją przypisaną **z definicji**
(`onDipole = -onCharge`), więc nie musi wynikać z żadnego potencjału. Gdyby
dołożona energia nie odpowiadała tej sile, reszta bilansu by **wzrosła**. Test
rozstrzygnął jednoznacznie:

| | bez członu | z członem |
|---|---|---|
| `reaction raw off` | \(3{,}9746716\cdot10^{-6}\) | \(\mathbf{1{,}642896\cdot10^{-6}}\) (\(-59\%\)) |
| `reaction raw LL/C` | \(4{,}4028823\cdot10^{-6}\) | \(\mathbf{2{,}0764777\cdot10^{-6}}\) (\(-53\%\)) |

Reszta zachowania energii **spadła o ponad połowę**. To brakująca energia, nie
sztuczny dodatek.

*Produkcja nietknięta:* mediana \(7368{,}8\) ps identycznie przed i po, w obu
kanałach; w trybie barierowym również bez zmian. Poprawka poprawnościowa, nie
zmiana fizyki obserwabli.

*Druga warstwa, uczciwie.* Audyt znalazł ten człon brakujący **w dwóch**
warstwach, więc dołożyłem go też do `dipoleAwarePeriapsis` — w punkcie
zwrotnym forma jest ścisła, bo tam \(v_r=0\) i \(v=L/\mu r\) stycznie, co
czyni go funkcją samego \(r\) i \(L\). Zmierzony na rzeczywistej trajektorii:
\(|U_{\rm so}|/|U_C|=2{,}52\cdot10^{-6}\) przy \(a_{Ps}\), czyli **większy
niż stojący obok człon dipol–dipol** (\(1{,}37\cdot10^{-6}\)); skaluje się jak
\(L/r^2\), a że \(L\) spada z \(1\hbar\) do \(\approx0{,}045\hbar\), przy
barierze sięga kilku procent członu kulombowskiego — nie \(34\%\), które
podawała tabela audytowa przy orientacji maksymalizującej i pełnym
\(2\mu_B\), podczas gdy zmierzone \(|(\boldsymbol\mu_2-\boldsymbol\mu_1)
\cdot\hat n|\) to \(0{,}38\,\mu_B\).

*Powyższe było błędne w obie strony.* Napisałem, że dla warstwy sekularnej
„nie ma takiego testu" i że „nie zmieniła niczego". Test **jest**, i to
oczywisty — policzyć formę sekularną i `chargeDipoleInteractionEnergy` na
**tym samym** stanie w punkcie zwrotnym i zażądać zgodności. Po prostu na to
nie wpadłem.

*Uruchomiony, znalazł dwa błędy* w ręcznie wyprowadzonej formie zamkniętej,
która tu stała:

| \(r/r^*\) | \(U\) z księgi [J] | \(U\) sekularne [J] | stosunek |
|---|---|---|---|
| 1 | \(-1{,}01895\cdot10^{-16}\) | \(+2{,}03781\cdot10^{-16}\) | \(-1{,}99992\) |
| 10 | \(-6{,}38061\cdot10^{-19}\) | \(+6{,}44441\cdot10^{-19}\) | \(-1{,}01\) |
| 547,5 | \(-2{,}90539\cdot10^{-23}\) | \(+2{,}90540\cdot10^{-23}\) | \(-1{,}00000\) |

Zły **znak** — bo \(\mathbf v\times\hat r=-v\hat L\), a nie \(+v\hat L\),
więc człon wchodził z odwrotnym znakiem na **każdej** trajektorii — oraz
niezgodność regularizacji krótkiego zasięgu, narastająca do **czynnika 2 przy
barierze**, bo księga obcina wektor separacji do podłogi, a forma zamknięta
nie.

*Naprawa.* Warstwa sekularna **woła** teraz `chargeDipoleInteractionEnergy` na
zbudowanym stanie w punkcie zwrotnym, zamiast powtarzać jej wzór. Tożsamość
zachodzi wtedy z konstrukcji. To ta sama lekcja, co przy `quantumFor`: jedna
recepta, w jednym miejscu — drugą kopię zawsze trzeba potem gonić.

*A wynik zerowy, który podałem jako uspokojenie, nie niósł żadnej informacji.*
Człon, którego czynnik orientacyjny ma losowy znak, a znak globalny jest
błędny, jest w agregacie **statystycznie nieodróżnialny** od poprawnego. Brak
efektu świadczył o tym, że test nie miał mocy, a nie o tym, że człon jest
nieszkodliwy.

*Zmierzone porządnie,* z poprawionym członem i przy 24/24 ukończonych
trajektoriach, więc porównanie jest sparowane, a nie ograniczone cenzurą:
mediana \(208{,}007\) ps i średnia \(223{,}308\pm37{,}95\) ps **identyczne**
z członem i bez, wiązanie terminalne identyczne, a **dokładnie jedna
trajektoria na 24** przechodzi z progu retardacyjnego na barierę Comptona.
Człon pozostaje więc niemal bezczynny — ale to jest teraz pomiar członu
**poprawnego**, a nie zepsutego.

(Po drodze odczyt pośredni \(238\to208\) ps wyglądał na efekt rzeczywisty;
okazał się artefaktem cenzury — 13 wobec 14 ukończeń. Przebieg z pełnym
ukończeniem go usunął.)

#### Orbitalny moment magnetyczny i sprzężenie spin-orbita

Pytanie: para i ortho różnią się nie tylko orientacją momentów, ale i tym, jak
te momenty wpływają na **orbitalny moment magnetyczny** pozytonium — czy jest
to uwzględnione?

*Pierwsza połowa odpowiedzi: nie ma na co wpływać.* W układzie CM

\[\boldsymbol\mu_L=\tfrac12\sum_i q_i\,\mathbf r_i\times\mathbf v_i
  =\frac{q_1m_2^2+q_2m_1^2}{2M^2}\cdot\frac{\mathbf L}{\mu_{\rm red}},\]

a dla pozytonium \(q_1=-e,\ q_2=+e,\ m_1=m_2\), więc licznik to
\(-em^2+em^2=0\). **Orbitalny moment magnetyczny pozytonium znika dokładnie.**
Nie jest to nowe ustalenie — `electrodynamics.hpp` już to notuje i sprawdził
numerycznie, przy okazji audytu kompletności multipolowej (to samo zerowanie
usuwa M1 i E2 z listy następnych rzędów, zostawiając M2/E3).

*Druga połowa: właściwy kanał istnieje i JEST w modelu.* Sprzężenie
spin-orbita nie potrzebuje orbitalnego momentu — działa przez pole magnetyczne
spinu partnera na **poruszający się ładunek**, \(q\,\mathbf v\times\mathbf
B_{\rm spin}\), plus reakcja. To `chargeDipoleForces`, zbudowane na pochodnych
Thomasa-BMT i wchodzące do `allExternalForces` oraz do gałęzi retardowanej.
Zweryfikowane jako domknięta para akcja-reakcja: \(|\mathbf F_1+\mathbf
F_2|=0\) dokładnie, na każdym zmierzonym promieniu.

*Rozmiar, i pozorna asymetria kanałów.* W geometrii z momentami wzdłuż
\(\mathbf L\) siła ta jest **czysto ortho**:

| \(r/r^*\) | para | ortho |
|---|---|---|
| 1 | **0** | 0,2415 |
| 2 | **0** | 0,1080 |
| 20 | **0** | 0,0038 |

(w jednostkach \(|F_{\rm Coulomb}|\); w geometrii z momentami wzdłuż separacji
ortho sięga 0,479). Wygląda to jak dopełnienie struktury M1 — M1 wzmacnia się
dla para i znika dla ortho, tu odwrotnie — czyli klasyczny odpowiednik reguły
C-parzystości.

*Ale znika po uśrednieniu.* Wielkością, która mogłaby ruszyć periapsis, jest
moment siły względem \(\mathbf L\). Uśredniony po izotropowych orientacjach i
fazach orbity, 20 000 próbek:

| \(r/r^*\) | kanał | \(\langle\tau_L\rangle\) | \(\langle|\tau_L|\rangle\) | stosunek |
|---|---|---|---|---|
| 1 | para | \(-1{,}43\cdot10^{-5}\) | \(6{,}57\cdot10^{-3}\) | \(-0{,}0022\) |
| 1 | ortho | \(-4{,}34\cdot10^{-5}\) | \(6{,}47\cdot10^{-3}\) | \(-0{,}0067\) |
| 2 | para | \(-2{,}76\cdot10^{-6}\) | \(1{,}665\cdot10^{-3}\) | \(-0{,}0017\) |
| 2 | ortho | \(-3{,}57\cdot10^{-6}\) | \(1{,}645\cdot10^{-3}\) | \(-0{,}0022\) |

Średnia to \(0{,}2\)–\(0{,}7\%\) wartości bezwzględnej przy szumie
\(1/\sqrt{20000}=0{,}7\%\), czyli **zero**. A \(\langle|\tau_L|\rangle\)
różni się między kanałami o \(1{,}5\%\) — **rozrzut jest ten sam**. Widoczna
asymetria była więc artefaktem szczególnej geometrii.

*Uwaga wycofana.* Napisałem tu wcześniej, że jest to ten sam washout co przy
M1, ale z różnicą: „kwantowanie spinu tu nie pomoże, bo uśrednienie idzie po
wspólnym kierunku momentów, którego kwantowanie nie dotyka — pomiar miał
\(\cos=\pm1\) ustalone i mimo to dał zero". **Obie części tego są błędne.**

Rozumowanie jest odwrócone. Wspólny kierunek wymywa się **całkowicie** — to
jest właśnie owo stłumienie \(1{,}6\cdot10^6\) — a to, co zostaje, zależy
**wyłącznie** od kąta wzajemnego, czyli dokładnie od tego, co kwantowanie
ustala. Powołany „pomiar" miał wprawdzie \(\cos=\pm1\), ale nie dał zera:
dał \(3{,}7\cdot10^{-9}\), które ukrył próg szumu Monte Carlo
(\(4{,}6\cdot10^{-5}\)).

*Kwantowanie pomaga, i najbardziej ortho.* Zmierzone na 40 trajektoriach na
kanał:

| | \(\langle\cos\rangle\) bez podłogi | zakres | z kwantowaniem | wzmocnienie |
|---|---|---|---|---|
| para | \(+0{,}713\) | \([0{,}516,\ 0{,}984]\) | \(+1\) dokładnie | \(\times1{,}40\) |
| ortho | \(-0{,}271\) | \([-0{,}991,\ +0{,}477]\) | \(-1\) dokładnie | \(\times3{,}69\) |

(analitycznie \(+0{,}75\) i \(-0{,}25\), bo \(\cos\) jest jednostajny na
\([-1,1]\), a cięcie idzie po \(0{,}5\)).

*Wniosek dla samej flagi i tak się broni — z zupełnie innego powodu.* Podłoga
zatrzymuje parę przy \(a_{Ps}\approx550\,r^*\), a człon spada stromo z
promieniem, więc para nigdy nie dociera tam, gdzie jest on duży. Zmierzone
skumulowane \(|dL|/L\): \(2\cdot10^{-14}\) z podłogą wobec
\(6{,}4\cdot10^{-3}\) bez niej. **Podłoga usuwa obszar, nie sprzężenie.**

*Jak czytać znak.* Na pojedynczej trajektorii znak skumulowanego \(dL\) idzie
za \(\text{sign}(\cos)\), jak proporcjonalność wymaga: 37/40 dla para, 31/40
dla ortho — a niezgodności to przebiegi, na których precesja BMT odsunęła
\(\cos\) od wartości początkowej. Znak sumy **zespołowej** nie jest tego
wiarygodnym odczytem, z tego samego powodu, co przy sumie energii fotonów:
dominuje ją jedna głęboka trajektoria. W jednej partii ortho ta dominująca
niosła \(dL=+0{,}205\) przy \(\cos=+0{,}414\) — dodatnim, co zakres ortho
\([-1,\,0{,}5)\) dopuszcza — przez co cała suma wyszła dodatnia, choć każdy
składnik z osobna szedł za własnym \(\cos\).

*Poprawka do powyższego, i wdrożenie.* Napisałem wyżej, że uśrednienie daje
zero — to było zero **w granicach rozdzielczości tamtej sondy**, która losowała
fazę orbity metodą Monte Carlo i miała próg szumu \(4{,}6\cdot10^{-5}\).
Pod właściwą kwadraturą po fazie reszta się rozdziela i **nie jest zerem**:

| \(r/r^*\) | para | ortho |
|---|---|---|
| 1 | \(+3{,}86\cdot10^{-9}\) | \(-3{,}69\cdot10^{-9}\) |
| 2 | \(+3{,}35\cdot10^{-10}\) | \(-3{,}30\cdot10^{-10}\) |
| 20 | \(+1{,}49\cdot10^{-14}\) | \(-1{,}47\cdot10^{-14}\) |

Uśrednienie po fazie tłumi moment siły \(1{,}6\) **miliona** razy — prawie
wszystko to drganie wewnątrzorbitalne — ale to, co zostaje, jest
**systematyczne**, nie rozrzutem.

*Doprecyzowane na deterministycznej siatce.* Powyższe liczby pochodziły z
losowania Monte Carlo po orientacjach i były antysymetryczne tylko z grubsza.
Na siatce iloczynowej \(48\times96\times64\) są **dokładne**:

```
<tau_L>    para = +3.69299e-09   ortho = -3.69299e-09   suma 8.6e-21
<|tau_L|>  para =  ortho co do ostatniego bitu, roznica wzgledna 0
```

Czyli oba twierdzenia z pierwszego podejścia wymagały poprawki, w przeciwne
strony. **Rozrzut jest identyczny dokładnie**, nie „w granicach 1,5%" — tamte
\(1{,}5\%\) były szumem próbkowania (błąd standardowy \(\approx0{,}7\%\),
więc różnica siedziała na dwóch sigmach i nie wolno było nazwać jej zerem).
Ale **antysymetria też jest dokładna**, a nie przybliżona — i to znaczy, że
teza „asymetria była artefaktem narzuconej geometrii" **była błędna**.

*Dlaczego.* Przemiatanie kąta wzajemnego pokazuje, że to, co przeżywa
uśrednienie, jest wprost proporcjonalne do \(\boldsymbol\mu_1\cdot
\boldsymbol\mu_2\):

| \(\cos\) | 1 | 0,75 | 0,5 | 0,25 | 0 | \(-0{,}5\) | \(-1\) |
|---|---|---|---|---|---|---|---|
| \(\langle\tau_L\rangle/\cos\) | \(3{,}69539\cdot10^{-9}\) | \(3{,}69539\) | \(3{,}69539\) | \(3{,}69539\) | — | \(3{,}69539\) | \(3{,}69539\) |

Do sześciu cyfr przy każdym niezerowym \(\cos\), a przy \(\cos=0\) średnia
wynosi \(-9{,}2\cdot10^{-21}\), czyli zero maszynowe. Wymywa się część
**liniowa** w momentach; część **dwuliniowa** \((\boldsymbol\mu_1\cdot
\boldsymbol\mu_2)\) nie wymywa się wcale. To, co zostaje, **jest**
dyskryminantą kanału.

*Konsekwencja produkcyjna, wcześniej przeoczona.* Kanały próbkuje się z
**zakresów** \(\cos\), nie z \(\pm1\): kierunki są izotropowe, więc
\(\cos\) jest jednostajny na \([-1,1]\), para bierze \(\cos\ge0{,}5\), a
ortho resztę. Stąd \(\langle\cos\rangle=+0{,}75\) dla para i
\(\mathbf{-0{,}25}\) dla ortho — człon jest **trzykrotnie silniejszy dla
para** i przeciwnego znaku. Pod kwantowaniem spinu
(`--ground-state-floor`) \(\cos\) staje się dokładnie \(\pm1\) i oba kanały
robią się symetryczne.

*Spadek z promieniem* nie jest pojedynczą potęgą: nachylenie lokalne wynosi
\(3{,}5\) między \(r^*\) a \(2r^*\) i \(4{,}35\) między \(2r^*\) a
\(20r^*\). Tak czy inaczej jest to efekt obszaru terminalnego — dokładnie
tam, gdzie reguła stopu o wszystkim decyduje.

Rozróżnienie, które to zmienia: zerowa średnia **po zespole** nie znaczy zerowy
efekt **na trajektorię**. Kierunek momentów jest wzdłuż jednego przebiegu
ustalony, więc każda trajektoria niesie własny systematyczny dryf \(L\),
znoszący się dopiero w średniej po wielu przebiegach.

*Wdrożone.* `azimuthAveragedSpinOrbitTorque` liczy ten moment siły
kwadraturą 64-węzłową po fazie orbity kołowej przy bieżącej półosi (mimośród
terminalny to \(0{,}001\)–\(0{,}019\), więc kołowość jest dobra znacznie
poniżej dokładności samego członu) i jest stosowany do \(L\) raz na
checkpoint, przez czas trwania przeskoku. Energii **nie** rusza: część
magnetyczna nie wykonuje pracy, a jedna mechanicznie całkowana orbita już
całkuje pełną siłę, więc dodanie tego do energii byłoby podwójnym liczeniem.

*Zbieżność kwadratury* sprawdzona, nie założona: skumulowane \(|dL|/L\) na
trajektorię wychodzi identyczne do sześciu cyfr przy 16, 32, 64, 256 i 1024
węzłach.

*Rozmiar efektu* — i to jest powód, dla którego członu nie wolno odrzucić po
medianie. Skumulowane \(|dL|/L\) na trajektorię:

| kanał | min | mediana | max |
|---|---|---|---|
| para | \(1{,}1\cdot10^{-8}\) | \(1{,}8\cdot10^{-5}\) | **0,255** |
| ortho | \(2{,}4\cdot10^{-9}\) | \(2{,}3\cdot10^{-5}\) | **0,205** |

Dla większości trajektorii to nic. Dla tych, które schodzą najgłębiej — czyli
tych przy ostrzu noża — to **ćwierć momentu pędu**. Rozkład rozciąga się na
cztery rzędy wielkości, co jest \(r^{-4{,}2}\) w działaniu.

*Wpływ na produkcję,* A/B wobec HEAD na tych samych ziarnach: mediana czasu
kolapsu **niezmieniona** (175,202 ps w obu), wiązanie terminalne niezmienione
do pięciu cyfr, a przyczyna stopu przenosi się dla 1 trajektorii na 27 z
bariery na próg retardacyjny. Liczba ukończeń pod stałym budżetem wallclock nie
drgnęła, więc 64 wywołania siły na checkpoint są przy mechanicznie całkowanej
orbicie niezauważalne.

#### Ostrze noża wobec obu członów magnetycznych: nadal nie, ale kanały jednak się różnią

Przeliczone po raz trzeci, bo od poprzedniego sprawdzenia energia
**ładunek-dipol** weszła wprost do `dipoleAwarePeriapsis` — czyli do wielkości
decydującej o ostrzu — i idzie jak \((\boldsymbol\mu_2-\boldsymbol\mu_1)\),
więc jest zerowa dla para i maksymalna dla ortho. Do tego naprawiono w tej
funkcji błąd znaku i błąd orbity kołowej. Tym razem oba człony bramkowane
osobno, sparowane para/ortho, 24 trajektorie, wszystkie ukończone:

| konfiguracja | flipy liczby fotonów | różnice przyczyny stopu |
|---|---|---|
| oba człony | **0** | 16 / 24 |
| tylko dipol–dipol | **0** | 16 / 24 |
| tylko ładunek–dipol | **0** | 4 / 24 |
| żaden | **0** | **0 / 24** |

**Odpowiedź na pytanie: nie.** Liczba fotonów nie zmienia się między kanałami w
żadnej z czterech konfiguracji. Ostrze \(2\) wobec \(3\) fotonów jest
nietknięte, więc \(45\%\) nie da się na te człony zrzucić.

**Ale wiersz ostatni mówi coś, czego wcześniej nie widziałem.** Z **żadnym** z
tych członów para i ortho dają **identyczne** przyczyny stopu — zero różnic na
24 trajektoriach. Cała zależność kanałowa reguły stopu jest więc niesiona
dokładnie przez te dwa człony sekularne: dipol–dipol dominująco (16/24),
ładunek–dipol słabiej (4/24). To nie jest nic, tylko trafia w atrybucję
granicy, a nie w obserwablę.

*I nie dociera do niczego mierzalnego.* Czasy życia zgadzają się między
kanałami **co do szóstej cyfry nawet tam, gdzie przyczyna stopu jest inna**
(\(398{,}513\) wobec \(398{,}513\); \(209{,}165\) wobec \(209{,}164\)),
mediany są identyczne (\(208{,}007\) ps), a wiązanie terminalne różni się o
\(0{,}03\%\) (\(2{,}5715\) wobec \(2{,}5722\) keV).

*Pułapka, w którą prawie wpadłem.* Na ziarnie 42 udział stopów na barierze
wyszedł \(41{,}7\%\) dla para wobec \(66{,}7\%\) dla ortho — 25 punktów
procentowych, co wygląda na duży efekt systematyczny. Na czterech ziarnach:

| ziarno | para | ortho | różnica |
|---|---|---|---|
| 42 | 41,7% | 66,7% | \(+25{,}0\) |
| 7 | 50,0% | 37,5% | \(\mathbf{-12{,}5}\) |
| 99 | 37,5% | 41,7% | \(+4{,}2\) |
| 1234 | 33,3% | 54,2% | \(+20{,}8\) |

Średnia \(+9{,}4\pm8{,}6\) punktu, czyli \(1{,}1\sigma\), znak zmienny.
Zgodne z zerem. Pojedyncze ziarno dałoby tu wniosek fałszywy — ten sam wzorzec,
co przy sumie energii fotonów i przy uśrednianiu momentu spin-orbita.

*Zestawienie.* Człony magnetyczne **rozdzielają kanały w regule stopu**
trajektoria po trajektorii (16 z 24 zmienia atrybucję), ale ani nie przerzucają
ostrza, ani nie zmieniają żadnej wielkości raportowanej, ani nie dają
systematycznego przesunięcia udziałów. Ostrze zostaje własnością ziarnistości
emisji.

#### Ostrze noża po wprowadzeniu członu spin-orbita: nadal nie, ale teraz wiadomo dlaczego

Poprzednia sekcja wykluczyła kanał rozpadu jako sprawcę ostrza, ale mierzyła
kod **sprzed** wprowadzenia sekularnego członu spin-orbita. Ten człon jest
dokładnie proporcjonalny do \(\boldsymbol\mu_1\cdot\boldsymbol\mu_2\),
trzykrotnie silniejszy dla para, przeciwnego znaku i sięga 25% momentu pędu na
najgłębszych trajektoriach — więc pytanie wymagało przeliczenia.

Sparowane para/ortho, 38 par ukończonych w obu konfiguracjach, ten sam
początkowy \((E,L)\):

| | zmiana liczby fotonów | \(|\Delta\text{PLC}/\text{PLC}|\) średnia | max |
|---|---|---|---|
| przed członem | **0** | \(3{,}16\cdot10^{-4}\) | \(2{,}33\cdot10^{-3}\) |
| po członie | **0** | \(5{,}49\cdot10^{-3}\) | \(8{,}52\cdot10^{-2}\) |

Człon wzmocnił różnicę kanałów **17× w średniej i 37× w maksimum**, do
\(8{,}5\%\) — a mimo to nie przerzucił ani jednej trajektorii. Przy
\(8{,}5\%\) pasmo flipu byłoby szerokie na ~13 jednostek PLC wobec zakresu
\(122\)–\(150\) populacji dwufotonowej, więc flipy **powinny** były wystąpić.
Nie wystąpiły, i rozbicie na populacje mówi dlaczego:

| populacja | PLC przy stopie | \(|\Delta\text{PLC}/\text{PLC}|\) |
|---|---|---|
| 2 fotony (**przy ostrzu**) | \(121{,}9\)–\(149{,}7\) | \(0\) dla 9 z 11, \(7\cdot10^{-6}\) dla reszty |
| 3 fotony (głęboko) | \(32{,}4\)–\(70{,}2\) | do \(8{,}5\cdot10^{-2}\) |

**Człon jest cztery rzędy wielkości większy tam, gdzie nie może nic
przerzucić.** Spadek jak \(r^{-4}\) działa przeciwko niemu: żyje przy małych
promieniach, czyli na trajektoriach trzyfotonowych, które siedzą przy
\(\text{PLC}=32\)–\(70\), daleko **pod** progiem 150 i już nic ich nie
uratuje. Populacja dwufotonowa, która o ostrzu decyduje, zatrzymuje się przy
\(10\)–\(16\,r^*\), gdzie człon jest praktycznie zerowy.

*Poprawka do wcześniejszego oszacowania.* W poprzedniej sekcji podałem
częstość flipu rzędu 1 na 250, biorąc maksimum \(|\Delta\text{PLC}|\) po
całej próbce. To było błędne o dwa rzędy: maksimum pochodzi z trajektorii
głębokich, a przy samym ostrzu przesunięcie wynosi \(\le7\cdot10^{-6}\), co
daje pasmo szerokości \(0{,}001\) w zakresie 28 jednostek, czyli **1 na
27 000**.

*Wniosek.* Kanał rozpadu wchodzi teraz do dynamiki orbitalnej realnie i
dokładnie proporcjonalnie do \(\boldsymbol\mu_1\cdot\boldsymbol\mu_2\), ale
**nie może** przerzucić ostrza — nie dlatego, że jest za mały, tylko dlatego,
że jest duży w niewłaściwym miejscu. Ostrze pozostaje własnością ziarnistości
emisji.

#### Czy ostrze noża jest efektem dipol-dipol? Nie — zmierzone i wykluczone

Podejrzenie było trafne co do miejsca: wniosek „para i ortho nie różnią się"
ustalono dla **budżetu energii**, a ostrze noża siedzi w **regule stopu**,
gdzie człon dipolowy wchodzi przez `dipoleAwarePeriapsis` i przy promieniu
terminalnym jest wart dziesiątki procent członu kulombowskiego, nie
\(10^{-5}\).

Konfiguracja nadaje się do tego idealnie: prędkości losują się **przed**
dipolami, więc para i ortho przy tym samym ziarnie mają identyczny początkowy
\((E,L)\) i różni je wyłącznie kąt wzajemny momentów. Sparowane, 40 przebiegów
na ziarno:

| ziarno | para: 2/3 fotony | ortho: 2/3 fotony | mediana PLC para | ortho |
|---|---|---|---|---|
| 42 | 11 / 25 | 11 / 25 | 66,1931 | 66,1933 |
| 7 | 7 / 30 | 7 / 30 | 48,6899 | 48,6841 |
| 99 | 8 / 26 | 8 / 26 | 57,5085 | 57,5085 |

Trajektoria po trajektorii, na 36 par ukończonych w obu konfiguracjach:

```
zmiana liczby fotonow:      0
|dPLC/PLC|:  srednia 3.3e-04,  max 2.3e-03
```

*Kontrola przeszła* — porównanie nie jest puste: człon dipolowy ma przeciwne
znaki, \(+0{,}310\) keV dla para i \(-0{,}337\) keV dla ortho, przy wiązaniu
terminalnym \(1{,}45\) keV, czyli **21%**. Konfiguracje różnią się mocno, a
ostrze i tak się nie rusza.

*Dlaczego tak mało.* To są **dwie różne wielkości**, i mylenie ich było moim
własnym błędem w postawieniu hipotezy. Te \(\pm0{,}31\) keV to wartość
**chwilowa**, dodawana do niezmiennika anihilacji \(W\). Periapsis czyta
wartość **uśrednioną azymutalnie** — tę, której resztę zmierzono w tej samej
sesji jako \(-5{,}07\pm3{,}33\) eV po zmianie wykładnika regularizacji z 6 na
12. Pięć eV wobec \(1450\) eV wiązania to \(0{,}35\%\), czyli rząd
zmierzonego przesunięcia PLC. Poprawka regularizacji zamyka więc sprawę także
tutaj.

*Ile mogłaby dać, gdyby dawała.* Przesunięcie \(\le0{,}23\%\) przerzuci przez
próg tylko trajektorię, która i tak ląduje w obrębie \(0{,}23\%\) od niego,
czyli przy \(\text{PLC}\in[149{,}65,\,150]\). Takie się zdarzają (zmierzone
maksima 149,677 i 149,98), ale to pasmo o szerokości \(0{,}35\) w zakresie
\(122\)–\(150\) populacji dwufotonowej, więc rzędu **1 na 250 trajektorii**.
Dla porównania: zmiana kwantu przesunęła medianę o 45% i liczbę fotonów z 21 na
29, czyli ruszyła stan po drugim fotonie o \(O(1)\), nie o dwie dziesiąte
procenta.

**Wniosek.** Ostrze noża jest własnością **ziarnistości emisji**, a nie kanału
rozpadu. Różnica dipolowa może je przerzucić, ale w skali \(10^{-3}\), i nie
może być mechanizmem stojącym za 45%. Wniosek „różnica kanałów nie wchodzi
przez budżet energetyczny" rozszerza się więc na regułę stopu — tym razem
zmierzony tam, gdzie człon dipolowy faktycznie jest duży, a nie tam, gdzie jest
zaniedbywalny.

#### Warunek początkowy wobec drabiny Bohra: dwa poziomy, jeden punkt zgodności

Podorbitalny start to konwencja SED/Bohra, więc uzgodnienie go z drabiną jest
pytaniem o model, nie o parametr. Rozpisane, pasmo produkcyjne to

```cpp
f_r ∈ [-0.10, +0.10],   f ∈ [0.88, 1.12]      (w jednostkach v_circ, r = a_Ps)
```

a z \(v_{\rm circ}^2=K/r\) wychodzi zamknięcie:

\[\frac{a_0}{a_{Ps}}=\frac{1}{2-f_r^2-f^2},\qquad
  n_E=\sqrt{a_0/a_{Ps}},\qquad
  e_0^2=1+f^2(f_r^2+f^2-2).\]

**Stan niesie więc dwa różne poziomy Bohra naraz**: momentu pędu
\(L/\hbar=f\) i energetyczny \(n_E\). Przyrównanie ich daje
\((f^2-1)^2=0\), czyli zgadzają się przy \(f=1\) i **nigdzie indziej w
paśmie**:

| \(f\) | \(a_0/a_{Ps}\) | \(n_E\) | \(L/\hbar\) | \(e_0\) |
|---|---|---|---|---|
| 0,88 | 0,816 | 0,903 | 0,88 | 0,226 |
| 1,00 | 1,000 | **1,000** | **1,00** | **0** |
| 1,12 | 1,341 | 1,158 | 1,12 | 0,254 |

Poza jednym punktem stan ma \(e_0\) do 0,28 i **nie jest stanem Bohra**. Samo
w sobie nie jest to wadą — to klasyczne orbity. Wadą staje się dopiero
dlatego, że inne części modelu stosują do nich pojęcia drabiny: `quantumFor`
czyta poziom **energetyczny**, `clampAboveGroundStateAngularMomentum` czyta
poziom **momentu pędu**.

*Konsekwencje, zmierzone.*

- \(n_E\) biegnie \(0{,}903\)–\(1{,}166\), a \(P(n_E<1)=0{,}493\) —
  około **połowa trajektorii startuje poniżej stanu podstawowego**. To jest
  dokładnie ten mechanizm, przez który `--ground-state-floor` daje zera.
- Połowa startuje też z \(L=f\hbar<\hbar\), poniżej podłogi momentu pędu tej
  samej flagi, więc pierwszy foton **podniósłby** L — emisja zwiększająca
  moment pędu.
- Pasmo jest centrowane na \(f=1\), ale **nie na \(n=1\)**: \(1/(2-f^2)\)
  jest wypukłe, więc \(E[a_0/a_{Ps}]=1{,}0288\) i \(E[n_E]=1{,}0117\).
  Zespół leży średnio 2,9% na zewnątrz promienia stanu podstawowego.
- Nie jest to drobny efekt. Przy ustalonym \(f\) średni czas kolapsu wynosi
  \(108{,}6\) ps dla \(f=0{,}88\), \(188{,}0\) dla \(f=1{,}00\) i
  \(367{,}2\) dla \(f=1{,}12\) — **czynnik 3,4** w poprzek pasma
  (klasycznie \(t\sim a^3\) przewiduje 4,6). Rozkładając wariancję:
  \(\sigma/\text{średnia}\) to \(1{,}127\) z pasmem wobec \(\approx0{,}93\)
  przy ustalonym \(f\), więc **pasmo daje ok. 32% wariancji** czasu kolapsu, a
  proces fotonowy pozostałe 68%.

*Szerokość pasma jest odziedziczona, nie wyprowadzona.* Komentarz w kodzie
mówi wprost, że przy poprawianiu centrowania zachowano „a comparable spread"
ze starego pasma \([0{,}72,\,0{,}97]\). Wyprowadzenie jej — jako szerokości
równowagi fluktuacyjno-dyssypacyjnej SED — to wątek `--zpf`, a nie stała, którą
da się tu wpisać.

*Ostra alternatywa — WDROŻONA.* \(f=1\), \(f_r=0\) dokładnie. Zweryfikowane:
\(n_0=1\) i \(e_0=0\) **co do cyfry** na każdej trajektorii, a `--level 2`
daje \(n_0=2\) dokładnie — przy \(r=n^2a_{Ps}\) orbita kołowa ma
\(L=n\hbar\) i \(E=-R/n^2\) tożsamościowo, więc `--level n` jest teraz
**ścisłym stanem Bohra**, a nie pasmem wokół niego. Stare pasmo zostaje jedną
flagą dalej (`CREM_INITIAL_BAND=1`) do porównań regresyjnych, tym samym
wzorcem, co poprawka harmoniczna.

Oba losowania są nadal **pobierane**, tylko nieużywane: usunięcie ich
przesunęłoby cały dalszy strumień losowy i zmieniłoby to w hurtową zmianę
trajektorii zamiast zmiany warunku początkowego, niszcząc porównywalność ze
wszystkim, co zmierzono wcześniej.

*Zmierzone skutki,* ta sama binarka, pasmo włączane flagą:

| | pasmo s42 | ostry s42 | pasmo s7 | ostry s7 |
|---|---|---|---|---|
| RMST [ps] | \(218{,}1\pm33{,}1\) | \(193{,}6\pm26{,}4\) | \(219{,}3\pm31{,}2\) | \(191{,}6\pm23{,}3\) |
| \(\sigma/\text{średnia}\) | 1,144 | **0,990** | 1,096 | **0,940** |
| stop na barierze | 19,3% | **47,2%** | 25,9% | **45,6%** |
| wiązanie term. [keV] | 1,584 | **2,515** | 1,847 | **4,602** |

RMST spada o \(\approx12\%\), zgodnie z przewidywaniem, a
\(\sigma/\text{średnia}\) z \(\approx1{,}12\) do \(\approx0{,}96\) — czyli
dokładnie o te \(32\%\) wariancji, które dawało pasmo. Rozrzut czasu kolapsu
jest teraz procesem emisji, a nie odziedziczoną szerokością pasma.

*Najważniejszy skutek jest jednak w trzecim wierszu.* **Bariera zatrzymuje
teraz \(\approx46\%\) trajektorii zamiast \(\approx22\%\)**, a wiązanie
terminalne rośnie z \(1{,}6\)–\(1{,}8\) do \(2{,}5\)–\(4{,}6\) keV. Para
dochodzi do zadeklarowanej granicy fizycznej modelu zamiast wypadać wcześniej
na progu numerycznym — czyli sekcja o regule stopu wyżej opisuje stan, który ta
zmiana w połowie naprawia.

Mechanizm jest ten sam, który kod już notuje przy warunku (b): orbity o
**niższym \(L\)**, czyli bardziej radialne, dosięgają progu
\(T/t_{\rm light}\) przy **większym** periapsis niż kołowe. Stare pasmo dawało
\(e_0\) do \(0{,}25\), nowe daje \(e_0=0\), więc orbity zostają kołowe i
docierają do bariery, zamiast wypaść na retardacji przy \(10\)–\(16\,r^*\).
(To wyjaśnienie, nie osobny pomiar — ale zgadza się i z kierunkiem wiązania
terminalnego, i z udokumentowanym w kodzie mechanizmem.)

*Skutek uboczny, spójny.* Przy ostrym przygotowaniu para startuje **dokładnie**
na \(n=1\), czyli na najniższym stanie, jaki podłoga dopuszcza — więc
`--ground-state-floor` na poziomie 1 nie ma dokąd kaskadować i wszystkie
trajektorie trafiają do strażnika `preparedBelowGroundState` (zmierzone: 15 z
15). To nie jest wada: tak właśnie wygląda stan podstawowy z podłogą emisji.
Sensowną konfiguracją pozostaje `--level 2 --ground-state-floor`. Usunięcie
pasma zamieniło więc **niespójność w warunek brzegowy**.

Warto zauważyć, do czego to prowadzi razem z podłogą: przy ostrym
przygotowaniu para startuje **dokładnie** na \(n=1\), więc pierwszy foton i
tak ją zatrzymuje. Podłoga ma sens wyłącznie ze startem wzbudzonym, i wtedy
kanoniczną konfiguracją jest `--level 2 --ground-state-floor` — kaskada
\(n=2\to1\) jako realny proces. Obie decyzje są więc jedną decyzją.

#### `--level 2 --ground-state-floor` produkcyjnie — i zakleszczenie, które to blokowało

Polecenie brzmiało: wprowadzić tę konfigurację jako domyślną. Pierwszy pomiar
wyglądał na twardy bloker — **0 z 16** trajektorii kończyło przy 60 s budżetu na
trajektorię, dochodząc do \(41{,}8\) ns czasu symulowanego. Ale to nie był
koszt obliczeń, tylko **błąd**.

*Ślad.* Po pierwszym fotonie stan przestawał się ruszać:

```
orbitsToSkip=157093 jumpParameter=0.3 skipHazard=0.153287 ... E_ref=8.52483e-19
orbitsToSkip=157092 jumpParameter=0.3 skipHazard=0.153286 ... E_ref=8.52485e-19
orbitsToSkip=157092 jumpParameter=0.3 skipHazard=0.153286 ... E_ref=8.52486e-19
```

`cumHazard` cyklicznie przekraczał próg i był odejmowany — czyli fotony
**leciały**, co dwa checkpointy — a energia nie drgnęła przez 460 kolejnych
checkpointów. Para stała przy \(n=1{,}367\) i paliła hazard w nieskończoność.

*Przyczyna.* `roomToFloor` przycinał foton do luki **energii wewnętrznej**,
\((E-E_{gs})\mu\). Ale foton o energii \(E_\gamma\) nie obniża energii
wewnętrznej o \(E_\gamma\): para się odrzuca i dokładnie
\(W_b=\sqrt{W_a^2+E_\gamma^2}+E_\gamma\), więc ubytek przewyższa foton o
energię odrzutu \(E_\gamma^2/2W_b\). Przycięcie do luki wewnętrznej sadza
więc parę **tuż poniżej** podłogi — o \(7{,}2\cdot10^{-7}\) względnie, czyli
**700×** więcej niż tolerancja \(10^{-9}\) strażnika poniżej. Strażnik,
którego własny komentarz głosił „this should not trigger", wykonywał się
**zawsze** — a ponieważ robi `break` już po zużyciu hazardu, foton przepadał
bez śladu.

*Poprawka.* Odwrócenie \(W_a^2=W_b^2-2W_bE_\gamma\) względem \(E_\gamma\)
sadzającego \(W_a\) dokładnie na podłodze daje

\[E_\gamma=\frac{W_b^2-W_a^2}{2W_b},\]

co jest mniejsze od luki wewnętrznej dokładnie o odrzut. Tolerancja strażnika
poszła z \(10^{-9}\) na \(10^{-6}\), żeby zwykłe zaokrąglenie przycięcia nie
mogło tego wskrzesić, a fizycznie istotne przestrzelenie nadal go wyzwalało.

*Skutek:*

| | przed | po |
|---|---|---|
| ukończone | **0 / 16** przy 60 s | **14 / 16** przy 45 s |
| mediana KM | — (\(S(t)=1\) po 41,8 ns) | **7161 ps** |
| stop | — | 100% na stanie podstawowym |
| wiązanie terminalne | — | \(6{,}80285\) eV |

*Wprowadzone domyślne:* `gInitialPrincipalLevel = 2`,
`gGroundStateEmissionFloor = true`, budżet wallclock \(20\to90\) s (kaskada
zajmuje \(\approx7\) ns czasu symulowanego wobec \(\approx150\) ps inspiralu
ograniczonego barierą). Standardowy przebieg to teraz **kaskada
\(n=2\to n=1\)** kończąca się na stanie podstawowym. Poprzednie zachowanie
jest pod `--no-ground-state-floor --level 1`.

*Zmiana znaczenia obserwabli, wypisywana teraz przez sam program.* Raportowany
czas **nie jest** klasycznym czasem inspiralu do bariery ani **czasem życia
anihilacji** — model nie ma dynamiki anihilacji, ani kanału kontaktowego, ani
rate'u. Jest czasem kaskady z przygotowanego poziomu do \(n=1\). Zdania w
raporcie, które twierdziły co innego, są teraz warunkowe i mówią, co
faktycznie zachodzi.

*Weryfikacja:* ścieżka bez podłogi nietknięta (zmiana w całości wewnątrz
`if(gGroundStateEmissionFloor)`), ortho zachowuje się identycznie (14/16,
mediana \(7{,}161\) ns), 0 awarii numerycznych, 39/39 PASS.

#### Czy włączyć `--ground-state-floor`? Jeszcze nie — i dlaczego (przebieg rozstrzygnięty wyżej)

Sekcja o regule stopu daje mocny argument ZA: bez podłogi 77% kolapsów kończy
zapas numeryczny, a stan terminalny jest niekontrolowanym przeskokiem. Z
podłogą znika i jedno, i drugie. Zmierzone:

```
stopped by             0 Compton barrier, 0 retardation limit,
                       20 ground-state floor (100%)
terminal binding       0.00680285 keV     = alpha^2 m_e c^2 / 4, co do cyfry
period/light-crossing  861.0 przy stopie, wobec progu 150
```

Model ani razu nie zbliża się do granicy własnej ważności, a stan końcowy jest
deterministyczny i **fizyczny**, zamiast być własnością ziarnistości.

**A mimo to flaga jest dziś zepsuta**, i to cicho. Mediana Kaplana–Meiera
wychodzi **0 ps**, bo:

| ziarno | ukończone | z czasem \(=0\) |
|---|---|---|
| 42 | 20 | **17** |
| 7 | 23 | **23** |

Przyczyna to zderzenie dwóch części modelu. Każda trajektoria startuje na
*promieniu* \(a_{Ps}\), ale z **podorbitalną** prędkością styczną, więc jej
półoś wielka zaczyna poniżej \(a_{Ps}\). Zmierzone: \(a_0/a_{Ps}\) biegnie
\(0{,}821\)–\(1{,}321\) z medianą \(0{,}94\)–\(1{,}02\), czyli poziom
Bohra \(n_0=\sqrt{a_0/a_{Ps}}\) biegnie \(0{,}906\)–\(1{,}149\) i **około
połowa trajektorii zaczyna poniżej \(n=1\)**. Podłoga mówi, że poniżej
\(n=1\) nic nie ma. Spotykają się w \(t=0\), podłoga ogłasza parę osiadłą,
zanim ta ruszy — a taki przebieg wchodził do próby przeżycia jako
zaobserwowany kolaps o długości zero i ściągał estymator do zera.

To jest ta sama granica, co \(93{,}4\%\) checkpointów poniżej \(n=1\) z
sekcji o kwancie emisji, tylko widziana w chwili startu.

*Naprawione.* Para przygotowana na poziomie podłogi albo poniżej jest teraz
oznaczana (`preparedBelowGroundState`), raportowana jako limit obserwacji z
własnym komunikatem i **nie wchodzi do próby przeżycia**. Zamiast fałszywej
mediany 0 ps przebieg mówi, co zaszło, i kieruje do `--level`.

*Konfiguracja, która ma sens.* `--level 2 --ground-state-floor`: \(n_0\)
wychodzi \(2{,}015\), zer nie ma, 100% stopów na stanie podstawowym,
wiązanie terminalne \(6{,}80285\) eV. Kaskada do \(n=1\) jest wtedy realnym
procesem, który model umie zmierzyć. **Koszt: \(\approx2002\) ps na kaskadę
wobec \(\approx100\) ps** dla przebiegu ograniczonego barierą, czyli około
dwudziestokrotnie dłużej — przy budżecie 15 s kończy 1 z 25, a przy
`--level 3` żadna.

*Warunki włączenia,* w tej kolejności: (1) uzgodnić warunek początkowy z
definicją podłogi, żeby domyślny start nie leżał pod nią; (2) zbudżetować
dwudziestokrotne spowolnienie, bo bez tego rozkład jest zdominowany przez
cenzurę. Do tego czasu flaga zostaje eksperymentem — ale nie dlatego, że jest
zła, tylko dlatego, że domyślny warunek początkowy jej przeczy.

#### Podłoga przy \(r^*\): rozważona, zmierzona, odrzucona

Po ustaleniu, że dolny szczebel drabiny to bariera przy
\((2/g)\alpha m_ec^2=3{,}7246\) keV, napraszała się podłoga emisji w tym
miejscu — odpowiednik `--ground-state-floor`, tylko o jedną potęgę \(\alpha\)
niżej. **Zmierzona jako bezczynna** i dlatego niewdrożona.

Rozkład wiązania terminalnego, cztery ziarna po 40 przebiegów:

| ziarno | \(n\) | min | mediana | max | powyżej 3,7246 keV |
|---|---|---|---|---|---|
| 42 | 38 | 0,227 | 1,451 | 3,878 | 1 (2,6%) |
| 7 | 38 | 0,226 | 2,207 | 3,710 | 0 |
| 99 | 37 | 0,226 | 1,588 | 3,658 | 0 |
| 1234 | 38 | 0,226 | 2,095 | 3,875 | 1 (2,6%) |

**2 trajektorie na 151 (1,3%)** w ogóle sięgają progu, i to ledwie. Dla nich
podłoga zapaliłaby się na ostatnim fotonie — w checkpoincie, w którym reguła
stopu i tak kończy przebieg, bo `periapsis <= comptonBarrierRadius` jest
**twardym stopem**. Mechanizm zapalałby się dokładnie wtedy, kiedy przestaje
mieć na co działać.

To jest różnica względem `--ground-state-floor`, i warto ją zapisać: stan
podstawowy leży przy 6,8 eV, czyli **tam, gdzie para startuje**, więc podłoga
ma przed sobą cały kolaps, w którym może działać. Przy \(r^*\) nie ma żadnego
obszaru — to koniec trajektorii, nie jej wnętrze. Dodatkowo 77% przebiegów
nigdy tam nie dociera, bo kończy na progu retardacyjnym przy 10–16 \(r^*\)
(sekcja o regule stopu wyżej).

Jedyny sposób, w jaki ta podłoga mogłaby cokolwiek zmienić, to przycięcie
ostatniego fotonu tak, żeby lądował na \(r^*\) — czyli droga odrzucona w
sekcji o regule stopu, bo wymyśla przejście, którego model nie ma, i ubiera
granicę ważności w fizykę. Odrzucenie stoi.

Drabina \(\alpha\) niżej pozostaje w mocy jako opis **bariery**; to, czego
nie ma, to mechanizm, który by przy niej parkował.

#### Poziom zerowy: drabina w potęgach \(\alpha\), nie na \(m_ec^2\)

Propozycja „poziom zerowy przy promieniu terminalnym **i** stanie energetycznym
co najmniej \(m_ec^2\)" ma dwie połowy, których **nie da się spełnić naraz**.
Promień terminalny nie jest wolnym parametrem: bariera to
\(r^*=g\hbar/4m_ec\), więc wiązanie tam wynosi

\[|E(r^*)|=\frac{K}{2r^*}=\frac{2\alpha m_ec^2}{g},\]

czyli \(3{,}724620\) keV wobec \(\alpha m_ec^2=3{,}728940\) keV — stosunek
\(0{,}998842=2/g\). Te \(0{,}12\%\) to **wyłącznie anomalny moment
magnetyczny elektronu**; poza nim równość jest dokładna. Żeby wiązanie sięgnęło
\(m_ec^2\), para musiałaby zejść do \(1{,}408970\) fm, czyli \(137{,}195\)
razy głębiej niż bariera — w obszar, który bariera istnieje po to, żeby
wykluczyć.

To, o co propozycja zahacza, model **już ma**, tylko szczeble idą w potęgach
\(\alpha\):

| szczebel | promień | energia | postać zamknięta | zmierzony stosunek |
|---|---|---|---|---|
| anihilacja | — | \(510{,}999\) keV | \(m_ec^2\) | — |
| bariera (terminal) | \(193{,}3035\) fm | \(3{,}7246\) keV | \((2/g)\,\alpha m_ec^2\) | \(0{,}998842\) |
| stan podstawowy | \(105{,}835\) pm | \(6{,}8028\) eV | \(\alpha^2m_ec^2/4\) | \(1{,}000000\) |

Dolny szczebel wychodzi **dokładnie**, do wszystkich cyfr, odkąd \(a_{Ps}\)
jest wyprowadzane z momentów magnetycznych pary, a nie wpisywane.

Wniosek dla modelu anihilacji: \(m_ec^2\) **nie jest poziomem orbitalnym**.
Jest tym, co anihilacja uwalnia, a nie tym, do czego para spada. Poziomem
zerowym orbity jest bariera — jedno \(\alpha\) poniżej masy spoczynkowej i
\(4/\alpha\) powyżej stanu podstawowego.

#### Kwant emisji: gdzie jest sztuczny i dlaczego poprawka została wycofana

`ħω_orb` jest prescription **korespondencyjną**, ważną asymptotycznie dla
dużych \(n\), gdzie \(\Delta E(n\to n{-}1)\to2R/n^3=\hbar\omega_{\rm orb}\).
Przy \(n=1\) żąda \(13{,}606\) eV wobec wiązania \(6{,}803\) eV —
**dwukrotnie więcej, niż para ma z czego zapłacić** — a największe rzeczywiste
przejście ze stanu podstawowego to \(1S\)–\(2S=5{,}102\) eV.

*Zasięg, zmierzony:* gałąź różnicy poziomów istniejąca w kodzie odcina się na
\(n\ge2\) i w produkcji **nigdy nie jest wykonywana** (\(0{,}0\%\)
checkpointów; służy przebiegom `--level`). Rozkład:

| zakres | udział | \(\hbar\omega/\Delta E\) |
|---|---|---|
| \(n\ge2\) | \(0{,}0\%\) | — |
| \(1\le n<2\) | \(6{,}6\%\) | mediana \(34\), do \(129\) |
| \(n<1\) | \(\mathbf{93{,}4\%}\) | drabina nie ma stanów |

*Poprawka, którą znam i której nie wdrożyłem.* W oknie \(1\le n<2\) jedynym
stanem końcowym jest \(n=1\), więc kwant to \(E(n)-E(1)=R(1-1/n^2)\) —
ograniczony przez \(R\), a więc zawsze wykonalny, i ciągły z gałęzią wyższą
(obie dają \(0{,}75R=5{,}102\) eV przy \(n=2\)). Naiwne przesunięcie progu
istniejącego wzoru **pogorszyłoby sprawę**: \(R(1/(n-1)^2-1/n^2)\) daje
\(3{,}56R\) już przy \(n=1{,}5\).

**Wycofana, bo skutek jest niewytłumaczony.** Podstawienie przesunęło produkcję
o \(45\%\) (mediana \(119{,}2\to62{,}9\) ps, wiązanie terminalne
\(2{,}667\to3{,}883\) keV) przy zmianie dotykającej \(6{,}6\%\)
checkpointów. Rachunek mówi, że tak być nie może: liczba fotonów skaluje się
jak \(1/\hbar\omega\), a każdy niesie \(\hbar\omega\), więc **energia
usunięta jest niezmiennicza** względem tego podstawienia. Jedno z dwóch jest
fałszywe i nie ustalono które.

Pierwsza hipoteza — że odświeżanie kaskady odbudowuje `hbar*omega_orb` i
rozspójnia hazard z energią — została **obalona pomiarem**: po ujednoliceniu
obu ścieżek wynik wyszedł identyczny co do każdej cyfry, więc ta gałąź nigdy
się nie wykonuje.

*Decyzja.* Cięcie zostaje na \(n\ge2\). Niewytłumaczony czynnik dwa w głównej
obserwabli jest gorszy niż kwant znany jako za duży w oknie obejmującym
\(6{,}6\%\) checkpointów i udokumentowany jako taki. Wycofanie przywraca bazę
w granicach szumu (RMST \(0{,}07\sigma\), średnia \(0{,}02\sigma\)).

*Co zostało naprawione mimo to.* Reguła istniała w **dwóch** kopiach — jawnej i
odbudowywanej w kaskadzie — i to ten sam wzorzec, który tego samego dnia
wyprodukował już jeden błąd produkcyjny (trzy ścieżki aktualizujące elementy, z
których jedna nie miała partnera dla momentu pędu). Obie kopie zastąpiono
jednym `quantumFor`.

*Gdzie naprawdę leży problem.* \(93{,}4\%\) kolapsu przebiega poniżej
\(n=1\), gdzie drabina nie ma stanów, więc i różnicy poziomów do podstawienia.
Sztuczność kwantu jest **objawem** przebywania poza dziedziną drabiny, nie jej
przyczyną — i to samo ograniczenie atakuje z drugiej strony
`--ground-state-floor`, nie poprawiając kwantu, tylko nie wpuszczając tam pary.

#### Progi energetyczne CREM wobec obserwowanych, i skąd bierze się \(a_{Ps}\)

| próg | CREM | obserwacja | status |
|---|---|---|---|
| wiązanie / jonizacja | \(6{,}802847\) eV | \(6{,}8028\) eV | konsekwencja momentu — patrz niżej |
| przejście \(n{=}2\to1\) | \(5{,}102135\) eV | \(5{,}1021\) eV | ta sama drabina |
| kwant emisji \(\hbar\omega_{\rm orb}\) | \(\mathbf{13{,}605693}\) eV | brak takiego progu | **niezgodność** |
| rozszczep nadsubtelny | \(7{,}57\) GHz | \(203{,}394\) GHz | \(3{,}7\%\) |
| skala anihilacyjna | \(186{,}74\) ps | \(125{,}16\) ps (p-Ps) | ten sam \(\alpha^5\), iloraz \(3/2\) |

*Wiersz trzeci jest najostrzejszy.* Kwant korespondencyjny wynosi
\(13{,}6057\) eV — **dwukrotność energii wiązania** i \(2{,}67\times\)
największe rzeczywiste przejście (\(1S\)–\(2S=5{,}10\) eV). Pozytonium nie ma
progu przy \(13{,}6\) eV; para o takiej energii jest już zjonizowana. Jest to
ta sama rzecz, którą model zapisuje jako sufit kinematyczny
(\(\hbar\omega/E_{\rm wiąz}=2\sqrt{a_{Ps}/a}\)), ale zestawiona z obserwacją
mówi ostrzej: **prescription \(\hbar\omega_{\rm orb}\) nie reprezentuje
żadnego rzeczywistego przejścia ze stanu podstawowego.** (Ciekawostka
algebraiczna: \(13{,}605693\) eV to dokładnie stała Rydberga — masa
zredukowana \(m_e/2\) połowi wiązanie, a \(\hbar\omega=2E_{\rm wiąz}\) je
przywraca.)

*Skąd bierze się \(a_{Ps}\).* Model ma **jedną** emergentną długość z własnej
elektrodynamiki: promień zrównania energii magnetycznej z kulombowską,

\[r_*^2=\frac{\mu_0}{4\pi}\frac{\mu_1\mu_2}{K}=\frac{\mu_1\mu_2}{c^2|q_1q_2|},\]

będący **tożsamościowo** tym samym co `comptonBarrierRadius` \(=g\hbar/4mc\),
bo \(\mu=(g/2)(q\hbar/2m)\) daje \(r_*=\mu/(|q|c)\) wprost. Zmierzone:
\(193{,}3035387648\) wobec \(193{,}3035388174\) fm — różnica \(2{,}7\cdot10^{-10}\),
czyli zaokrąglenie CODATA między stabelaryzowanym momentem a \(\hbar\).

Wyeliminowanie \(\hbar\) między \(r_*\) a relacją Bohra daje

\[a_{\rm pary}=\frac{16\,(m_1+m_2)\,\mu_1\mu_2}{g_1g_2\,|q_1q_2|\,K},\]

**bez jawnego \(\hbar\) i bez \(c\)**. `pairBohrRadius` liczy się teraz tą
drogą; zgodność ze starą postacią \(\hbar^2/(\mu_{\rm red}K)\) wynosi
\(10^{-16}\) na czterech parach (e⁺e⁻, µ⁺µ⁻, p⁺p̄, p+e⁻), więc zmiana jest
numerycznie zerowa.

> **Czym to nie jest.** Nie jest klasycznym wyprowadzeniem promienia Bohra:
> \(\mu\) niesie \(\hbar\), więc relacja eliminuje \(\hbar\) między dwiema
> wielkościami, które obie je zawierają. Czym jest: dowodem, że \(a_{\rm pary}\)
> i moment magnetyczny **nie są niezależnymi wejściami**. Raportowane wiązanie
> wynika ze zmierzonego momentu przez klasyczną elektrodynamikę, a nie jest
> drugą stałą kwantową wpisaną obok — i po tej zmianie obie nie mogą się już
> rozjechać, bo są tym samym wyrażeniem.
>
> Wcześniejsze partie tego dokumentu nazywały wiązanie \(6{,}8028\) eV
> „tożsamością, bo \(a_{Ps}\) jest wejściem". Ściślej: jest konsekwencją
> **jednego** wejścia kwantowego, nie dwóch.

#### Błąd: energia zabierana bez momentu pędu

Znaleziony przy okazji badania punktów zwrotnych, i poważniejszy od pytania,
które do niego doprowadziło.

*Objaw.* Przy zamrożonych \((E,L)\) pełny potencjał nie dawał pasma radialnego
dla \(62\%\) checkpointów. Wyglądało to na efekt dipolowy albo usterkę sondy;
było jednym i drugim. Margines energii nad orbitą kołową wynosi
\(\tfrac12K^2e^2/L^2\), a **we wszystkich** przypadkach bez pasma wychodził
dokładnie zero — czyli \(e=0\) co do bitu, co bierze się z klamry
\(\max(0,\cdot)\) w `eccentricitySquared`.

*Przyczyna.* Surowy dyskryminant \(1+2\varepsilon l^2/K^2\) **przed** klamrą:

| | \(n\) | \(\mathrm{disc}\le0\) | mediana | minimum |
|---|---|---|---|---|
| z pasmem | \(466\) | \(0{,}0\%\) | \(+2{,}03\cdot10^{-3}\) | \(+1{,}1\cdot10^{-16}\) |
| bez pasma | \(770\) | \(\mathbf{99{,}0\%}\) | \(-7{,}03\cdot10^{-3}\) | \(-6{,}95\cdot10^{-2}\) |

Ujemny dyskryminant znaczy \(L\) większe od wartości kołowej dla danego
\(E\) — para nie opisuje **żadnej** orbity. Powód: ścieżka produkcyjna
dopisywała stratę energii z mierzonej orbity **nie ruszając momentu pędu**.
Obniżanie \(|E|\) przy stałym \(L\) uokrągla orbitę, a po przekroczeniu
granicy kołowej dyskryminant przechodzi przez zero. Klamra raportowała takie
stany jako orbity kołowe, więc nic tego nie widziało.

Ścieżka fotonowa nigdy tego nie miała (jej komentarz stwierdza spójność „by
construction"), a deterministyczna ścieżka masowa paruje skok energii z
\(L\mathrel{*}=\)`energyGrowth`\(^{\text{angularExponent}}\). Tylko ten
dodatkowy kredyt, dołożony później, poszedł bez partnera.

*Naprawa: oba elementy z elektrodynamiki.* Moment pędu bierze się teraz z
**tej samej kwadratury strumienia dalekiego pola**, z której idzie energia
(`radiatedAngularMomentum` obok `orbitalRadiatedEnergy`), a nie z prawa
Keplera.

Pierwszy kandydat na „pomiar" był pułapką i warto to zapisać:
`deltaAngularMomentumPerOrbit`, różnica przebiegu i tła, jest w trybie
stochastycznym **zerem z konstrukcji** — ten model nie ma ciągłej siły
reakcji, więc przebieg i tło są mechanicznie prawie identyczne. Zmierzone:
\(2{,}5\cdot10^{-12}\) wobec fizycznych \(10^{-4}\), osiem rzędów za mało.
Strumień jest odporny na tę pułapkę, bo nie zależy od tego, czy zastosowano
siłę reakcji — i to jest powód, dla którego strona energetyczna od dawna z
niego korzysta.

*Walidacja krzyżowa.* Strumień wobec prawa \(k(e)=-(1-e^2)/(2+e^2)\), które
stało tu jako pierwsza naprawa: iloraz ma **medianę \(1{,}0000\)** przy
zakresie \(10\)–\(90\) percentyla \(0{,}9923\)–\(1{,}0077\). Prawo było więc
trafne do \(0{,}8\%\), a strumień dokłada retardację, Darwina i dipol, których
ono nie niesie. Prawo zostaje jako fallback przy niedostępnym strumieniu.

*Skutek, zmierzony trójstronnie* (\(200\) zdarzeń, to samo ziarno):

| | A: z błędem | B: prawo \(k(e)\) | C: strumień |
|---|---|---|---|
| średnia ukończonych | \(164{,}04\pm11{,}85\) ps | \(161{,}71\pm11{,}87\) | \(161{,}64\pm11{,}66\) |
| promień terminalny | \(282{,}565\) fm | \(269{,}212\) fm | \(269{,}96\) fm |
| wiązanie terminalne | \(2{,}548\) keV | \(2{,}674\) keV | \(2{,}667\) keV |

\(A\to B\) jest istotne: promień terminalny \(-4{,}7\%\), wiązanie
\(+5{,}0\%\). Błąd zawyżał więc promień końcowy i zaniżał wiązanie o \(5\%\),
co przenosi się na niezmiennik anihilacyjny (\(\approx-126\) eV w \(W\),
\(-63\) eV w linii). \(B\to C\) to \(0{,}28\%\) — różnica \(0{,}8\%\) na
checkpoint **nie kumuluje się**.

*Zastrzeżenie.* Strumień niesie kanał M1, którego `orbitalRadiatedEnergy` nie
niesie. Zmierzone zanieczyszczenie: najwyżej \(9{,}3\cdot10^{-5}\) przy
promieniu terminalnym, \(2{,}2\cdot10^{-15}\) przy \(a_{Ps}\) — rząd poniżej
rozrzutu między obiema drogami, więc świadomie zaakceptowane.

#### Czy warstwa harmonik, okresu i mimośrodu powinna przestać być keplerowska

Poprawka peryapsis usunęła jedną niespójność i odsłoniła drugą: reguła
zatrzymania rozwiązuje pełny potencjał, a mimośród, okres i dobór harmonik
pozostały keplerowskie. Człon \(1/r^3\) łamie twierdzenie Bertranda, więc
orbita precesuje, \(\omega_r\neq\omega_\varphi\), a widmo przenosi się z
całkowitych wielokrotności jednej częstości na kombinacje
\(m\omega_\varphi+n\omega_r\) — czego maszyneria z całkowitym
`harmonicNumber` nie umie wyrazić.

*Zmierzone wzdłuż rzeczywistych trajektorii* (diagnostyka `CREM_APSIDAL`,
\(757\) checkpointów z pełnym zestawem):

| wielkość | mediana | \(90\%\) | \(>1\%\) checkpointów |
|---|---|---|---|
| kąt apsydalny \(|\Delta\varphi/\pi-1|\) | \(2{,}0\cdot10^{-5}\) | \(1{,}1\cdot10^{-4}\) | \(1{,}72\%\) |
| okres \(|T_{\rm pełny}-T_K|/T_K\) | \(0{,}0\) | \(3{,}7\cdot10^{-9}\) | \(1{,}85\%\) |
| mimośród \(|e_{\rm pełny}-e_K|\) | \(1{,}9\cdot10^{-4}\) | \(1{,}9\cdot10^{-3}\) | \(5{,}02\%\) |

Funkcja jest zwalidowana na znanej odpowiedzi: dla dipoli zerowych zwraca
\(1{,}000000000\) przy mimośrodach od \(0{,}8\) do \(0{,}001\). Opis
keplerowski jest więc dokładny do \(10^{-4}\) przez \(90\%\) kolapsu i
załamuje się dopiero w końcowych \(2\)–\(5\%\).

*Dlaczego mimo to nie da się tego tanio naprawić.* Naturalnym pomysłem jest
odmrozić \((E,L)\) i nadać elementom pełny potencjał. Jest to architektonicznie
słuszne — usuwa też artefakt, przez który przy zamrożonych \((E,L)\) i orbicie
prawie kołowej \(56{,}6\%\) checkpointów nie miało zdefiniowanego pasma
radialnego. Blokadą nie są jednak trzydzieści jeden miejsc odczytu elementów,
tylko **analityczny skip**:

```cpp
energyGrowth   = pow(1-jumpParameter, -2.0/3.0);
angularExponent = -(1-e^2)/(2+e^2);
integralFactor = (3/J)*(1-pow(1-J,1.0/3.0));
```

Wykładnik \(-2/3\) pochodzi wprost z prawa \(a^3\), a te formy zamknięte
pokrywają do \(200\,000\) orbit na checkpoint i to one czynią model
obliczalnym. Dla \(-K/r+C/r^3\) nie mają odpowiedników analitycznych.

*Sprawdzone obejście, i dlaczego nie działa.* Gdyby skip degenerował się do
jednej orbity tam, gdzie potencjał przestaje być keplerowski, oba opisy dałoby
się rozdzielić wzdłuż naturalnej granicy. Zmierzone:

| półoś | mediana `orbitsToSkip` | \(|\Delta\varphi/\pi-1|>1\%\) |
|---|---|---|
| \(>10000\) fm | \(9529\) | \(0{,}0\%\) |
| \(1000\)–\(10000\) | \(458\) | \(3{,}2\%\) |
| \(500\)–\(1000\) | \(1079\) | \(100\%\) |
| \(<500\) fm | \(433\) | \(77{,}8\%\) |

Skip **nie** degeneruje się: przy \(a<500\) fm nadal pokrywa \(433\) orbity,
czyli obwiednia keplerowska pracuje pełną parą dokładnie tam, gdzie jest
nieważna. Granicy do rozdzielenia nie ma.

*Rachunek zaburzeń w \(C\) — wykonany, i jego granica.* Dla ustalonej
orientacji dipoli uśredniona po azymucie energia jest dokładnie \(C/r^3\) ze
stałym \(C=-\tfrac{\mu_0}{4\pi}\big[\tfrac12\boldsymbol\mu_1\!\cdot\!\boldsymbol\mu_2-\tfrac32(\boldsymbol\mu_1\!\cdot\!\hat L)(\boldsymbol\mu_2\!\cdot\!\hat L)\big]\),
więc problem jest klasyczną perturbacją Keplera. Przez akcję radialną
\(S_r=\oint p_r\,dr\):

\[\delta S_r=-2C\int_{r_-}^{r_+}\frac{dr}{r^3\sqrt{X}}\]

i kluczowe uproszczenie \(dr/\sqrt X=dt=r^2d\varphi/L\) redukuje to do
\((1/L)\int d\varphi/r\), co dla \(r=p/(1+e\cos\varphi)\) daje \(\pi/(Lp)\) na
pół obiegu. Stąd

\[\delta S_r=-\frac{2\pi CK}{L^3}.\]

**Nie zależy od \(E\)** — a ponieważ \(T_r=\partial S_r/\partial E\), wynika z
tego \(\delta T=0\): poprawka pierwszego rzędu do okresu **znika
tożsamościowo**. To wyjaśnia zmierzoną medianę odchylenia okresu równą dokładnie
\(0{,}0\) i \(90\%\) na \(3{,}7\cdot10^{-9}\). Okres wypada więc z pytania o
niekeplerowskość — nie z powodu małości, tylko struktury.

Z \(\Phi=-\partial S_r/\partial L\) wychodzi natomiast precesja:

\[\frac{\Delta\varphi}{\pi}=1-\frac{3CK}{L^4}.\]

*Zweryfikowane numerycznie:*

| \(a\) [fm] | numerycznie | teoria \(1\). rzędu | iloraz |
|---|---|---|---|
| \(31\,697\) | \(6{,}885\cdot10^{-5}\) | \(6{,}887\cdot10^{-5}\) | \(0{,}9997\) |
| \(4138\) | \(4{,}081\cdot10^{-3}\) | \(4{,}041\cdot10^{-3}\) | \(1{,}0099\) |
| \(1000\) | \(2{,}919\cdot10^{-1}\) | \(6{,}920\cdot10^{-2}\) | \(4{,}22\) |
| \(282{,}6\) | \(-1{,}912\cdot10^{-1}\) | \(+8{,}667\cdot10^{-1}\) | \(-0{,}22\) |

**Wzór jest potwierdzony do \(0{,}03\%\) przy \(31\,697\) fm i do \(1\%\) przy
\(4138\) fm — i załamuje się dokładnie tam, gdzie poprawka staje się duża.**
Przy \(1000\) fm myli się czterokrotnie, przy \(282\) fm ma zły znak, bo
przyciągający \(C/r^3\) wtedy orbitę niszczy, a nie zaburza.

*Wniosek dla obwiedni.* Rachunek zaburzeń działa tam, gdzie \(C\) nie ma
znaczenia, i zawodzi tam, gdzie ma. Obwiednia sekularna psuje się przy małych
\(a\), czyli w obszarze rozbieżności szeregu — perturbacyjna droga nie naprawi
więc skipu w reżimie terminalnym. Nadanie elementom pełnego potencjału
pozostaje możliwe tylko przez numeryczną obwiednię, której skip miał właśnie
unikać. Diagnostyka `CREM_APSIDAL` zostaje w kodzie, żeby każda przyszła próba
miała przeciw czemu się walidować.

*Zmierzony skutek — znacznie mniejszy, niż zapowiadał test wrażliwości.* A/B
sekwencyjne, \(200\) zdarzeń, to samo ziarno:

| | przed | po |
|---|---|---|
| ukończenia | \(186/200\) | \(184/200\) |
| mediana KM | \(120{,}355\) ps | \(123{,}618\) ps |
| RMST | \(181{,}78\pm13{,}24\) ps | \(183{,}73\pm13{,}38\) ps |
| średnia ukończonych | \(163{,}09\pm11{,}74\) ps | \(164{,}04\pm11{,}85\) ps |
| promień terminalny | \(281{,}078\) fm | \(282{,}565\) fm |
| wiązanie terminalne | \(2{,}5615\) keV | \(2{,}54802\) keV |

RMST przesuwa się o \(0{,}15\sigma\), średnia o \(0{,}08\sigma\), ukończenia
mieszczą się w szumie Poissona. Promień terminalny rośnie o \(0{,}53\%\), a
wiązanie spada o tyle samo — spójnie, bo \(r\propto1/E\).

**Test wrażliwości przeceniał skutek**, i warto wiedzieć dlaczego: liczył
przesunięcie przy **zamrożonych** \((E,L)\) dla jednej orbity terminalnej,
podczas gdy w przebiegu poprawiona reguła działa na **każdym** checkpoincie, więc
warstwa sekularna dochodzi do stanu końcowego stopniowo i sama się dostosowuje.
Pozycja \(36{,}4\%\) „brak dozwolonej orbity" nie przekłada się więc na
\(36{,}4\%\) zatrzymanych trajektorii.

#### Naprawa: wykładnik regulatora 6 → 12

Skoro resztka jest proporcjonalna do \(w'(r)\), naprawą jest regulator, którego
gradient jest pomijalny tam, gdzie pomijalna jest już jego wartość. Przy
wykładniku \(6\) waga przy barierze wynosiła \(0{,}998\) — czyli regulator
prawie nie zmieniał *wartości* — ale jego *gradient* dawał \(0{,}79\%\)
potencjału kulombowskiego.

*Wyprowadzenie, nie dobór.* Współczynnik poprzeczny profilu to
\(w[n(1-w)-1]/r^3\), a przy podstawieniu \(u=(r_{\rm reg}/r)^n\) jego ekstremum
spełnia

\[(n-1)(3-n)u^2+(n^2+4n-6)u-3=0.\]

Dla \(n=6\) daje to \(5u^2-18u+1=0\), \(u=(9+2\sqrt{19})/5\), czyli szczyt przy
\((9-2\sqrt{19})^{1/6}\) — **dokładnie stałą, którą kod już nosił**, co
weryfikuje wyprowadzenie na znanym przypadku przed użyciem go do nowego. Dla
\(n=12\): \(33u^2-62u+1=0\), \(u=(31+4\sqrt{58})/33\), szczyt przy
\([33/(31+4\sqrt{58})]^{1/12}=0{,}94949269\).

| | wykł. 6 | wykł. 12 |
|---|---|---|
| \(K_n\) (szczyt \(|U|r_{\rm reg}^3/(\mu_0\mu_1\mu_2/4\pi)\)) | \(1{,}5244265\) | \(2{,}7783644\) |
| \(r_{\rm reg}\) | \(68{,}472\) fm | \(83{,}639\) fm |
| \(w\)(bariera) | \(0{,}998028\) | \(0{,}999957\) |
| sufit energii dipolowej | \(0{,}5\) | \(0{,}5\) |

Analityczne \(K_{12}\) zgadza się z niezależną maksymalizacją numeryczną
(\(4\cdot10^6\) punktów) do \(1{,}19\cdot10^{-13}\), a energia w analitycznym
szczycie trafia w sufit \(0{,}5\,mc^2\) z odchyłką \(1{,}7\cdot10^{-16}\).

*Efekt, zmierzony.* Resztka przy barierze, \(N=10^6\) orientacji:

```
wykladnik 6:   -58,9 eV   przy 7-9 sigma
wykladnik 12:   -5,07 eV  +- 3,33   ->  1,52 sigma
```

czyli z jednoznacznie wykrywalnej staje się **nieodróżnialną od zera**, zgodnie
z przewidywaniem analitycznym \(-2{,}57\) eV (rozbieżność \(0{,}75\sigma\)).

*Koszt dla produkcji: żaden.* A/B sekwencyjne na pustej maszynie, \(200\)
zdarzeń, to samo ziarno — wyniki **identyczne co do cyfry**: mediana KM
\(120{,}355\) ps, RMST \(181{,}78\pm13{,}2435\) ps, średnia ukończonych
\(163{,}086\pm11{,}7395\) ps, promień terminalny \(281{,}078\) fm, ukończenia
\(186/200\) w obu. Trajektorie nie docierają tam, gdzie regulator się różni, a
przy barierze nowy regulator jest **mniej** natrętny także co do wartości.

> Walidator złapał przy tej zmianie dwie rzeczy. `static_assert` przypinający
> \(r_{\rm reg}\) wywalił kompilację, bo zapomniałem go przeliczyć. Potem
> `measuredDipoleCurlPeak` rozminął się o \(3{,}0\cdot10^{-4}\) przy tolerancji
> \(10^{-4}\) — ale winna była **siatka, nie stała**: szczyt przy wykładniku
> \(12\) jest węższy w skali logarytmicznej, a \(4000\) punktów na sześć dekad
> daje \(0{,}35\%\) odstępu w promieniu. Rozstrzygnął to test przy
> analitycznym promieniu (odchyłka \(1{,}7\cdot10^{-16}\)), więc siatka
> została zagęszczona do \(16\,000\) punktów zamiast rozluźnienia tolerancji —
> ta ostatnia broni przed cichym powrotem produkcji do skalarnego surogatu
> \(w/r^3\) i osłabianie jej byłoby naprawianiem miernika zamiast pomiaru.

Wniosek dla ewentualnego modelu anihilacji trzeba więc osłabić: różnica kanałów
w budżecie energetycznym **istnieje przy barierze**, ale jest o \(1{,}6\%\)
Coulomba i wymaga \(\sim10^4\) trajektorii, żeby ją zobaczyć. W QED różnica
kanałów jest regułą wyboru na stanie końcowym i daje trzy rzędy w czasie życia
— czego \(1{,}6\%\) w energii wiązania nie odtworzy. Model odtwarza
*strukturę* (ścisłe kasowanie M1 w jednym kanale), ale nie *wielkość*.

*Granica całego kierunku.* Model z trzema podłogami odtwarza stan podstawowy
co do cyfry — ale odtwarza go, bo podano mu **trzy rzeczy**: \(a_{Ps}\) przez
podłogę energetyczną, \(\hbar\) przez podłogę na momencie pędu i
\(S\in\{0,1\}\) przez kwantowanie spinu. Zgodność promienia terminalnego z
\(a_{Ps}\) jest **tautologiczna**. To, co model realnie wnosi, to droga do
tego stanu i kinematyka wyjścia z niego — nie sam stan.

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

**Ten sam audyt znalazł też rozbieżność w precesji spinu** — poniżej opisaną
w czasie, w jakim ją wtedy widziano, jako nierozstrzygniętą i nienaprawioną
celowo (wymagającą osobnej derywacji, nie zgadywania, która strona jest
błędna). *Została rozstrzygnięta i naprawiona* — pełny przebieg dochodzenia,
łącznie z dwiema błędnymi naprawami po drodze, jest w akapitach
„**Naprawione**" i „**Naprawa właściwa**" niżej w tej samej sekcji. Zdania w
czasie teraźniejszym w trzech kolejnych akapitach należy więc czytać jako
zapis stanu ówczesnego, nie bieżącego. Kod ma dwie implementacje: `thomasBmtEffectiveField`
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
relatywistycznych prędkościach było to wtedy otwarte pytanie. *(Zamknięte —
patrz „Naprawa właściwa" niżej: `advanceCovariantBmt` nie jest już ścieżką
produkcyjną i nosi `[[maybe_unused]]`, a `properDipolePrecessionRate` stosuje
wzór Jacksona wprost do momentu spoczynkowego. Bramkowany test
`bmt-precession-invariant` pilnuje tego od tamtej pory.)*

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

Naprawa (pierwsze podejście, **też błędna** — opis i korekta niżej): nowa
`advanceThomasBmtDipole`/`properDipolePrecessionRate`
([modules/electrodynamics.hpp](modules/electrodynamics.hpp)) omijała
czterowektor całkowicie. Zamiast tego: (1) rozwiń dipol laboratoryjny
(`lorentzBoostDipole`) z bieżącego `properDipole`; (2) policz jego tempo ze
wzoru Jacksona zastosowanego do `μ_lab`; (3) odwzoruj tempo z powrotem na
`properDipole` odwrotnym boostem w zamkniętej formie.

**Korekta drugiego audytu: krok (2) stosował równanie Jacksona do
niewłaściwego wektora.** Jackson (3. wyd., 11.170) definiuje `s` jako spin
**w układzie spoczynkowym**; `state.firstProperDipole` *jest* tym momentem
spoczynkowym ([modules/state.hpp](modules/state.hpp)), a `state.firstDipole`
to dopiero jego boost. Równanie należało więc zastosować wprost do
`properDipole`, bez żadnego boostu w żadną stronę. Złożenie
`L⁻¹[L(μ)×B_eff]` nie jest rotacją `μ_proper`: zmierzone
`μ_proper·dμ_proper/dt` było niezerowe na poziomie względnym `γ−1` —
`2,3·10⁻⁵` przy β pozytonium, `4,3·10⁻³` przy β=0,1, `0,63` przy β=0,9 —
**niezależnie od `g`**, więc w odróżnieniu od błędu czterowektorowego
**nie znikało przy g=2**.

Dryf normy, jaki to wytwarzało, był pochłaniany co krok przez renormalizację
w `advanceThomasBmtDipole` i opisany tam jako obcięcie RK4 `O(dt⁵)`. To było
błędne rozpoznanie: po usunięciu renormalizacji dryf **na jednostkę czasu
jest stały**, `7,16·10⁸ 1/s`, płaski dla `dt` od `10⁻¹⁸` do `10⁻²²` (pięć
dekad), podczas gdy człon `O(dt⁵)` spadałby `10⁴×` na dekadę. Był to dryf
pierwszego rzędu, czyli człon fizycznej wielkości, a nie błąd numeryczny.
Przy jednym `dt` obie rzeczy wyglądają identycznie — rozróżnia je dopiero
skan po `dt`, którego sam próg dryfu normy wykonać nie mógł.

Naprawa właściwa: `properDipolePrecessionRate` to teraz dosłownie

\[
\frac{d\boldsymbol\mu_{\rm proper}}{dt}
=\frac{q}{m}\,\boldsymbol\mu_{\rm proper}\times\mathbf B_{\rm BMT},
\]

czyli iloczyn wektorowy z `thomasBmtEffectiveField` i nic więcej — bez
czterowektora, bez boostu, bez projekcji. Jest to precesja z samej
konstrukcji, więc zachowuje `|μ_proper|` dokładnie, tak jak norma momentu
spoczynkowego zachowywać się musi. Renormalizacja w `advanceThomasBmtDipole`
zostaje, ale teraz jest już wyłącznie ubezpieczeniem od rzeczywistego
obcięcia RK4: bez niej dryf normy schodzi do `1,1·10⁻¹⁶` przy `dt=10⁻¹⁸` i do
zera maszynowego niżej.

**Test certyfikujący poprzednią „naprawę" był tautologią.** `BMT vs eff
field` porównywał `d(state.firstDipole)/dt` z `μ_lab×B_eff` — dokładnie tę
tożsamość, którą ówczesna implementacja była zbudowana spełniać. Zwracał
zero maszynowe z konstrukcji i nie mógł zwrócić niczego innego; spadek
`5,8·10⁻⁶ → 2,2·10⁻¹⁶` nie był potwierdzeniem poprawki. Sprawdzenie, którego
podmiot i odniesienie są tym samym równaniem, nie mierzy niczego.

Zastąpiony parą sond w [modules/maxwell_validation.hpp](modules/maxwell_validation.hpp):
`BMT vs eff field` zostaje jako jawnie oznaczone sprawdzenie okablowania
(łapie podmienione `g`, `q/m`, pole albo zepsute RK4 — nie wybór równania),
a nowy, **bramkowany** test `bmt-precession-invariant` mierzy
`|μ·dμ/dt|/(|μ||dμ/dt|)`, czyli ortogonalność tempa do momentu. To
niezmiennik precesji, niezależny od tego, jaki wzór jest podpięty, więc
zachowuje moc rozróżniającą po dowolnym przepisaniu sektora: dla trasy przez
dipol laboratoryjny czytałby `4,3·10⁻³` (β=0,1) i `0,63` (β=0,9), dla
poprawnej precesji czyta zero maszynowe.

Zweryfikowane po korekcie: `positronium_validation` **38/38** (`bmt-precession-invariant`
dodany do sekcji „algebraic identities"), `BMT precession: 0 / 0`,
`BMT vs eff field: 0 / 0`. Niezależne potwierdzenie z testu, którego nie
dotykano: residuum `covariant BMT` — mierzące kowariancję Lorentza ewolucji
dipola przez `advanceCovariantBmt`, czyli zupełnie inne sformułowanie
(czterowektorowe), zasilane polami z historii **ewoluowanych produkcyjnie** —
spadło z `2,7568153·10⁻⁷` do `7,7877557·10⁻¹¹`, czyli 3540×. Gdyby zmiana
szła w złą stronę, ta liczba by wzrosła. `--diagnose` (e⁺e⁻, seed 12345)
daje `trajectory: PASS`, a trajektoria jest praktycznie nietknięta —
zmiany na poziomie `10⁻⁹` względnie, zgodnie z oczekiwaniem: przy β≈0,0037
pary e⁺e⁻ czynnik `γ−1≈6,7·10⁻⁶`, a kanał dipolowy jest i tak podrzędny.
**Żaden opublikowany wynik e⁺e⁻ się nie zmienia**; poprawka ma znaczenie dla
par relatywistycznych i wysoko-anomalnych (`--pair` z protonami, sondy β=0,9)
oraz blisko kolapsu. `advanceCovariantBmt` pozostaje w kodzie (oznaczona
`[[maybe_unused]]`), podobnie jak `inverseTensorBoostMagnetic`, której
produkcja już nie używa — obie służą testom samospójności pod boostem w
`maxwell_validation.hpp`.

**Domknięcie: całkowanie zastąpione rozwiązaniem ścisłym.** Po naprawie
zostało pytanie, skąd bierze się resztkowy dryf normy `μ_proper` na pełnej
trajektorii (`max |mu| drift`, mierzony w `--diagnose` względem stałej
`(g/2)·magneton`). Odpowiedź: nie z fizyki i nie z obcięcia RK4, lecz z samej
**renormalizacji**. `advanceThomasBmtDipole` kończył się przeskalowaniem
`result*(targetNorm/result.norm())`, czyli pierwiastkiem i dzieleniem w każdym
podkroku; ich zaokrąglenia kumulują się przez ~10⁶ podkroków trajektorii
szybciej niż zaokrąglenia samego obrotu.

Ponieważ `applyDipolePrecession` zamraża prędkość i pole na czas całego
wywołania, `B_eff` jest w podkroku **wektorem stałym**, więc równanie to
\(\dot{\boldsymbol\mu}=\boldsymbol\omega\times\boldsymbol\mu\) ze stałym
\(\boldsymbol\omega\), którego rozwiązaniem ścisłym jest sztywny obrót o kąt
\(|\boldsymbol\omega|\,dt\). `advanceThomasBmtDipole` liczy teraz ten obrót
wzorem Rodriguesa zamiast całkować go RK4. Zmierzone:

| kąt kroku | RK4 + renormalizacja | wzór Rodriguesa |
|---|---|---|
| 10⁻⁵ rad (skala produkcyjna) | 1,49·10⁻¹³ | 1,49·10⁻¹³ |
| 0,1 rad | 7,92·10⁻⁸ | 4,25·10⁻¹² |
| 1 rad | 5,60·10⁻³ | 4,13·10⁻¹³ |
| 2 rad | 1,00·10⁻¹ | 1,57·10⁻¹² |

Przy kątach, które produkcja faktycznie wykonuje, obie metody są nierozróżnialne
(obie na poziomie zaokrąglenia) — dawny komentarz twierdzący, że RK4 jest tu
„dokładne niezależnie od długości kroku", był więc praktycznie prawdziwy, choć
formalnie nie. Różnica pojawia się dopiero przy dużych kątach, gdzie RK4 myli
się o 10%, a obrót pozostaje ścisły; to zapas bezpieczeństwa, nie bieżąca
poprawka. Realne zyski są trzy:

- `max |mu| drift` na trajektorii produkcyjnej: `2,55·10⁻¹²%` → `1,73·10⁻¹²%`
  (1,47×). Skromnie, bo w produkcji renormalizacja nie jest jedynym źródłem
  zaokrągleń; w długim teście syntetycznym (2·10⁶ kroków) różnica to 15×
  (`4,24·10⁻¹¹` → `2,79·10⁻¹²`).
- Koszt: **2,1× szybciej** (2·10⁷ wywołań: 4,79 s → 2,25 s) — jedno wyliczenie
  `thomasBmtEffectiveField` zamiast czterech, bez pierwiastka i dzielenia
  renormalizacji.
- Najważniejsze, strukturalnie: **norma jest zachowana z konstrukcji, a nie
  naprawiana.** Renormalizacja była dokładnie tym, co ukryło pierwotny błąd —
  pochłaniała co krok człon pierwszego rzędu. Bez niej każdy przyszły defekt
  strukturalny w tym sektorze musi ujawnić się jako dryf, zamiast zostać po
  cichu usunięty.

`properDipolePrecessionRate` zostaje (`[[maybe_unused]]`) jako **specyfikacja**
tego, co obrót ma rozwiązywać, i to ona jest podmiotem testów. Dzięki temu
`BMT vs eff field` przestał być tautologią: porównuje teraz rozwiązanie
zamknięte z numerycznym scałkowaniem specyfikacji przez 0,1 rad (4096
podkroków RK4), co wychwyci błędną oś obrotu, przekręcony znak albo zgubiony
czynnik — czego wariant różnicowy pierwszego rzędu wychwycić nie mógł, bo brał
krok tak krótki, że obie strony zgadzały się trywialnie. Zmierzone
`2,36·10⁻¹⁵`. Dołożony `BMT norm (2 rad)`: 64 obroty po 2 rad, bez żadnej
renormalizacji, dryf `1,33·10⁻¹⁵`. Oba wchodzą do bramki
`bmt-precession-invariant`.

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
