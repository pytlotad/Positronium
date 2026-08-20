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

\[a=\sqrt[3]{\frac{\mu_0}{4\pi}\frac{\mu_1\mu_2}{E}}=47{,}22\,\mathrm{fm},\]

czyli separacją, przy której klasyczna energia oddziaływania dipol-dipol osiąga
\(E\) — energię spoczynkową **lżejszego** składnika pary. Wchodzą oba momenty,
\(\mu_1\mu_2\), a nie kwadrat jednego: to ta sama liczba wyłącznie dla cząstki
i jej antycząstki, bo moment protonu jest 658 razy mniejszy od elektronowego,
więc p+e⁻ potrzebuje 5,43 fm tam, gdzie e⁺e⁻ potrzebuje 47,22. Promień jest
wyprowadzany z wybranej pary, a nie wpisany jako stała.

Momenty są przy tym tymi, które niesie dynamika, czyli \((g/2)\mu_B\) dla
elektronu, a nie gołym magnetonem Bohra. Wzór z \(\mu_B^2\) i masą elektronu
dawał 47,18 fm i był poprawny wyłącznie dla pozytonium. Profil \(w(r)/r^3=r^3/(r^6+a^6)\) ma
maksimum dokładnie w \(r=a\), więc ten wybór ogranicza klasyczną energię
dipolową do \(m_ec^2/2\) w całej dziedzinie. Model klasyczny nie powinien
generować energii oddziaływania powyżej progu kreacji par.

Poprzednia wartość \(0{,}5\cdot\)`nuclearCutoff` \(=5\,\mathrm{fm}\) była
dobrana wyłącznie tak, aby wygładzanie pozostało wewnątrz raportowanej
dziedziny. Dopuszczała energię dipolową \(2{,}1\cdot10^8\,\mathrm{eV}\),
czyli 420 mas spoczynkowych elektronu, i czyniła sektor magnetyczny 106 razy
bardziej punktowym niż sektor ładunkowy, który już zakłada źródło o skończonym
promieniu `chargeCloudRestRadius`. Podniesienie promienia do 47 fm pozostawia
każdą raportowaną orbitę bit-identyczną: przy dziesiątkach pikometrów waga
różni się od jedności o mniej niż \(10^{-19}\); sprawdzono to na trybie
diagnostycznym oraz na eksperymencie 3.

To **nie** usuwa krótkozasięgowej bariery dipolowej. Bariera leży tam, gdzie
nieregularyzowane energie dipolowa i kulombowska się zrównują,
\(r^\*=\sqrt{(\mu_0/4\pi)\mu_1\mu_2/(ke^2)}=193\,\mathrm{fm}\), czyli w
połowie zredukowanej długości Comptona — daleko poza jakimkolwiek promieniem
wygładzania, który nie rusza raportowanych orbit. Poniżej tej skali klasyczny
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
- reszta skończonego strumienia promieniowania po boostcie \(0{,}2084\)
  i skumulowanego czteropędu promieniowania \(0{,}2084\) — najsłabszy wynik
  kowariancyjny całego audytu i najostrzejsze ograniczenie interpretacji
  strumienia w układach ruchomych. Scharakteryzowany pomiarowo, choć przyczyna
  nie jest jeszcze wskazana:

  - **nie jest to zbyt krótka historia opóźniona.** Promień powierzchni
    kontrolnej \(10^6a_0\) wymaga czasu opóźnionego \(1{,}8\cdot10^{-13}\) s,
    a historia sondy obejmuje \(7{,}2\cdot10^{-19}\) s, czyli 245 000 razy
    mniej. Zjechanie z promieniem do \(1\,a_0\), gdzie czas opóźniony mieści
    się w historii z czterokrotnym zapasem, zbija resztę zaledwie z
    \(0{,}2084\) do \(0{,}1695\);
  - **nie jest to rozdzielczość kwadratury.** Zmiana liczby kierunków
    Lebiediewa z 50 na 770 zmienia resztę dopiero na piątej cyfrze;
  - **jest to efekt geometrii boostu, pierwszego rzędu w \(\beta\).** Reszta
    wynosi \(0{,}2084\), \(0{,}0917\), \(0{,}0414\) i \(0{,}0202\) przy
    \(\beta=0{,}35\), \(0{,}20\), \(0{,}10\) i \(0{,}05\); iloraz
    \(\text{reszta}/\beta\) dąży do \(0{,}40\), a lokalny wykładnik do
    \(1{,}04\). Znika więc liniowo z prędkością układu.

  Dalej wykluczone: **nie są to źle transformowane przyspieszenia.**
  Przyspieszenie własne, będące niezmiennikiem Lorentza, zgadza się między
  układami do \(2{,}6\cdot10^{-8}\) przy \(\beta=0{,}35\). Sektor dynamiczny
  jest kowariantny; niekowariantne jest samo liczenie strumienia.

  Rozbicie residuum wskazuje winowajcę. Składowa **energii** jest kwadratowa
  w \(\beta\) (\(1{,}96\cdot10^{-1}\), \(1{,}42\cdot10^{-2}\),
  \(3{,}52\cdot10^{-3}\) przy \(\beta=0{,}35\), \(0{,}10\), \(0{,}05\)),
  czyli zachowuje się poprawnie. Dominuje składowa **pędu**, liniowa w
  \(\beta\): wypromieniowany pęd w układzie ruchomym jest konsekwentnie
  \(1{,}40\) raza większy, niż przewiduje boost, a czynnik ten dąży do stałej
  przy \(\beta\to0\) (\(1{,}334\), \(1{,}394\), \(1{,}398\)).

  Silna poszlaka co do mechanizmu: dla promieniowania dipolowego
  \(dP/d\Omega\propto\sin^2\theta\) zachodzi
  \(\langle\cos^2\theta\rangle=\tfrac15\), a \(1+2\langle\cos^2\theta\rangle
  =1{,}4\) — zgodność z \(1{,}398\) na poziomie \(0{,}15\%\). Sugeruje to, że
  reguła użyta do porównania, „zboostuj tempo i podziel przez \(\gamma\)",
  pomija zależny od kształtu wiązki człon aberracji: tempo przez sferę przy
  ustalonym czasie **nie jest** czterowektorem, bo \(dt\) nie jest
  niezmiennikiem. Nie jest jeszcze rozstrzygnięte, czy wada leży w regule
  porównania w teście, czy w samym integratorze strumienia — obie dają ten sam
  podpis.

  Wniosek: wypromieniowany czteropęd **nie jest kowariantny w pierwszym rzędzie
  w \(\beta\)** i wada siedzi w składowej pędowej, nie energetycznej. Ponieważ produkcja liczy w jednym układzie, nie
  unieważnia to bezpośrednio raportowanych czasów kolapsu, ale zabrania
  interpretowania strumienia jako wielkości kowariantnej;
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

Bilans energii jest **tożsamością algebraiczną**, nie testem zachowania:
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
rekurencyjnego obliczania pola BMT.

Dostępny jest alternatywny model `coherentElectricDipole`: wspólne pole reakcji
Abrahama–Lorentza jest wyznaczane z tej samej trzeciej pochodnej elektrycznego
momentu dipolowego pary. Siły na elektron i pozyton są dokładnie przeciwne;
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
w(r)=\frac{1}{1+(r_m/r)^6},\qquad r_m=0{,}5r_{\rm cutoff}.
\]

Przejście regulatora leży wewnątrz wyłączonego obszaru modelu punktowego,
więc nie deformuje raportowanej części trajektorii na skali atomowej.
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
liczonym ze wszystkich użytecznych trajektorii.

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
dwukrotne niedopromieniowanie. Podana wartość pochodzi z przebiegu N=1000
o pełnej kompletacji (1000/1000 trajektorii dotarło do granicy, zero cenzury),
przy rozrzucie sigma/średnia = 0,41.

Zmierzony czas życia p-Ps to 125,1 ps, więc pozostała rozbieżność to czynnik
\(\approx3{,}7\). Świadomie nie próbujemy jej domknąć: silnik odtwarza teraz
zamknięty wzór klasycznej inspirali z dokładnością 2%
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
oddziaływania. Każda z \(N\) prób ma wagę \(A/N\). Wyświetlane są:

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
kątowym. Zdarzenie dochodzące do \(10^{-14}\,\mathrm m\) trafia do osobnego
kanału `reach cutoff`; jego przekrój nie jest przekrojem anihilacji QED.
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
| `N_1_2_annihilation_time.pdf` | Widmo czasu anihilacji z **mierzonej** stałej rozpadu, rysowane analitycznie. Pionowa linia = średni czas kolapsu CREM, dla porównania skal. |
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
> samych miarach czyste, na poziomie \(10^{-4}\) energii orbity.
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
| `5_1_2_collision_energy.pdf` | Losowana energia w układzie środka masy, z zaznaczeniem podpróbki związanej. |
| `5_1_3_impact_parameter.pdf` | Losowany parametr zderzenia (rozkład półnormalny), z podpróbką związaną. |
| `5_1_4_dipole_alignment.pdf` | Wyrównanie dipoli w stanach związanych, z progiem para/orto przy `+0,5`. Stosunek 1:3 wynika tu z **geometrii** progu na izotropowej sferze, a nie ze statystyki spinowej — zbieżność liczbowa jest przypadkowa. |
| `5_1_5_annihilation_time_para.pdf` | Widmo czasu anihilacji z mierzonej stałej rozpadu p-Ps, zawężone do zdarzeń sklasyfikowanych jako para. |
| `5_1_6_annihilation_time_ortho.pdf` | To samo dla klasy orto. |
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
| `particle_species.hpp` | tablica gatunków (masa, ładunek, `g`), para integrowana przez przebieg i skale, które z niej wynikają: masa zredukowana, promień Bohra pary, energia wiązania, promień regularyzacji dipola oraz wyszukiwanie gatunku dla `--pair` |
| `electrodynamics.hpp` | **prawa sił**: pola opóźnione, Darwin, sprzężenie dipolowe, reakcja promieniowania, strumień dalekiego pola |
| `crem_engine.hpp` | **numeryka**: rekonstrukcja historii przyczynowej i adaptacyjny integrator (sonda błędu, podział kroku, retencja historii) |
| `crem_trajectory.hpp` | warstwa trajektorii: `Frame`, wspólna pętla całkowania, sampler warunków początkowych i klasyfikacja zjawiska |
| `crem_collapse.hpp` | estymator kolapsu: całkowanie sekularne z uśrednianiem po orbitach oraz zamknięte odniesienia elektrodynamiczne |
| `statistics_archive.hpp` | katalog wartości zmierzonych i teoretycznych |
| `root_export.hpp` | atomowy zapis PDF |
| `maxwell_validation*.hpp` | zestaw testów budowany do `positronium_validation` |

Trzy moduły `crem_*` nie zawierają **żadnego kodu ROOT** — cała prezentacja
pozostaje w `positronium.cpp`. Podział przeprowadzono jako przeniesienie
tekstu: wyniki po nim są bit-identyczne (walidacja 28/28, diagnostyka i czas
kolapsu bez zmian).

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
wiązki e⁺e⁻. Tryb, wybór eksperymentu i ziarno można podać bez interakcji:

```bash
./positronium --mode visual --visual-style line --phenomenon 2 --seed 42
./positronium --mode visual --visual-style dot --phenomenon 2 --seed 42
./positronium --mode statistical --phenomenon 1 --runs 1000 --seed 42
./positronium --mode statistical --phenomenon 4 --runs 1000 --seed 42 \
    --beam-energy-ev 20 --theta-min-deg 5 --angle-bins 10
```

Parametry wiązki można dodatkowo kontrolować przez `--bmax-pm` i
`--matching-radius-pm`; gdy ich nie podano, program dobiera je z energii i
akceptancji kątowej.

Cztery opcje sterują samą fizyką i kosztem eksperymentów związanych:

| Opcja | Domyślnie | Znaczenie |
| --- | --- | --- |
| `--zpf`, `--zpf-band` | `0` (wyłączone) | **Eksperyment, nie część modelu.** Klasyczne pole punktu zerowego elektrodynamiki stochastycznej: losowe fale płaskie o widmie \(\rho(\omega)=\hbar\omega^3/2\pi^2c^3\), 64 mody o równej energii, orientacje i fazy z ziarna `--seed`. `--zpf` skaluje **amplitudę** (1 = poziom fizyczny, moc pochłaniana rośnie jak kwadrat), `--zpf-band lo,hi` ustala pasmo w jednostkach częstości orbitalnej pary (domyślnie `0.3,3`). To jest fluktuacyjna połowa pary fluktuacja–dyssypacja; dyssypacyjną, czyli reakcję promieniowania, model ma od zawsze. Wchodzi w te same trzy miejsca co pole jednorodne, ale próbkowane osobno dla każdej cząstki, bo zależy od położenia i czasu. **Nie odtwarza stanu podstawowego SED — patrz niżej.** |
| `--external-field` | brak (pytanie na starcie) | Jednorodne zewnętrzne pole magnetyczne w mikroteslach; `0` wyłącza. Orientacja jest losowana izotropowo z ziarna `--seed`, więc odtwarza się razem z resztą przebiegu, i jest wypisywana na starcie. Gdy opcji nie podano, a przebieg jest interaktywny, program pyta o to **przed wszystkimi pozostałymi pytaniami** i oferuje 50 µT (skala pola ziemskiego). Przebieg wsadowy z podanym `--mode` i `--phenomenon` nigdy nie pyta i domyślnie nie ma pola. Pole wchodzi w sumę sił chwilowych, w sumę sił retardowanych oraz w pole lokalne widziane przez obie cząstki, przez co obejmuje precesję Thomasa-BMT. Przy 50 µT tempo cyklotronowe \(eB/m\) wynosi 8,8·10⁶ rad/s wobec tempa orbitalnego rzędu 3·10¹⁵ rad/s, więc orbita pozostaje nietknięta, a widocznym kanałem jest precesja dipoli — około 3·10⁻⁴ rad w ciągu 35 ps kolapsu. |
| `--pair` | `electron,positron` | Para cząstek, którą całkuje przebieg, podana jako `pierwsza,druga`. Dostępne gatunki: `electron`, `positron`, `muon`, `antimuon`, `proton`, `antiproton`. Para musi być przyciągająca i nieść przeciwne ładunki elementarne, inaczej opcja jest odrzucana. Wybrana para jest wypisywana na starcie wraz z masą zredukowaną, promieniem Bohra pary i energią wiązania. Honoruje ją także `./positronium_validation`. |
| `--radiation-reaction` | `individual` | Model reakcji promieniowania ładunku: `disabled`, `coherent` (Abraham-Lorentz na dipolu elektrycznym pary), `individual` (Landau-Lifszyc zredukowanego rzędu, osobno dla każdej cząstki) albo `automatic` (mieszanka obu). Przy `disabled` żaden kanał nie odbiera energii orbitalnej, więc klasyczna inspirala nie zachodzi i eksperymenty 1/2 zgłaszają brak zaniku. Wybrany model jest wypisywany na starcie. |

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
panele `annihilation_time_*` z danych pozytonium, co dla innej pary jest mylące.

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
| `--crem-wallclock-budget-s` | `20` | Budżet zegarowy na jedno zdarzenie w eksperymentach 1/2. Po jego wyczerpaniu trajektoria jest cenzurowana prawostronnie. Wartość domyślna daje niską frakcję ukończenia przy `N = 1000`; do wniosków ilościowych trzeba ją podnieść tak, aby kompletacja zbliżyła się do 100%. |

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

Literatura uzasadnia użyte elementy elektrodynamiki klasycznej i wejścia
generatora zaniku, ale nie waliduje fenomenologicznego utożsamienia orientacji
klasycznych dipoli ze stanami para-/ortopozytonium w animacji.
