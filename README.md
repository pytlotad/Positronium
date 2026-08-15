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

Wszystkie obliczenia wykonywane są w jednostkach SI. Elektron i pozyton mają tę
samą masę \(m_e\), przeciwne ładunki \(-e\) i \(+e\) oraz momenty magnetyczne o
stałej wartości magnetonu Bohra \(\mu_B\). Początkowa odległość wynosi jeden
promień Bohra \(a_0\), a środek masy układu spoczywa.

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
trwałego momentu \(|\boldsymbol\mu|=\mu_B\) przy klasycznym \(g=1\) jest
założeniem modelu, nie klasycznym wyprowadzeniem własności elektronu. Rzeczywisty
moment elektronu, spin 1/2, splątanie, wymiana, energia wiązania, anihilacja i
poprawki radiacyjne należą do QED i nie mogą zostać dokładnie odtworzone przez
ten integrator klasyczny.

Bieżący `positronium_validation` potwierdza własności algebraiczne i numeryczne,
ale jego końcowe `PASS` oznacza wyłącznie przejście ustawionych progów regresji.
Nie jest certyfikatem kompletności fizycznej. W audycie z 14 sierpnia 2026 r.
uzyskano między innymi:

- resztę kowariancji pojedynczego pola Liénarda–Wiecherta
  \(9{,}89\cdot10^{-16}\), lecz resztę skończonego strumienia promieniowania
  po boostcie \(0{,}2105\) i skumulowanego czteropędu promieniowania
  \(0{,}2109\);
- rozbieżność skumulowanej pracy indywidualnej reakcji LL względem
  koherentnego strumienia ładunkowego równą \(0{,}4884\) energii
  wypromieniowanej;
- resztę siły po boostcie \(1{,}10\cdot10^{-6}\), ale resztę ewolucji tensora
  dipolowego \(3{,}86\cdot10^{-3}\);
- dokładne numeryczne domknięcie raportowanego bilansu cząstki–promieniowanie–
  pole związane, ponieważ `boundField*` jest obliczane jako brakująca reszta
  bilansu po każdym kroku. Jest to księgowanie, a nie niezależne obliczenie
  energii i pędu pola bliskiego. To domknięcie jest **tożsamością algebraiczną**:
  wykresy `diagnostic_*_balance` i pozycje `identity resid` / `identity |dP|` /
  `identity |dJ|` w trybie diagnostycznym mierzą wyłącznie zaokrąglenia i błąd
  interpolacji punktu końcowego, więc nie są testem zachowania. Niezależnymi
  miarami są `|E_bound|/E_rad` (na orbitach związanych 1,2–1,4) oraz
  `|dE_LL-vs-flux|/E_rad` (0,03–0,09); pierwsza z nich wyznacza faktyczną
  granicę wiarygodności raportowanej energii wypromieniowanej.

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

Model zachowuje klasyczną wartość orbitalną \(g=1\), czyli \(a=-1/2\), bez
kwantowej anomalii momentu elektronu. Uwzględnia precesję Thomasa, wpływ pola
elektrycznego i podłużnej składowej pola magnetycznego. Znaki precesji
elektronu i pozytonu są przeciwne. Obrót jest wykonywany wzorem Rodriguesa w
dwóch symetrycznych półkrokach i zachowuje \(|\boldsymbol\mu_i|=\mu_B\).

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
2. `Statistical analysis` wykonuje domyślnie 100 kalibracji CREM dla kanałów 1/2
   oraz szybki podgląd 20/100 trajektorii odpowiednio dla kanałów 3/4. Zestaw
   paneli zależy od wybranego eksperymentu; program nie wymusza tych samych
   czterech histogramów dla zjawisk o innej fizyce.

### Związane pozytonium

Dla przygotowanego p-Ps lub o-Ps eksperymentalną charakterystyką rozpadu jest
\(\Gamma=1/\tau\), a nie przekrój czynny. Obecny CREM nie zawiera kwantowego
operatora anihilacji. Statistical 1/2 oblicza zamiast tego operacyjny klasyczny
czas kolapsu: dla każdej realizacji krótki odcinek wielu orbit jest całkowany
pełnym modelem CREM, a zmierzona średnia moc promieniowania kalibruje równanie

\[
\frac{dE}{dt}=P_{\rm CREM},\qquad
E=-\frac{k e^2}{2a},\qquad P(a)\propto a^{-4}.
\]

Równanie sekularne jest całkowane do `0.01*a0`. Wynik nie jest kwantowym czasem
anihilacji. Zewnętrzne pomiary nie są wejściem CREM; tworzą wyłącznie krzywe
porównawcze

\[
t_{\rm collapse}=\frac{|E_0|}{3\langle P_{\rm CREM}\rangle}
\left[1-\left(\frac{0.01a_0}{a}\right)^3\right].
\]

Jest to przybliżenie uśrednione po orbicie: pełna trajektoria kalibracyjna trwa
`1.5 fs`, natomiast dalsze tysiące orbit nie są całkowane krok po kroku.
Domyślne 100 trajektorii zajęło w pomiarze referencyjnym około 2,6 minuty na
eksperyment przy czterech wątkach; czas zależy od procesora.

\[
P_{\rm exp}(t)=\tau_{\rm exp}^{-1}e^{-t/\tau_{\rm exp}},
\]

z parametrami pochodzącymi z cytowanych pomiarów. Osobny generator produktów
zachowuje energię i pęd w układzie spoczynkowym pozytonium.

Dla para-pozytonium wyświetlane są:

- rozkład ekstrapolowanego czasu kolapsu CREM i krzywa eksperymentalna;
- idealna linia energii dwóch fotonów około \(511\,\mathrm{keV}\);
- rozkład \(\cos\theta_\gamma\) dla niepolaryzowanego źródła;
- karta kinematyczna kanału \(2\gamma\), w którym fotony są przeciwbieżne.

Dla orto-pozytonium wyświetlane są:

- rozkład ekstrapolowanego czasu kolapsu CREM i krzywa eksperymentalna;
- inkluzywne widmo energii fotonów w przybliżeniu Ore’a–Powella;
- dwuwymiarowy wykres Dalitza \(E_{max}\)–\(E_{mid}\);
- kąt między dwoma najbardziej energetycznymi fotonami kanału \(3\gamma\).

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
czterowymiarowym równaniem Thomasa–BMT dla klasycznego \(g=1\). Krok RK4
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
distributions/2_1_3_three_photon_dalitz.pdf
distributions/3_1_3_energy_loss_cross_section.pdf
distributions/4_1_1_differential_cross_section.pdf
```

Diagnostyki drugiej strony również są zapisywane jako osobne pliki. Każdy PDF
zawiera dokładnie jeden pad i jeden odpowiadający mu rozkład; plik nazwany
`diagnostic_energy_closure` nie zawiera pozostałych trzech diagnostyk. Zbiorcze
PDF całych stron Statistical nie są tworzone. Kolejne uruchomienie atomowo
nadpisuje pliki o tych samych nazwach; program najpierw sprawdza poprawność
nowego PDF, aby nie utracić poprzedniego obrazu przy błędzie renderowania.
Eksport działa także z `--no-gui`.

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
./positronium --mode statistical --phenomenon 1 --runs 100 --seed 42
./positronium --mode statistical --phenomenon 4 --runs 100 --seed 42 \
    --beam-energy-ev 20 --theta-min-deg 5 --angle-bins 10
```

Parametry wiązki można dodatkowo kontrolować przez `--bmax-pm` i
`--matching-radius-pm`; gdy ich nie podano, program dobiera je z energii i
akceptancji kątowej. `--no-gui` wykonuje serię i wypisuje podsumowanie bez
otwierania okna ROOT, ale nadal renderuje i zapisuje pliki PDF, co jest
przydatne w obliczeniach wsadowych:

```bash
./positronium --mode statistical --phenomenon 2 --runs 100 \
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
  Statistical 1/2 raportuje odrębny, operacyjny czas kolapsu CREM otrzymany
  przez kalibrację pełnej trajektorii i ekstrapolację sekularną, natomiast
  zewnętrzne czasy życia służą wyłącznie do porównania;
- odpowiedź detektora, oddziaływania pozytonium w materiale, pick-off,
  quenching, poprawki radiacyjne i rzadsze kanały rozpadu;
- relatywistyczne amplitudy QED niezwiązanej pary, w tym rozpraszanie Bhabhy
  i przekrój anihilacji e⁺e⁻;
- pełne, wzajemnie opóźnione pola Liénarda–Wiecherta i kompletna dynamika pola
  Maxwella;
- poprawki oddziaływania ładunków wyższe niż rząd \(v^2/c^2\);
- relatywistyczna precesja BMT, precesja Thomasa, transformacja dipoli między
  układami odniesienia oraz kwantowe sprzężenie spin–orbita;
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
